subroutine matrix_rotate(thetax, thetay, thetaz, R)
    implicit none
    real, intent(in) :: thetax, thetay, thetaz
    real, intent(out), DIMENSION(:,:), allocatable ::  R
    allocate(R(3, 3))
    R(1,1) = cos(thetay)*cos(thetaz)
    R(1,2) = cos(thetay)*sin(thetaz)
    R(1,3) = -sin(thetay)
    R(2,1) = sin(thetax)*sin(thetay)*cos(thetaz) - cos(thetax)*sin(thetaz)
    R(2,2) = sin(thetax)*sin(thetay)*sin(thetaz) + cos(thetax)*cos(thetaz)
    R(2,3) = sin(thetax)*cos(thetay)
    R(3,1) = cos(thetax)*sin(thetay)*cos(thetaz) + sin(thetax)*sin(thetaz)
    R(3,2) = cos(thetax)*sin(thetay)*sin(thetaz) - sin(thetax)*cos(thetaz)
    R(3,3) = cos(thetax)*cos(thetay)
    
end subroutine matrix_rotate

subroutine matrix_multiply(A, B, OUT)
    ! This subroutine multiplies the matrices A and B.
    ! 
    ! INPUT:
    !   A(M,N)  --  A 2D matrix of 64 bit floats.
    !   B(N,P)  --  A 2D matrix of 64 bit floats,
    ! 
    ! OUTPUT:
    !   OUT(M,P)  --  The matrix that is the result of (AB).
    ! 
    USE ISO_FORTRAN_ENV, ONLY: REAL64 ! <- Get a float64 type.
    IMPLICIT NONE  ! <- Make undefined variable usage raise errors.
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:) :: A, B
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(A,1),SIZE(B,2)) :: OUT
  
    ! Compute the matrix multiplication of A and B.
    OUT(:,:) = MATMUL(A,B)
  
END SUBROUTINE matrix_multiply

function sensor_line_generator(x_boundary, Nxf) result(pos_vector)
    integer, intent(in) :: x_boundary
    integer, intent(in) :: Nxf
    real, dimension(3, Nxf) :: pos_vector
    real, dimension(Nxf) :: x
    integer :: i
    write(*,*)"NXf", Nxf
    do i = 1, Nxf
        x(i) = 2.0 * x_boundary * (i - 1) / (Nxf - 1) - x_boundary
    end do
    pos_vector(1, :) = x
    pos_vector(2, :) = 0.0
    pos_vector(3, :) = 0.0
end function sensor_line_generator

subroutine trimmer_index(pos_vectors, tol, first_and_last)
  implicit none
    real, intent(in) :: pos_vectors(:, :)
    real, intent(in) :: tol
    integer, intent(out), dimension(2) :: first_and_last
    real, dimension(:), allocatable :: factor
    integer :: n, i, first_index, last_index
    logical, allocatable :: mask(:)
    real :: eddy_range
    ! Calculate the factor
    eddy_range = 3.0
    factor = pos_vectors(1, :)**2 + pos_vectors(2, :)**2 + pos_vectors(3, :)**2
    n = size(factor)
    ! Create a logical mask for the condition
    allocate(mask(n))
    do i = 1, n
      mask(i) = sqrt(factor(i)) < eddy_range
    end do
    !write(*,*) mask
    ! Initialize indices
    first_index = 0
    last_index = 0
    ! Find the first and last indices that satisfy the condition
    do i = 1, n
      if (mask(i) .and. first_index == 0) then
        first_index = i
      end if
    end do
    last_index = first_index + nint(2*eddy_range/tol)
    first_and_last(1) = first_index
    first_and_last(2) = last_index
    !write(*,*) factor(first_index:last_index)
end subroutine trimmer_index

subroutine trimmer(first_and_last, xaxis, xaxis_trimmed)
  implicit none
  integer, intent(in) :: first_and_last(:)
  real, intent(in), allocatable :: xaxis(:, :)
  real, intent(out), dimension(:, :), allocatable :: xaxis_trimmed
  allocate(xaxis_trimmed(3, first_and_last(2) - first_and_last(1) + 1))
  xaxis_trimmed = xaxis(:, first_and_last(1):first_and_last(2))
end subroutine trimmer


subroutine velocity_calc(xaxis_trimmed, u, v, w)
  implicit none
  real, intent(in) :: xaxis_trimmed(:, :)
  real, intent(out), allocatable :: u(:), v(:), w(:)
  real, dimension(:), allocatable :: factor
  integer :: N
  N = size(xaxis_trimmed, 2)
  allocate(factor(N), u(N), v(N), w(N))
  factor = xaxis_trimmed(1, :)**2 + xaxis_trimmed(2, :)**2 + xaxis_trimmed(3, :)**2
  u = -xaxis_trimmed(2, :) * exp(-factor * 2.0)
  v = xaxis_trimmed(1, :) * exp(-factor * 2.0)
  w = 0.0d0
end subroutine velocity_calc

subroutine vector_sums(u, v, w, u_total, v_total, w_total, first_index, last_index)
  implicit none
  real, intent(in) :: u(:), v(:), w(:)
  real, intent(out), dimension(:) :: u_total, v_total, w_total
  integer, intent(in) :: first_index, last_index
  u_total(first_index:last_index) = u_total(first_index:last_index) + u
  v_total(first_index:last_index) = v_total(first_index:last_index) + v
  w_total(first_index:last_index) = w_total(first_index:last_index) + w
end subroutine vector_sums

subroutine main_calculation(input_ints, a_list, theta_list, velocity_total)
  implicit none
  integer, intent(in), dimension(:) :: input_ints
  real, intent(in), dimension(:,:) :: a_list, theta_list
  real, intent(out), dimension(:,:), allocatable :: velocity_total
  real, dimension(:,:), allocatable :: pos_vector_translated, pos_vector
  real, dimension(:), allocatable :: u, v, w, u_total_0, v_total_0, w_total_0, u_total, v_total, w_total
  real, dimension(:,:), allocatable :: R, R_inv, xaxis_trimmed, xaxis_rotated, velocities, a_list_T, velocity_total_trans, velocity_pre
  integer, dimension(2) :: first_and_last
  integer :: i , j, N, x_boundary, Nxf, N_E, first_index, last_index
  real :: tol

  interface
    subroutine velocity_calc(xaxis_trimmed, u, v, w)
      implicit none
      real, intent(in) :: xaxis_trimmed(:, :)
      real, intent(out), allocatable :: u(:), v(:), w(:)
    end subroutine velocity_calc

    function sensor_line_generator(x_boundary, Nxf) result(pos_vector)
      integer, intent(in) :: x_boundary
      integer, intent(in) :: Nxf
      real, dimension(3, Nxf) :: pos_vector
    end function sensor_line_generator
  end interface

  interface
    subroutine matrix_rotate(thetax, thetay, thetaz, R)
      implicit none
      real, intent(in) :: thetax, thetay, thetaz
      real, intent(out), dimension(:,:), allocatable :: R
    end subroutine matrix_rotate
  end interface

  interface
    subroutine vector_sums(u, v, w, u_total, v_total, w_total, first_index, last_index)
      implicit none
      real, intent(in) :: u(:), v(:), w(:)
      real, intent(out), dimension(:) :: u_total, v_total, w_total
      integer, intent(in) :: first_index, last_index
    end subroutine vector_sums
  end interface

  interface
    subroutine trimmer_index(pos_vectors, tol, first_and_last)
      implicit none
      real, intent(in) :: pos_vectors(:, :)
      real, intent(in) :: tol
      integer, intent(out), dimension(2) :: first_and_last
    end subroutine trimmer_index
  end interface

  interface
    subroutine trimmer(first_and_last, xaxis, xaxis_trimmed)
      implicit none
      integer, intent(in) :: first_and_last(:)
      real, intent(in), allocatable :: xaxis(:, :)
      real, intent(out), dimension(:, :), allocatable :: xaxis_trimmed
    end subroutine trimmer
  end interface
  write(*,*) "Input Variables", input_ints
  x_boundary = input_ints(1)
  Nxf = input_ints(2)
  N_E = input_ints(3)
  tol = 0.05
  write(*,*) "Maximum x value in a_list", maxval(a_list(1, :))
  pos_vector = sensor_line_generator(x_boundary, Nxf)
  write(*,*) "Sensor Line Generated"
  allocate(pos_vector_translated(3, Nxf))
  allocate(u_total_0(Nxf), v_total_0(Nxf), w_total_0(Nxf))
  do i = 1, Nxf
    u_total_0(i) = 0.0
    v_total_0(i) = 0.0
    w_total_0(i) = 0.0
  end do
  u_total = u_total_0
  v_total = v_total_0
  w_total = w_total_0
  write(*,*) "Initial Variables Allocated"
  a_list_T = transpose(a_list)
  write(*,*) shape(a_list_T)
  write(*,*) shape(a_list)
  do i = 1, N_E
!    write(*,*) a_list(:, i)
    if (mod(i,1000) == 0) then
      write(*,*) "Eddy Loop", i
    end if
    do j = 1, Nxf
        pos_vector_translated(:, j) = pos_vector(:,j) - a_list(:, i)
    end do
    call trimmer_index(pos_vector_translated, tol, first_and_last)
    call trimmer(first_and_last, pos_vector_translated, xaxis_trimmed)
    N = first_and_last(2) - first_and_last(1) + 1
    first_index = first_and_last(1)
    last_index = first_and_last(2)
    if (first_index == 0) then
      cycle
    end if
    call matrix_rotate(theta_list(i, 1), theta_list(i, 2), theta_list(i, 3), R)
    R_inv = transpose(R)
    xaxis_rotated = matmul(R, xaxis_trimmed)
    call velocity_calc(xaxis_rotated, u, v, w)
    allocate(velocity_pre(3, N))
    velocity_pre(1, :) = u
    velocity_pre(2, :) = v
    velocity_pre(3, :) = w
    velocities = matmul(R_inv, velocity_pre)
    deallocate(velocity_pre)
    u_total(first_index:last_index) = u_total(first_index:last_index) + velocities(1, :)
    v_total(first_index:last_index) = v_total(first_index:last_index) + velocities(2, :)
    w_total(first_index:last_index) = w_total(first_index:last_index) + velocities(3, :)
  end do
  write(*,*) "Iteration Complete"
  allocate(velocity_total_trans(3, Nxf))
  velocity_total_trans(1,:) = u_total
  velocity_total_trans(2,:) = v_total
  velocity_total_trans(3,:) = w_total
  velocity_total = transpose(velocity_total_trans)
  write(*,*) "Velocity output shape", shape(velocity_total)
  write(*,*) "Component shapes:", shape(velocity_total(:,1)), shape(velocity_total(:,2)), shape(velocity_total(:,3))
end subroutine main_calculation