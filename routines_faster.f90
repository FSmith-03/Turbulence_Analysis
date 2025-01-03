
subroutine matrix_rotate(thetax, thetay, thetaz, R)
  implicit none
  real, intent(in) :: thetax, thetay, thetaz
  real, intent(out), DIMENSION(:,:), allocatable ::  R
  if (.not. allocated(R)) then
    allocate(R(3, 3))
  end if
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
  eddy_range = 3.5
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


subroutine velocity_calc(xaxis_trimmed, ARP, ARN, L_e, u, v, w)
implicit none
real, intent(in) :: xaxis_trimmed(:, :)
integer, intent(in) :: ARP, ARN, L_e
real, intent(out), allocatable :: u(:), v(:), w(:)
real, dimension(:), allocatable :: factor
integer :: N
N = size(xaxis_trimmed, 2)
allocate(factor(N), u(N), v(N), w(N))
factor = -(((xaxis_trimmed(1, :)**2 + xaxis_trimmed(2, :)**2) * 2.0) / (ARP * L_e)**2 + (xaxis_trimmed(3, :)**2 / (ARN * L_e)**2))
u = -xaxis_trimmed(2, :) * exp(factor)
v = xaxis_trimmed(1, :) * exp(factor)
w = 0.0d0
end subroutine velocity_calc

subroutine main_calculation(input_ints, a_list, theta_list)
  implicit none

  ! Input parameters
  integer, intent(in), dimension(:) :: input_ints
  real, intent(in), dimension(:,:) :: a_list, theta_list

  ! Output parameter
  real, dimension(:,:), allocatable :: velocity_total

  ! Local variables
  real, dimension(:,:), allocatable :: pos_vector_translated, pos_vector
  real, dimension(:), allocatable :: u, v, w, u_total_0, v_total_0, w_total_0, u_total, v_total, w_total
  real, dimension(:,:), allocatable :: R, R_inv, xaxis_trimmed, xaxis_rotated, velocities, a_list_T, velocity_total_trans, velocity_pre
  integer, dimension(2) :: first_and_last
  integer :: i, j, N, x_boundary, Nxf, N_E, first_index, last_index
  integer :: ARP, ARN, L_e
  real :: tol

  ! Interface declarations
  interface
    subroutine velocity_calc(xaxis_trimmed, ARP, ARN, L_e, u, v, w)
      implicit none
      real, intent(in) :: xaxis_trimmed(:, :)
      integer, intent(in) :: ARP, ARN, L_e
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

  ! Print input variables
  write(*,*) "Input Variables", input_ints

  ! Initialize variables
  x_boundary = input_ints(1)
  Nxf = input_ints(2)
  N_E = input_ints(3)
  ARP = 1
  ARN = 1
  L_e = 1
  tol = 0.05

  write(*,*) "Maximum x value in a_list", maxval(a_list(1, :))

  ! Generate sensor line
  pos_vector = sensor_line_generator(x_boundary, Nxf)
  write(*,*) "Sensor Line Generated"

  ! Allocate arrays
  allocate(pos_vector_translated(3, Nxf))
  allocate(u_total_0(Nxf), v_total_0(Nxf), w_total_0(Nxf))

  ! Initialize total velocity arrays
  u_total_0 = 0.0
  v_total_0 = 0.0
  w_total_0 = 0.0

  u_total = u_total_0
  v_total = v_total_0
  w_total = w_total_0

  write(*,*) "Initial Variables Allocated"

  ! Transpose a_list
  a_list_T = transpose(a_list)
  write(*,*) shape(a_list_T)
  write(*,*) shape(a_list)

  ! Loop over eddies
  do i = 1, N_E
    if (mod(i, 1000) == 0) then
      write(*,*) "Progress", real(i) / real(N_E) * 100.0, "%"
    end if

    ! Translate position vectors
    do j = 1, Nxf
      pos_vector_translated(:, j) = pos_vector(:, j) - a_list(:, i)
    end do

    ! Trim position vectors
    call trimmer_index(pos_vector_translated, tol, first_and_last)
    call trimmer(first_and_last, pos_vector_translated, xaxis_trimmed)

    N = first_and_last(2) - first_and_last(1) + 1
    first_index = first_and_last(1)
    last_index = first_and_last(2)

    if (first_index == 0) then
      cycle
    end if

    ! Rotate position vectors
    call matrix_rotate(theta_list(i, 1), theta_list(i, 2), theta_list(i, 3), R)
    R_inv = transpose(R)
    xaxis_rotated = matmul(R, xaxis_trimmed)

    ! Calculate velocities
    call velocity_calc(xaxis_rotated, ARP, ARN, L_e, u, v, w)

    ! Allocate and populate velocity_pre array
    if (.not. allocated(velocity_pre)) then
      allocate(velocity_pre(3, N))
    end if
    velocity_pre(1, :) = u
    velocity_pre(2, :) = v
    velocity_pre(3, :) = w

    ! Rotate velocities back
    velocities = matmul(R_inv, velocity_pre)

    ! Update total velocities
    u_total(first_index:last_index) = u_total(first_index:last_index) + velocities(1, :)
    v_total(first_index:last_index) = v_total(first_index:last_index) + velocities(2, :)
    w_total(first_index:last_index) = w_total(first_index:last_index) + velocities(3, :)
    
    OPEN(UNIT=10, FILE="debug.dat", ACTION='WRITE')
    WRITE(10, *) i, first_index, last_index, N
    write(10, *) velocity_pre
    CLOSE(10)
  end do

  write(*,*) "Iteration Complete"

  ! Allocate and populate velocity_total_trans array
  if (.not. allocated(velocity_total_trans)) then
    allocate(velocity_total_trans(3, Nxf))
  end if
  velocity_total_trans(1, :) = u_total
  velocity_total_trans(2, :) = v_total
  velocity_total_trans(3, :) = w_total

  ! Transpose velocity_total_trans to get velocity_total
  velocity_total = transpose(velocity_total_trans)

  ! Print output shapes and average velocities
  write(*,*) "Velocity output shape", shape(velocity_total)
  write(*,*) "Component shapes:", shape(velocity_total(:, 1)), shape(velocity_total(:, 2)), shape(velocity_total(:, 3))
  write(*,*) "Average Velocity", sum(velocity_total(:, 1)**2) / Nxf, sum(velocity_total(:, 2)**2) / Nxf, sum(velocity_total(:, 3)**2) / Nxf

  ! Deallocate arrays
  deallocate(pos_vector_translated, pos_vector, u, v, w, u_total_0, v_total_0, w_total_0, u_total, v_total, w_total, R, R_inv, xaxis_trimmed, xaxis_rotated, velocities, a_list_T, velocity_total_trans)

  write(*,*) "Output file name: output.dat"
  OPEN(UNIT=10, FILE="output.dat", STATUS='REPLACE')
  WRITE(10, *) velocity_total
  CLOSE(10)
end subroutine main_calculation
