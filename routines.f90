
subroutine matrix_rotate_f(thetax, thetay, thetaz, R)
    implicit none
    real, intent(in) :: thetax, thetay, thetaz
    real, intent(out), DIMENSION(3,3) ::  R
    R(1,1) = cos(thetay)*cos(thetaz)
    R(1,2) = cos(thetay)*sin(thetaz)
    R(1,3) = -sin(thetay)
    R(2,1) = sin(thetax)*sin(thetay)*cos(thetaz) - cos(thetax)*sin(thetaz)
    R(2,2) = sin(thetax)*sin(thetay)*sin(thetaz) + cos(thetax)*cos(thetaz)
    R(2,3) = sin(thetax)*cos(thetay)
    R(3,1) = cos(thetax)*sin(thetay)*cos(thetaz) + sin(thetax)*sin(thetaz)
    R(3,2) = cos(thetax)*sin(thetay)*sin(thetaz) - sin(thetax)*cos(thetaz)
    R(3,3) = cos(thetax)*cos(thetay)
    
end subroutine matrix_rotate_f

subroutine matrix_multiply_f(A, B, OUT)
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
  
END SUBROUTINE matrix_multiply_f

function sensor_line_generator_f() result(pos_vector)
    use types
    implicit none
    type (input_data) :: id
    real :: x_boundary
    integer :: Nxf
    real, dimension(3, Nxf) :: pos_vector
    real, dimension(Nxf) :: x
    Nxf = id%Nxf
    x_boundary = id%x_boundary

    integer :: i
    x = [(2.0d0 * x_boundary * (i - 1) / (Nxf - 1) - x_boundary, i = 1, Nxf)]
    pos_vector(1, :) = x
    pos_vector(2, :) = 0.0d0
    pos_vector(3, :) = 0.0d0
    write(*,*) pos_vector
end function sensor_line_generator_f

subroutine eddy_range_f(pos_vectors, first_index, last_index)
    real(kind=8), intent(in) :: pos_vectors(:, :)
    integer, intent(out) :: first_index, last_index
    real(kind=8), allocatable :: factor(:)
  
    integer :: n, i
    logical, allocatable :: mask(:)
  
    ! Calculate the factor
    n = size(pos_vectors, 2)
    allocate(mask(n))
    factor = pos_vectors(1, :)**2 + pos_vectors(2, :)**2 + pos_vectors(3, :)**2
  
    ! Create a logical mask for the condition
    mask = sqrt(factor) < 4.0d0
  
    ! Initialize indices
    first_index = 0
    last_index = 0
  
    ! Find the first and last indices that satisfy the condition
    do i = 1, n
      if (mask(i) .and. first_index == 0) then
        first_index = i
      end if
      if (mask(n + 1 - i) .and. last_index == 0) then
        last_index = n + 1 - i
      end if
      if (first_index /= 0 .and. last_index /= 0) exit
    end do
  
    ! Deallocate the logical array
    deallocate(mask)
  end subroutine eddy_range_f

subroutine velocity_generator_f(xaxis_trimmed, factor, u, v, w)
    real(kind=8), intent(in) :: xaxis_trimmed(:, :), factor(:)
    real(kind=8), intent(out) :: u(:), v(:), w(:)

    u = -xaxis_trimmed(2, :) * exp(-factor * 2.0d0)
    v = xaxis_trimmed(1, :) * exp(-factor * 2.0d0)
    w = 0.0d0
end subroutine velocity_generator_f


subroutine total_velocities_f(Nx, x_boundary, N_E, theta_list, a_list, u_total, v_total, w_total)
    integer, intent(in) :: Nx, N_E
    real(kind=8), intent(in) :: x_boundary
    real(kind=8), intent(in) :: theta_list(3, N_E)
    real(kind=8), intent(in) :: a_list(3, N_E)
    real(kind=8), intent(out) :: u_total(:), v_total(:), w_total(:)

    interface
        subroutine eddy_range_f(pos_vectors, first_index, last_index, factor)
            real(kind=8), intent(in) :: pos_vectors(:, :)
            integer, intent(out) :: first_index, last_index
            real(kind=8), intent(out) :: factor(:)
        end subroutine eddy_range_f
    end interface

    interface
        subroutine velocity_generator_f(xaxis_trimmed, factor, u, v, w)
            real(kind=8), intent(in) :: xaxis_trimmed(:, :), factor(:)
            real(kind=8), intent(out) :: u(:), v(:), w(:)
        end subroutine velocity_generator_f
    end interface

    real(kind=8), dimension(Nx, 3) :: xaxis, xaxis_translated
    real(kind=8), dimension(:, :), allocatable :: eddy_pos_translated, eddy_pos_trimmed, eddy_pos_rotated
    real(kind=8), dimension(3, 3) :: R, R_inv
    real(kind=8), dimension(:), allocatable :: factor, u, v, w
    integer :: i, first_index, last_index

    xaxis = sensor_line(x_boundary, Nx)
    u_total = 0.0d0
    v_total = 0.0d0
    w_total = 0.0d0

    do i = 1, N_E
        xaxis_translated = xaxis - spread(a_list(i, :), dim=1, ncopies=Nx)
        call eddy_range_f(xaxis_translated, first_index, last_index, factor)
        allocate(u(last_index - first_index + 1), v(last_index - first_index + 1), w(last_index - first_index + 1))

        eddy_pos_trimmed = (xaxis(:, first_index:last_index))
        call rotation_total_f(theta_list(i, 1), theta_list(i, 2), theta_list(i, 3), R)
        R_inv = transpose(R)

        eddy_pos_rotated = matmul(R, eddy_pos_trimmed)
        call velocity_generator_f(eddy_pos_rotated, factor, u, v, w)

        u_total(first_index:last_index) = u_total(first_index:last_index) + matmul(R_inv(:, 1), reshape(u, [1, size(u)]))
        v_total(first_index:last_index) = v_total(first_index:last_index) + matmul(R_inv(:, 2), reshape(v, [1, size(v)]))
        w_total(first_index:last_index) = w_total(first_index:last_index) + matmul(R_inv(:, 3), reshape(w, [1, size(w)]))
    end do
end subroutine total_velocities_f



