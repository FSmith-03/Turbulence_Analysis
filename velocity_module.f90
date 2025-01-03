module velocity_module
    implicit none
    contains
  
    function sensor_line(x_boundary, Nxf) result(pos_vector)
      real(kind=8), intent(in) :: x_boundary
      integer, intent(in) :: Nxf
      real(kind=8), dimension(3, Nxf) :: pos_vector
      real(kind=8), dimension(Nxf) :: x
  
      integer :: i
      x = [(2.0d0 * x_boundary * (i - 1) / (Nxf - 1) - x_boundary, i = 1, Nxf)]
      pos_vector(1, :) = x
      pos_vector(2, :) = 0.0d0
      pos_vector(3, :) = 0.0d0
    end function sensor_line
  
    subroutine eddy_range(pos_vectors, first_index, last_index, factor)
      real(kind=8), intent(in) :: pos_vectors(:, :)
      integer, intent(out) :: first_index, last_index
      real(kind=8), intent(out) :: factor(:)
      integer :: n
      real(kind=8), dimension(size(pos_vectors, 2)) :: mask
  
      n = size(pos_vectors, 2)
      factor = pos_vectors(1, :)**2 + pos_vectors(2, :)**2 + pos_vectors(3, :)**2
      mask = sqrt(factor) < 4.0d0
  
      first_index = 0
      last_index = 0
      do i = 1, n
        if (mask(i) .and. first_index == 0) first_index = i
        if (mask(n + 1 - i) .and. last_index == 0) last_index = n + 1 - i
        if (first_index /= 0 .and. last_index /= 0) exit
      end do
    end subroutine eddy_range
  
    subroutine velocity_generator(xaxis_trimmed, factor, u, v, w)
      real(kind=8), intent(in) :: xaxis_trimmed(:, :), factor(:)
      real(kind=8), intent(out) :: u(:), v(:), w(:)
  
      u = -xaxis_trimmed(2, :) * exp(-factor * 2.0d0)
      v = xaxis_trimmed(1, :) * exp(-factor * 2.0d0)
      w = 0.0d0
    end subroutine velocity_generator
  
    subroutine rotation_total(thetax, thetay, thetaz, R)
      real(kind=8), intent(in) :: thetax, thetay, thetaz
      real(kind=8), intent(out) :: R(3, 3)
  
      R(1, 1) = cos(thetay) * cos(thetaz)
      R(1, 2) = -cos(thetay) * sin(thetaz)
      R(1, 3) = sin(thetay)
  
      R(2, 1) = cos(thetax) * sin(thetaz) + sin(thetax) * sin(thetay) * cos(thetaz)
      R(2, 2) = cos(thetax) * cos(thetaz) - sin(thetax) * sin(thetay) * sin(thetaz)
      R(2, 3) = -sin(thetax) * cos(thetay)
  
      R(3, 1) = sin(thetax) * sin(thetaz) - cos(thetax) * sin(thetay) * cos(thetaz)
      R(3, 2) = sin(thetax) * cos(thetaz) + cos(thetax) * sin(thetay) * sin(thetaz)
      R(3, 3) = cos(thetax) * cos(thetay)
    end subroutine rotation_total
  
    subroutine total_velocities(Nx, x_boundary, N_E, theta_list, a_list, u_total, v_total, w_total)
      integer, intent(in) :: Nx, N_E
      real(kind=8), intent(in) :: x_boundary
      real(kind=8), intent(in) :: theta_list(3, N_E)
      real(kind=8), intent(in) :: a_list(3, N_E)
      real(kind=8), intent(out) :: u_total(:), v_total(:), w_total(:)
  
      real(kind=8), dimension(3, Nx) :: xaxis
      real(kind=8), dimension(3, :) :: eddy_pos_translated, eddy_pos_trimmed, eddy_pos_rotated
      real(kind=8), dimension(3, 3) :: R, R_inv
      real(kind=8), dimension(:) :: factor, u, v, w
      integer :: i, first_index, last_index
  
      xaxis = sensor_line(x_boundary, Nx)
      u_total = 0.0d0
      v_total = 0.0d0
      w_total = 0.0d0
  
      do i = 1, N_E
        call eddy_range(xaxis - reshape(a_list(:, i), [3, 1]), first_index, last_index, factor)
        eddy_pos_trimmed = (xaxis(:, first_index:last_index))
        call rotation_total(theta_list(1, i), theta_list(2, i), theta_list(3, i), R)
        R_inv = transpose(R)
  
        eddy_pos_rotated = matmul(R, eddy_pos_trimmed)
        call velocity_generator(eddy_pos_rotated, factor, u, v, w)
  
        u_total(first_index:last_index) = u_total(first_index:last_index) + matmul(R_inv(:, 1), u)
        v_total(first_index:last_index) = v_total(first_index:last_index) + matmul(R_inv(:, 2), v)
        w_total(first_index:last_index) = w_total(first_index:last_index) + matmul(R_inv(:, 3), w)
      end do
    end subroutine total_velocities
  
  end module velocity_module
  