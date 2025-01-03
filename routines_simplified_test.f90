program main
    implicit none
    interface
      subroutine main_calculation(input_ints, a_list, theta_list, u_total, v_total, w_total)
        implicit none
        integer, intent(in), dimension(:) :: input_ints
        real, dimension(:,:) :: a_list, theta_list
        real, intent(out), dimension(:), allocatable :: u_total, v_total, w_total
      end subroutine main_calculation
    end interface
    integer, dimension(3) :: input_ints
    real, dimension(:,:), allocatable :: a_list, theta_list
    real, dimension(:), allocatable :: u_total, v_total, w_total
    integer :: i
    input_ints = [10, 1000, 10]
    write(*,*) "Input Variables", input_ints
    allocate(a_list(3, input_ints(3)), theta_list(3, input_ints(3)))
    do i = 1, input_ints(3)
      a_list(:, i) = [1.0, 1.0, 1.0]
      theta_list(:, i) = [0.0, 0.0, 0.0]
    end do
    call main_calculation(input_ints, a_list, theta_list, u_total, v_total, w_total)
    write(*,*) "Total Velocities", u_total, v_total, w_total
  end program main