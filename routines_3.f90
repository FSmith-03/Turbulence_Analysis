
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
  real, dimension(:,:), allocatable :: R, xaxis_trimmed, xaxis_rotated, velocities, a_list_T, velocity_total_trans, velocity_pre
  integer :: i, j, N, x_boundary, Nxf, N_E, first_index, last_index, ARP, ARN, L_e, test
  real :: tol, thetax, thetay, thetaz
  real, dimension(:), allocatable :: N_half_list, factor, x
  real, dimension(3,3) :: R_inv
  integer, dimension(6) :: ARN_list

  ! Print input variables
  write(*,*) "Input Variables", input_ints

  ! Initialize variables
  x_boundary = input_ints(1)
  Nxf = input_ints(2)
  N_E = input_ints(3)
  tol = 0.05
  ARP = 1
  ARN_list = (/1, 2, 5, 10, 20, 40/)
  test = input_ints(4)
  L_e = 1
  write(*,*) "Maximum x value in a_list", maxval(a_list(1, :))

  ! Generate sensor line
  if (.not. allocated(pos_vector)) then
    allocate(pos_vector(3, Nxf))
  end if
  if (.not. allocated(x)) then
    allocate(x(Nxf))
  end if
  do i = 1, Nxf
    x(i) = 2.0 * x_boundary * (i - 1) / (Nxf - 1) - x_boundary
  end do
  pos_vector(1, :) = x
  pos_vector(2, :) = 0.0
  pos_vector(3, :) = 0.0

  ! Allocate arrays
  if (.not. allocated(pos_vector_translated)) then
    allocate(pos_vector_translated(3, Nxf))
  end if
  if (.not. allocated(u_total_0)) then
    allocate(u_total_0(Nxf), v_total_0(Nxf), w_total_0(Nxf))
  end if

  ! Initialize total velocity arrays

  write(*,*) "Initial Variables Allocated"

  ! Transpose a_list
  a_list_T = transpose(a_list)
  write(*,*) 'N_E', N_E
  if (.not. allocated(N_half_list)) then
    allocate(N_half_list(N_E))
  end if
  N_half_list = ((x_boundary + a_list_T(:,1)) * Nxf) / (2 * x_boundary)
  write(*,*) 'Min N_half', minval(N_half_list)
  
    do i = 1, Nxf
      u_total_0(i) = 0.0
      v_total_0(i) = 0.0
      w_total_0(i) = 0.0
    end do
    
    ! Set total velocities to initial values of zero
    u_total = u_total_0
    v_total = v_total_0
    w_total = w_total_0

    ! ARN Value from ARN_list
    ARN = ARN_list(test)
    write(*,*) "ARN", ARN


    write(*,*) "Eddy Loop Start"

    ! Loop over eddies
    do i = 1, N_E

      ! Print progress for every 1000th iteration
      if (mod(i, 1000) == 0) then
        write(*,*) "Progress", real(i) / real(N_E) * 100.0, "%"
      end if

      ! Translate position vectors
      do j = 1, Nxf
        pos_vector_translated(:, j) = pos_vector(:, j) - a_list(:, i)
      end do

      ! Trim position vectors
      first_index = int(N_half_list(i) - (3.5/tol) * ARN)
      if (first_index < 0) then
        first_index = 0
      end if
      last_index = int(N_half_list(i) + (3.5/tol) * ARN)
      if (last_index > Nxf) then
        last_index = Nxf
      end if
      N = last_index - first_index + 1

      if (.not. allocated(velocity_pre)) then
        allocate(velocity_pre(3, N))
      end if

      xaxis_trimmed = pos_vector_translated(:, first_index:last_index)
      
      ! Rotate position vectors
      if (.not. allocated(R)) then
        allocate(R(3, 3))
      end if
      thetax = theta_list(1, i)
      thetay = theta_list(2, i)
      thetaz = theta_list(3, i)
      R(1,1) = cos(thetay)*cos(thetaz)
      R(1,2) = cos(thetay)*sin(thetaz)
      R(1,3) = -sin(thetay)
      R(2,1) = sin(thetax)*sin(thetay)*cos(thetaz) - cos(thetax)*sin(thetaz)
      R(2,2) = sin(thetax)*sin(thetay)*sin(thetaz) + cos(thetax)*cos(thetaz)
      R(2,3) = sin(thetax)*cos(thetay)
      R(3,1) = cos(thetax)*sin(thetay)*cos(thetaz) + sin(thetax)*sin(thetaz)
      R(3,2) = cos(thetax)*sin(thetay)*sin(thetaz) - sin(thetax)*cos(thetaz)
      R(3,3) = cos(thetax)*cos(thetay)
      R_inv = transpose(R)
      xaxis_rotated = matmul(R, xaxis_trimmed)
      
      ! Calculate velocities
      N = size(xaxis_trimmed, 2)
      if (.not. allocated(u)) then
        allocate(u(N), v(N), w(N))
      end if
      ! Factor for spherical eddies
      if (ARN == 1) then
      factor = -2 * (xaxis_trimmed(1, :)**2 + xaxis_trimmed(2, :)**2 + xaxis_trimmed(3, :)**2)
      ! Factor for pancake and needle eddies
      else
      factor = -2 * ( ((xaxis_trimmed(1, :)**2 + xaxis_trimmed(2, :)**2) / (ARP * L_e)**2) + ((xaxis_trimmed(3, :)**2) / (ARN * L_e)**2))
      end if
      u = -xaxis_trimmed(2, :) * exp(factor)
      v = xaxis_trimmed(1, :) * exp(factor)
      w = 0.0d0
      !u_0 = -y_r * np.exp(-2 * (((x_r**2 + y_r**2) / (ARP * L**2)) + (z_r**2 / (ARN * L)**2)))
      !v_0 = x_r * np.exp(-2 * (((x_r**2 + y_r**2) / (ARP * L**2)) + (z_r**2 / (ARN * L)**2)))

      velocity_pre(1, :) = u
      velocity_pre(2, :) = v
      velocity_pre(3, :) = w

      ! Rotate velocities back to base frame
      velocities = matmul(R_inv, velocity_pre)

      ! Update total velocities
      u_total(first_index:last_index) = u_total(first_index:last_index) + velocities(1, :)
      v_total(first_index:last_index) = v_total(first_index:last_index) + velocities(2, :)
      w_total(first_index:last_index) = w_total(first_index:last_index) + velocities(3, :)
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
    ! Fortran code
    write(*,*) "Output file name: ", "output.dat"
    OPEN(UNIT=10, FILE="output.dat", STATUS='REPLACE')
    WRITE(10, *) velocity_total
    CLOSE(10)
end subroutine main_calculation
