module types

    type input_data

        real :: x_boundary
        integer :: Nxf
        real, dimension(:), allocatable :: first_and_last
    
    end type input_data

end module types