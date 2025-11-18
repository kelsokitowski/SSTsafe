module array_dimensions
    implicit none

    ! This module defines array dimensions at compile time
    ! KLENGTH_SIZE is set via preprocessor directive -DKLENGTH_SIZE=<value>
    ! Example: gfortran -DKLENGTH_SIZE=144 ...

#ifndef KLENGTH_SIZE
#error "KLENGTH_SIZE must be defined at compile time using -DKLENGTH_SIZE=<value>"
#endif

    integer, parameter :: dp = 8
    integer, parameter :: KLENGTH_PARAM = KLENGTH_SIZE

end module array_dimensions
