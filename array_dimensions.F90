module array_dimensions
    implicit none

    ! This module defines array dimensions at compile time
    ! kLength is set via preprocessor directive -DkLength=<value>
    ! Example: gfortran -DkLength=144 ...

#ifndef kLength
#error "kLength must be defined at compile time using -DkLength=<value>"
#endif

    integer, parameter :: dp = 8
    integer, parameter :: KLENGTH_PARAM = kLength

end module array_dimensions
