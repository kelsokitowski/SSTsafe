module array_dimensions
    implicit none

    ! This module defines array dimensions at compile time
    ! kLength is set via preprocessor directive -DkLength=<value>
    ! Example: mpifort -DkLength=144 ...

    ! Note: If kLength is not defined, compilation will fail with error about kLength
    integer, parameter :: KLENGTH_PARAM = kLength
    integer, parameter :: dp = 8

end module array_dimensions
