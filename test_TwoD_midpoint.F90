program test_TwoD_Midpoint
    use EDQNMstratifiedModules
    implicit none

    integer, parameter :: dp = 8
    integer :: kLength, kj, pj, qj
    double precision :: N, nu, D
    double precision :: t, t0
    double precision :: k
    double precision, allocatable :: pVals(:)
    double precision, allocatable :: E(:), ET(:), EHdir(:), EHpol(:), ETH(:), F(:)
    double precision, allocatable :: mu1(:), mu3(:)
    double precision, allocatable :: HPOL(:), HDIR(:), HT(:)
    double precision, allocatable :: ExtForcing(:)
    double precision, allocatable :: weight(:,:,:)
    double precision, allocatable :: centroidX(:,:,:), centroidY(:,:,:)
    integer, allocatable :: Q11(:,:,:,:)
    integer, allocatable :: triadFlag(:,:,:)
    integer, allocatable :: insideCutCell(:,:,:), outsideCutCell(:,:,:)
    double precision :: du, ySquaredAvg
    integer :: qjMin_dummy, qjMax_dummy   ! unused but required by interface
    double precision :: forcing(6)

    ! -----------------------
    ! Basic scalar parameters
    ! -----------------------
    kLength = 5
    N  = 1.0d0
    nu = 0.01d0
    D  = 0.02d0
    k  = 2.0d0
    t  = 1.0d0
    t0 = 0.1d0
    du = 1.0d0
    ySquaredAvg = 0.0d0
    qjMin_dummy = 1
    qjMax_dummy = kLength

    ! -----------------------
    ! Allocate arrays
    ! -----------------------
    allocate(pVals(kLength))
    allocate(E(kLength), ET(kLength), EHdir(kLength), EHpol(kLength), ETH(kLength), F(kLength))
    allocate(mu1(kLength), mu3(kLength))
    allocate(HPOL(kLength), HDIR(kLength), HT(kLength))
    allocate(ExtForcing(kLength))

    allocate(weight(kLength,kLength,kLength))
    allocate(centroidX(kLength,kLength,kLength))
    allocate(centroidY(kLength,kLength,kLength))

    allocate(Q11(2,kLength,kLength,kLength))
    allocate(triadFlag(kLength,kLength,kLength))
    allocate(insideCutCell(kLength,kLength,kLength))
    allocate(outsideCutCell(kLength,kLength,kLength))

    ! -----------------------
    ! Manufacture simple test data
    ! -----------------------
    do pj = 1, kLength
        pVals(pj) = dble(pj)   ! 1,2,3,4,5
        E(pj)     = 0.1d0*pj
        ET(pj)    = 0.2d0*pj
        EHdir(pj) = 0.01d0*pj
        EHpol(pj) = 0.015d0*pj
        ETH(pj)   = 0.005d0*pj
        F(pj)     = 0.03d0*pj
        ExtForcing(pj) = 0.001d0 * pj
    end do

    weight = 1.0d0
    centroidX = 1.0d0
    centroidY = 1.0d0
    Q11 = 1
    triadFlag = 1      ! everything active
    insideCutCell = 0
    outsideCutCell = 0

    ! -----------------------
    ! Compute mu1, mu3 from manufactured E, ET
    ! -----------------------
    call getMu(E, pVals, mu1)
    call getMu(ET, pVals, mu3)

    ! -----------------------
    ! Compute H-values (HPOL, HDIR, HT)
    ! -----------------------
    call getHvals(kLength, E, ET, EHdir, EHpol, ETH, HT, HDIR, HPOL)

    ! -----------------------
    ! Now test TwoD_MidpointTest for each kj
    ! -----------------------
    print *, "============================================="
    print *, "Testing TwoD_MidpointTest with kLength=", kLength
    print *, "============================================="

    do kj = 1, kLength
        call TwoD_MidpointTest( N, nu, D, kj, k, pVals, HPOL, HDIR, HT, E, F, ET, &
                                 mu1, mu3, t, weight, du, qjMin_dummy, qjMax_dummy, ySquaredAvg, &
                                 triadFlag, outsideCutCell, Q11, centroidX, centroidY, insideCutCell, &
                                 ExtForcing, t0, forcing )

        print *, "------ kj = ", kj, " ------"
        print *, "forcing = ", forcing
    end do

    print *, "Completed single-rank test."

end program test_TwoD_Midpoint
