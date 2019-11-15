!!<summary>A module for constructing Hamiltonian matrix in k space representation.</summary>
module ClassHamiltonianConstructor
    use ClassNanotubeBundle
    use Class_kBase
    use Lattice2d
    use LatticeVector
    use ClassVectorInteger
    use ClassVectorComplex8
    use ClassVectorIntegerArray1d
    use ClassInterlayerAttenuation

    !!<member name="PI">\pi.</member>
    doubleprecision, parameter, private :: PI = 4.0 * atan(1.0)
    !!<summary>A type for constructing Hamiltonian matrix in k space representation.</summary>
    type HamiltonianConstructor
       !!<member name="hbarV">Constant. hbar times band velocity.</member>
        doubleprecision, private :: hbarV = 0.525205
        !!<member name="u0">Constant coefficient for u(K, \delta(r)).</member>
        doubleprecision, private :: u0 = 0.103
        !!<member name="nCutOff">The cutoff limit for n.</member>
        integer, private :: nCutoff = 4
        !!<member name="xCutoff">The range where a k point of one tube couples with another of the other tube.</member>
        doubleprecision, private :: xCutoff = 4.0
    contains
        procedure, public :: set_xCutoff
        procedure, public :: sethbarV
        procedure, public :: setU0
        procedure, public :: set_nCutOff
!        procedure, private :: getDisplacement
!        procedure, private :: getGaussianCenter
        procedure, private :: flipReducedCoordinateX
        procedure, public :: constructHamiltonianMatrix
        procedure, private :: couplingStrength
        procedure, private :: sumOfPhase
        procedure, public :: get_xCutoff
    end type HamiltonianConstructor
contains
    subroutine set_xCutoff(this, xCutoff)
        class(HamiltonianConstructor) :: this
        doubleprecision, intent(in) :: xCutoff

        this%xCutoff = xCutoff
    end subroutine

    function get_xCutoff(this) result(xCutoff)
        class(HamiltonianConstructor) :: this
        doubleprecision :: xCutoff

        xCutoff = this%xCutoff
    end function

    !!<summary>Fourier transform of a Gaussian not centering at 0 will results in a phase factor including a Gaussian. This is the summation of 6 phases for a infinite bundle.</summary>
    function sumOfPhase(this, kx, kx2, chiralVectorLength) result(s)
        class(HamiltonianConstructor) :: this
        doubleprecision, intent(in) :: kx, kx2
        doubleprecision, intent(in) :: chiralVectorLength
        complex(8) :: s
        doubleprecision, dimension(6) :: x, c

        integer :: i

        x = chiralVectorLength * (/ -1.0 / 3.0, 0.0, 1.0 / 3.0, -1.0 / 3.0, 0.0, 1.0 / 3.0 /)
        do i = 1, 6
            c(i) = chiralVectorLength * dble(2 * i - 1) / 12.0
        end do

        s = 0
        do i = 1, 6
            s = s + exp(cmplx(0.0, c(i) * kx + x(i) * kx2, kind=kind(1.0d0)))
        end do
    end function

    !!<summary>Fourier transformation of t(r).</summary>
    function couplingStrength(this, r) result(u)
        class(HamiltonianConstructor) :: this
        doubleprecision, intent(in) :: r
        doubleprecision :: u

        u = -0.0000514674 * r**3 + 0.00444809 * r**2 - 0.129352 * r + 1.27828
!        u = 0.11
    end function

    !!<summary>Set the value of hbar times band velocity.</summary>
    subroutine sethbarV(this, hbarV)
        class(HamiltonianConstructor) :: this
        doubleprecision, intent(in) :: hbarV

        this%hbarV = hbarV
    end subroutine sethbarV

    !!<summary>Set u_0, the constant coefficient for u(K, \delta(r)).</summary>
    subroutine setU0(this, u0)
        class(HamiltonianConstructor) :: this
        doubleprecision, intent(in) :: u0

        this%u0 = u0
    end subroutine setU0

    !!<summary>Set the cutoff limit for n.</summary>
    subroutine set_nCutOff(this, nCutOff)
        class(HamiltonianConstructor) :: this
        integer, intent(in) :: nCutoff

        this%nCutoff = nCutoff
    end subroutine set_nCutOff

    !!<summary>Calculate the displacement of two planes when two tubes are flattened.</summary>
!    function getDisplacement(this, bundle, bondIndex) result(displacement)
!        class(HamiltonianConstructor) :: this
!        type(NanotubeBundle), intent(in) :: bundle
!        integer, intent(in) :: bondIndex
!        doubleprecision, dimension(2) :: displacement
!
!        ! local variables
!        type(Bond) :: tmpBond
!        type(SingleWallNanotube) :: nanotube1, nanotube2
!        doubleprecision :: theta1
!        doubleprecision :: theta2
!
!        tmpBond = bundle%getBond(bondIndex)
!        nanotube1 = bundle%getNanotube(tmpBond%indexTube1)
!        nanotube2 = bundle%getNanotube(tmpBond%indexTube2)
!        theta1 = nanotube1%orientation
!        theta2 = nanotube2%orientation
!
!        displacement = (/ (theta1 + theta2 - PI - 2.0 * tmpBond%angleFromX) / (2.0 * PI) * bundle%getChiralVectorLength(), nanotube2%shiftZ - nanotube1%shiftZ /)
!    end function

    !!<summary>Calculate the center of Gaussian function describing the intertube interaction area.</summary>
!    function getGaussianCenter(this, bundle, bondIndex) result(center)
!        class(HamiltonianConstructor) :: this
!        type(NanotubeBundle), intent(in) :: bundle
!        integer, intent(in) :: bondIndex
!        doubleprecision :: center
!
!        ! local variables
!        type(Bond) :: tmpBond
!        type(SingleWallNanotube) :: nanotube1
!
!        tmpBond = bundle%getBond(bondIndex)
!        nanotube1 = bundle%getNanotube(tmpBond%indexTube1)
!
!        center = (nanotube1%orientation - tmpBond%angleFromX) / (2.0 * PI) * bundle%getChiralVectorLength()
!    end function

    !!<summary>Flip a k base reduced coordinate in x direction.</summary>
    function flipReducedCoordinateX(this, reducedCoordinate) result(res)
        class(HamiltonianConstructor) :: this
        integer, dimension(2), intent(in) :: reducedCoordinate
        integer, dimension(2) :: res

        res = (/reducedCoordinate(1), reducedCoordinate(1) - reducedCoordinate(2)/)
    end function

    !!<summary>Given a set of bases in k space, construct the hamiltonian matrix with respect to these bases. DO NOT FORGET TO DEALLOCATE!</summary>
    subroutine constructHamiltonianMatrix(this, bundle, base1, base2, xi, matrix, interlayerInteraction, couplingComponentList)
        implicit none
        class(HamiltonianConstructor) :: this
        !!<parameter name="bundle">Targeted double nanotube bundle.</parameter>
        type(NanotubeBundle), intent(in) :: bundle
        !!<parameter name="base">K bases.</parameter>
        type(kBase) :: base1, base2
        !!<parameter name="xi">\xi indicates whether K point(+1) or K' point(-1).</parameter>
        integer, intent(in) :: xi
        !!<parameter name="interlayerInteraction">Whether to switch on the interlayer interaction.</parameter>
        doubleprecision, intent(in) :: interlayerInteraction
        !!<parameter name="matix">result.<\parameter>
        complex(8), dimension(:, :), allocatable, intent(inout) :: matrix
        !!<parameter name="couplingComponentList">If it contains (/m, n/), coupling is active between k and k + m * G^M_1 + n * G^M_2<\parameter>
        type(VectorIntegerArray1d), intent(in) :: couplingComponentList

        ! local variables
        type(VectorInteger) :: kVector
!        type(Bond) :: tmpBond
        type(SingleWallNanotube) :: tmpNanotube
        type(InterlayerAttenuation) :: attenuation
        doubleprecision, dimension(2, 2) :: b
        doubleprecision, dimension(2) :: kk, kPoint, displacement, tmpVector, b1, b2
        doubleprecision :: theta, gaussianCenter
        complex(8) :: h_AB, u
        integer :: i, ii, j, k, row, column, y
        !!<local name="l">Number of bases in base.<local>
        integer :: l1, l2
        !!<local name="ratio">The ratio of |G^M_2| to 2 * pi / |C|, equaling n_2 in this case.</local>
!        doubleprecision :: ratio
        integer :: tmpIndex
        integer, dimension(2) :: tmpIntArray
        !!<local name="interlayerCoefficient">Calculate the multiplication of several constant in advance to accelerate.<local>
        doubleprecision :: interlayerCoefficient
        !!<local name="interlayerExponentialCoefficient">Calculate the multiplication of several constant in the exponent in advance to accelerate.</local>
        doubleprecision :: interlayerExponentialCoefficient
        doubleprecision :: umax, tmpkx, tmpkx2
        !!<local name="omega">A phase factor.</local>
        complex(8), parameter :: omega = cmplx(cos(2.0 * PI / 3.0), sin(2.0 * PI / 3.0))
!        ratio = bundle%get_n()
        tmpNanotube = bundle%getNanotube(1)
        interlayerCoefficient = interlayerInteraction / sqrt(tmpNanotube%getChiralVectorLength())
        umax = bundle%getHalfWidth() / sqrt(tmpNanotube%getChiralVectorLength())
        tmpNanotube = bundle%getNanotube(2)
        interlayerCoefficient = interlayerCoefficient / sqrt(tmpNanotube%getChiralVectorLength())
        interlayerExponentialCoefficient = - (bundle%getHalfWidth() ** 2) / 2.0
        umax = umax / sqrt(tmpNanotube%getChiralVectorLength())
!        print *, __FILE__, __LINE__, 'xCutOff:', this%xCutoff
!        this%xCutoff = 0.01
!        print *, __FILE__, __LINE__, umax, tmpNanotube%getChiralVectorLength()

        call reciprocalLatticeVector(a1, a2, b1, b2)


        attenuation = bundle%getAttenuation()
!        print *, __FILE__, __LINE__, interlayerCoefficient, bundle%getHalfWidth(), tmpNanotube%getChiralVectorLength()

        l1 = base1%getQuantity()
        l2 = base2%getQuantity()
        row = 2 * (l1 + l2)
        column = 2 * (l1 + l2)

        allocate(matrix(row, column))
        ! Initialize the matrix to zero
        do i = 1, row
            do j = 1, column
                matrix(i, j) = cmplx(0, 0)
            end do
        end do

        ! check whether \xi = +1 or -1
        if ((xi /= 1) .and. (xi /= -1)) then
            print *, __FILE__, __LINE__, "ERROR! \xi must be +1 or -1."
            stop
        end if


        !------------------------------------------------------------------------!
        !         |k, A_1> |k, B_1> |k, A_2> |k, B_2>  . . .  |k, A_n> |k, B_n>  !
        !         +----------------+-----------------+-------+-----------------+ !
        ! |k, A_1>|                |                 |       |                 | !
        !         |      H_1       | \dagger{U_{12}} | . . . | \dagger{U_{n1}} | !
        ! |k, B_1>|                |                 |       |                 | !
        !         +----------------+-----------------+       +-----------------+ !
        ! |k, A_2>|                |                 |       |                 | !
        !         |     U_{12}     |       H_2       | . . . | \dagger{U_{n2}} | !
        ! |k, B_2>|                |                 |       |                 | !
        !         +----------------+-----------------+       +-----------------+ !
        !     .   |       .                 .           .              .       | !
        !     .   |       .                 .             .            .       | !
        !     .   |       .                 .               .          .       | !
        !         +----------------+-----------------+       +-----------------+ !
        ! |k, A_n>|                |                 |       |                 | !
        !         |     U_{n1}     |      U_{n1}     | . . . |       H_n       | !
        ! |k, B_n>|                |                 |       |                 | !
        !         +----------------+-----------------+-------+-----------------+ !
        !------------------------------------------------------------------------!

        gaussianCenter = 0
        !------------------------!
        ! Intralayer Hamiltonian !
        !------------------------!
        tmpNanotube = bundle%getNanotube(1)
        kPoint = tmpNanotube%get_k()
        theta = tmpNanotube%getTheta()
        do j = 1, l1
            kk = base1%getCoordinate(j)
            h_AB = -this%hbarV * cmplx(kk(1) - kPoint(1), kPoint(2) - kk(2)) * cmplx(cos(-theta), sin(-theta))
            matrix(j, l1 + j) = h_AB
!            matrix(l1 + j, j) = conjg(h_AB)
        end do

        tmpNanotube = bundle%getNanotube(2)
        kPoint = tmpNanotube%get_k()
        theta = tmpNanotube%getTheta()
        do j = 1, l2
            kk = base2%getCoordinate(j)
            h_AB = -this%hbarV * cmplx(kk(1) - kPoint(1), kPoint(2) - kk(2)) * cmplx(cos(-theta), sin(-theta))
            matrix(2 * l1 + j, 2 * l1 + l2 + j) = h_AB
!            matrix(2 * l1 + l2 + j, 2 * l1 + j) = conjg(h_AB)
        end do

        !------------------------!
        ! Interlayer Hamiltonian !
        !------------------------!
!        print *, __FILE__, __LINE__, interlayerCoefficient, interlayerExponentialCoefficient, bundle%getHalfWidth()
!        print *, __FILE__, __LINE__, bundle%get_b1()
!        print *, __FILE__, __LINE__, bundle%get_b2()
!        do ii = 1, couplingComponentList%getLength()
!            tmpIntArray = couplingComponentList%get(ii)
!            tmpVector = tmpIntArray(1) * bundle%get_b1() + tmpIntArray(2) * bundle%get_b2()
!            print *, __FILE__, __LINE__, ii, tmpIntArray(1), tmpIntArray(2), bundle%get_kBar() + tmpVector, modulus2d(bundle%get_kBar() + tmpVector), this%couplingStrength(modulus2d(bundle%get_kBar() + tmpVector))
!        end do
        umax = 0
        do j = 1, l1
            kk = base1%getCoordinate(j)
!            print *, __FILE__, __LINE__, j, base1%getIndex2(kk)

            do ii = 1, couplingComponentList%getLength()
                tmpIntArray = couplingComponentList%get(ii)
                tmpVector = tmpIntArray(1) * bundle%get_b1() + tmpIntArray(2) * bundle%get_b2()
                this%u0 = this%couplingStrength(modulus2d(kk + tmpVector))
                tmpVector = tmpIntArray(1) * bundle%get_gm1() + tmpIntArray(2) * bundle%get_gm2()
                kVector = base2%kPointInRangeX(kk(1) + dot_product((/1.0, 0.0/), tmpVector), this%xCutoff)
                do k = 1, kVector%getLength()
                    tmpkx = base2%getCoordinateX(kVector%get(k)) - kk(1) - dot_product((/1.0, 0.0/), tmpVector)
                    tmpkx2 = base2%getCoordinateX(kVector%get(k)) + dot_product((/1.0, 0.0/), tmpIntArray(1) * bundle%get_b1() + tmpIntArray(2) * bundle%get_b2())
                    u = this%u0 * interlayerCoefficient * attenuation%attenuationFunctionFT(tmpkx) * &
                        (1 + bundle%isInfiniteInt() * (this%sumOfPhase(tmpkx, tmpkx2, tmpNanotube%getChiralVectorLength()) - 1))
                    tmpIndex = base2%getIndex2((/base2%getCoordinateX(kVector%get(k)), kk(2) + &
                                                                      dot_product((/0.0, 1.0/), tmpVector)/))
                    if(tmpIndex /= 0) then
                        ! U_{A_1A_2}
                        matrix(j     , tmpIndex + 2 * l1     ) = &
                        matrix(j     , tmpIndex + 2 * l1     ) + u
                        ! U_{A_1B_2}
                        matrix(j     , tmpIndex + 2 * l1 + l2) = &
                        matrix(j     , tmpIndex + 2 * l1 + l2) + u * omega**xi
                        ! U_{B_1A_2}
                        matrix(j + l1, tmpIndex + 2 * l1     ) = &
                        matrix(j + l1, tmpIndex + 2 * l1     ) + u * omega**(-xi)
                        ! U_{B_1B_2}
                        matrix(j + l1, tmpIndex + 2 * l1 + l2) = &
                        matrix(j + l1, tmpIndex + 2 * l1 + l2) + u
                    end if
                end do
            end do

!            ! \delta(k)
!            kVector = base2%kPointInRangeX(kk(1), this%xCutoff)
!            this%u0 = this%couplingStrength(modulus2d(kk))
!            do k = 1, kVector%getLength()
!                u = this%u0 * interlayerCoefficient * exp(interlayerExponentialCoefficient * &
!                    (base2%getCoordinateX(kVector%get(k)) - kk(1))**2)
!!                y = floor((kk(2) - dot_product((/0.0, 1.0/), base2%getShift() + base2%getCenter())) / base2%get_yStep() + 0.5)
!                tmpIndex = base2%getIndex2((/base2%getCoordinateX(kVector%get(k)), kk(2)/))
!                if(tmpIndex /= 0) then
!                    ! U_{A_1A_2}
!                    matrix(j     , tmpIndex + 2 * l1     ) = &
!                    matrix(j     , tmpIndex + 2 * l1     ) + 6.0 * conjg(u * exp(cmplx(0, -kk(1) * gaussianCenter, kind=kind(1.0d0))))
!                    ! U_{A_1B_2}
!                    matrix(j     , tmpIndex + 2 * l1 + l2) = &
!                    matrix(j     , tmpIndex + 2 * l1 + l2) + 6.0 * conjg(u * exp(cmplx(0, -kk(1) * gaussianCenter, kind=kind(1.0d0))))
!                    ! U_{B_1A_2}
!                    matrix(j + l1, tmpIndex + 2 * l1     ) = &
!                    matrix(j + l1, tmpIndex + 2 * l1     ) + 6.0 * conjg(u * exp(cmplx(0, -kk(1) * gaussianCenter, kind=kind(1.0d0))))
!                    ! U_{B_1B_2}
!                    matrix(j + l1, tmpIndex + 2 * l1 + l2) = &
!                    matrix(j + l1, tmpIndex + 2 * l1 + l2) + 6.0 * conjg(u * exp(cmplx(0, -kk(1) * gaussianCenter, kind=kind(1.0d0))))
!                end if
!            end do
!
!            ! \delta(k + G^M_2)
!            kVector = base2%kPointInRangeX(kk(1) + dot_product((/1.0, 0.0/), bundle%get_gm2()), this%xCutoff)
!            this%u0 = this%couplingStrength(modulus2d(kk + bundle%get_gm2()))
!            do k = 1, kVector%getLength()
!                u = this%u0 * interlayerCoefficient * exp(interlayerExponentialCoefficient * &
!                    (base2%getCoordinateX(kVector%get(k)) - kk(1) - dot_product((/1.0, 0.0/),bundle%get_gm2()))**2)
!!                y = floor((kk(2) - dot_product((/0.0, 1.0/), base2%getShift()  + &
!!                                                              base2%getCenter() - &
!!                                                              bundle%get_gm2())) / base2%get_yStep() + 0.5)
!                tmpIndex = base2%getIndex2((/base2%getCoordinateX(kVector%get(k)), kk(2) + dot_product((/0.0, 1.0/), bundle%get_gm2())/))
!!                print *, __FILE__, __LINE__, kk
!!                if(tmpIndex == 0) then
!!                    print *, __FILE__, __LINE__, 'a'
!!                else
!!                    print *, __FILE__, __LINE__, base2%getCoordinate(tmpIndex)
!!                end if
!!                print *, __FILE__, __LINE__, bundle%get_gm2()
!                if(tmpIndex /= 0) then
!                    tmpVector = bundle%get_gm2()
!                    ! U_{A_1A_2}
!                    matrix(j     , tmpIndex + 2 * l1     ) = &
!                    matrix(j     , tmpIndex + 2 * l1     ) + 6.0 * u * exp(cmplx(0, -(kk(1) + tmpVector(1)) * gaussianCenter, kind=kind(1.0d0))) * exp(cmplx(0, xi * dot_product(b1, displacement), kind=kind(1.0d0)))
!                    ! U_{A_1B_2}
!                    matrix(j     , tmpIndex + 2 * l1 + l2) = &
!                    matrix(j     , tmpIndex + 2 * l1 + l2) + 6.0 * u * exp(cmplx(0, -(kk(1) + tmpVector(1)) * gaussianCenter, kind=kind(1.0d0))) * exp(cmplx(0, xi * dot_product(b1, displacement), kind=kind(1.0d0))) * omega**xi
!                    ! U_{B_1A_2}
!                    matrix(j + l1, tmpIndex + 2 * l1     ) = &
!                    matrix(j + l1, tmpIndex + 2 * l1     ) + 6.0 * u * exp(cmplx(0, -(kk(1) + tmpVector(1)) * gaussianCenter, kind=kind(1.0d0))) * exp(cmplx(0, xi * dot_product(b1, displacement), kind=kind(1.0d0))) * omega**(-xi)
!                    ! U_{B_1B_2}
!                    matrix(j + l1, tmpIndex + 2 * l1 + l2) = &
!                    matrix(j + l1, tmpIndex + 2 * l1 + l2) + 6.0 * u * exp(cmplx(0, -(kk(1) + tmpVector(1)) * gaussianCenter, kind=kind(1.0d0))) * exp(cmplx(0, xi * dot_product(b1, displacement), kind=kind(1.0d0)))
!                end if
!            end do
!
!            ! \delta(k + G^M_1 + G^M_2)
!            kVector = base2%kPointInRangeX(kk(1) + dot_product((/1.0, 0.0/), bundle%get_gm2()) &
!                                                 + dot_product((/1.0, 0.0/), bundle%get_gm1()), this%xCutoff)
!            this%u0 = this%couplingStrength(modulus2d(kk + bundle%get_gm1() + bundle%get_gm2()))
!            do k = 1, kVector%getLength()
!                u = this%u0 * interlayerCoefficient * exp(interlayerExponentialCoefficient * &
!                    (base2%getCoordinateX(kVector%get(k)) - kk(1) - dot_product((/1.0, 0.0/),bundle%get_gm1()) - dot_product((/1.0, 0.0/),bundle%get_gm2()))**2)
!                if(abs(u) > umax) then
!                    umax = abs(u)
!                end if
!!                y = floor((kk(2) - dot_product((/0.0, 1.0/), base2%getShift()  + &
!!                                                              base2%getCenter() - &
!!                                                              bundle%get_gm1()  - &
!!                                                              bundle%get_gm2()))/ base2%get_yStep() + 0.5)
!                tmpIndex = base2%getIndex2((/base2%getCoordinateX(kVector%get(k)), kk(2) + dot_product((/0.0, 1.0/), bundle%get_gm1() + bundle%get_gm2())/))
!                if(tmpIndex /= 0) then
!                    tmpVector = bundle%get_gm1() + bundle%get_gm2()
!                    ! U_{A_1A_2}
!                    matrix(j     , tmpIndex + 2 * l1     ) = &
!                    matrix(j     , tmpIndex + 2 * l1     ) + 6.0 * u * exp(cmplx(0, -(kk(1) + tmpVector(1)) * gaussianCenter, kind=kind(1.0d0))) * exp(cmplx(0, xi * dot_product(b1 + b2, displacement), kind=kind(1.0d0)))
!                    ! U_{A_1B_2}
!                    matrix(j     , tmpIndex + 2 * l1 + l2) = &
!                    matrix(j     , tmpIndex + 2 * l1 + l2) + 6.0 * u * exp(cmplx(0, -(kk(1) + tmpVector(1)) * gaussianCenter, kind=kind(1.0d0))) * exp(cmplx(0, xi * dot_product(b1 + b2, displacement), kind=kind(1.0d0))) * omega**(-xi)
!                    ! U_{B_1A_2}
!                    matrix(j + l1, tmpIndex + 2 * l1     ) = &
!                    matrix(j + l1, tmpIndex + 2 * l1     ) + 6.0 * u * exp(cmplx(0, -(kk(1) + tmpVector(1)) * gaussianCenter, kind=kind(1.0d0))) * exp(cmplx(0, xi * dot_product(b1 + b2, displacement), kind=kind(1.0d0))) * omega**xi
!                    ! U_{B_1B_2}
!                    matrix(j + l1, tmpIndex + 2 * l1 + l2) = &
!                    matrix(j + l1, tmpIndex + 2 * l1 + l2) + 6.0 * u * exp(cmplx(0, -(kk(1) + tmpVector(1)) * gaussianCenter, kind=kind(1.0d0))) * exp(cmplx(0, xi * dot_product(b1 + b2, displacement), kind=kind(1.0d0)))
!                end if
!            end do
        end do
!        print *, __FILE__, __LINE__, umax
    end subroutine constructHamiltonianMatrix
end module ClassHamiltonianConstructor
