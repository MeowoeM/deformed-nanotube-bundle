module ClassNanotubeBundle
    use LatticeVector
    use Lattice2d
    use ClassSingleWallNanotube
    use ClassInterlayerAttenuation

    !!<member name="PI">\pi.</member>
    doubleprecision, parameter, private :: PI = 4.0 * atan(1.0)

    !!<summary>Nanotube bundle class. All the nanotube can be cut and flattened to a rectangular graphene ribbon whose width is the chiral vector. when rolled up again and two sides of each tube sealed, let the sealed line face towards y dircetion. The chiral vectors of all the tubes in this bundle are supposed to be (m, n) or (m + n, -n)</summary>
    type NanotubeBundle
        !!<member name-"isInfinite">Whether it is an infinite bundle or one consisting of only 2 tubes.</member>
        logical, private :: isInfinite
        !!<member name="nanotubeArray">An array storing nanotubes.</member>
        type(SingleWallNanotube), dimension(2), private :: nanotubeArray
        !!<member name="theta">The angle from the zigzag direction to the chiral vector.</member>
        doubleprecision, private :: theta
        !!<member name="lm1">Superlattice vector.</member>
        doubleprecision, dimension(2), private :: lm1
        !!<member name="lm2">Superlattice vector.</member>
        doubleprecision, dimension(2), private :: lm2
        !!<member name="gm1">Reciprocal superlattice vector.</member>
        doubleprecision, dimension(2), private :: gm1
        !!<member name="gm2">Reciprocal superlattice vector.</member>
        doubleprecision, dimension(2), private :: gm2
        !!<member name="kBar">Average K point, (K_1 + K_2) / 2.</member>
        doubleprecision, dimension(2), private :: kBar
        !!<member name="kPrimeBar">Average K' point, (K'_1 + K'_2) / 2.</member>
        doubleprecision, dimension(2), private :: kPrimeBar
        !!<member name="angleFromX">the angle from x direction to the line connecting the centers of two tubes.</member>
        doubleprecision, private :: angleFromX
        !!<member name="halfWidth">Half width of the Gaussian function which controls the interaction strength between two layers.</member>
        doubleprecision, private :: halfWidth
        !!<memeber name="interlayerAttenuation">Coupling strength decaying between tubes.</member>
        type(InterlayerAttenuation) :: interlayerAttenuation


        doubleprecision, dimension(2), private :: b1, b2
    contains
        procedure, public :: initializeNanotubeBundle
        procedure, public :: isInfiniteInt
        procedure, public :: getNanotube
        procedure, public :: get_lm1
        procedure, public :: get_lm2
        procedure, public :: get_gm1
        procedure, public :: get_gm2
        procedure, public :: get_kBar
        procedure, public :: get_kPrimeBar
        procedure, public :: calcShift
        procedure, public :: set_gm1
        procedure, public :: get_b1
        procedure, public :: get_b2
        procedure, public :: setHalfWidth
        procedure, public :: getHalfWidth
        procedure, public :: setHalfWidthAuto
        procedure, public :: initAttenuation => initializeInterBundleAttenuation
        procedure, public :: genterate_kTable => genterate_kTableWrapper
        procedure, public :: getAttenuation
    end type

!    interface NanotubeBundle
!        procedure :: constructNanotubeBundle
!    end interface NanotubeBundle
contains
    function newNanotubeBundle(tube1, tube2, angleFromX, isInfinite) result(this)
        type(NanotubeBundle) :: this
        class(SingleWallNanotube) :: tube1, tube2
        doubleprecision :: angleFromX
        logical, intent(in) :: isInfinite

        call this%initializeNanotubeBundle(tube1, tube2, angleFromX, isInfinite)
    end function

    !!<summary>Constructor of NanotubeBundle type.</summary>
    subroutine initializeNanotubeBundle(this, tube1, tube2, angleFromX, isInfinite)
        class(NanotubeBundle) :: this
        class(SingleWallNanotube) :: tube1, tube2
        doubleprecision :: angleFromX
        logical, intent(in) :: isInfinite

        ! local variables
        ! G_i = (R^{-1} - R) b_i
        doubleprecision, dimension(2) :: b1, b2
        doubleprecision :: numerator1, numerator2, denominator
        !!<local name="tmp">[R(\theta_1) - R(-\theta_2)]^{-1}.</local>
        doubleprecision, dimension(2,2) :: tmp

        this%nanotubeArray(1) = tube1
        this%nanotubeArray(2) = tube2
        this%angleFromX = angleFromX
        this%isInfinite = isInfinite

        this%theta = tube1%getTheta() + tube2%getTheta()
        this%kBar = 1.0 / 2.0 * (tube1%get_k() + tube2%get_k())
        this%kPrimeBar = 1.0 / 2.0 * (tube1%get_kPrime() + tube2%get_kPrime())

        ! the superlattice vector is given by L_i^M = [R(\theta_1) - R(-\theta_2)]^{-1} a_i
        numerator1 = cos(tube1%getTheta()) - cos(tube2%getTheta())
        numerator2 = sin(tube1%getTheta()) + sin(tube2%getTheta())
        denominator = numerator1 ** 2 + numerator2 ** 2
        tmp = reshape((/numerator1 / denominator, &
                        numerator2 / denominator, &
                       -numerator2 / denominator, &
                        numerator1 / denominator/), shape(tmp))
        this%lm1 = matmul(tmp, a1)
        this%lm2 = matmul(tmp, a2)

        call reciprocalLatticeVector(a1, a2, b1, b2)
!        print *, tube1%getTheta(), tube2%getTheta()
        this%gm1 = rotate2d(b1, -tube1%getTheta()) - rotate2d(b1, -tube2%getTheta())
        this%gm2 = rotate2d(b2, -tube1%getTheta()) - rotate2d(b2, -tube2%getTheta())

        call reciprocalLatticeVector(rotate2d(a1, -tube1%getTheta()), &
                                     rotate2d(a2, -tube1%getTheta()), &
                                     b1, b2)
        this%b1 = b1
        this%b2 = b2
        print *, __FILE__, __LINE__, 'b1', b1
        print *, __FILE__, __LINE__, 'b2', b2

        call this%setHalfWidthAuto(tube1%getChiralVectorLength())
    end subroutine

    !!<summary>If infinite bundle, return 1 else 0.</summary>
    function isInfiniteInt(this) result(isInfinite)
        class(NanotubeBundle) :: this
        integer :: isInfinite

        if(this%isInfinite) then
            isInfinite = 1
        else
            isInfinite = 0
        end if
    end function

    !!<summary>Configuring attenuation function describing coupling strength decaying between two tubes.</summary>
    subroutine initializeInterBundleAttenuation(this, n, nw)
        class(NanotubeBundle) :: this
        integer, intent(in) :: n
        doubleprecision, intent(in) :: nw

        call this%interlayerAttenuation%init(n, nw, this%halfWidth)
    end subroutine

    !!<summary>Wrapper of InterlayerAttenuation%genterate_kTable</summary>
    subroutine genterate_kTableWrapper(this, kMin, kMax, nk)
        class(NanotubeBundle) :: this
        doubleprecision, intent(in) :: kMin, kMax
        integer, intent(in) :: nk

        call this%interlayerAttenuation%genterate_kTable(kMin, kMax, nk)
    end subroutine

    function getAttenuation(this) result(attenuation)
        class(NanotubeBundle) :: this
        type(InterlayerAttenuation) :: attenuation

        attenuation = this%interlayerAttenuation
    end function

    !!<summary>Set the half width of the Gaussian function which controls the interaction strength between two layers automatically.</summary>
    subroutine setHalfWidthAuto(this, chiralVectorLength)
        class(NanotubeBundle) :: this
        doubleprecision, intent(in) :: chiralVectorLength

        this%halfWidth = sqrt(chiralVectorLength / 2.0 / PI * 0.045 / 2.0)
    end subroutine

    !!<summary>Set the half width of the Gaussian function which controls the interaction strength between two layers.</summary>
    subroutine setHalfWidth(this, halfWidth)
        class(NanotubeBundle) :: this
        doubleprecision, intent(in) :: halfWidth

        this%halfWidth = halfWidth
    end subroutine setHalfWidth

    function getHalfWidth(this) result(halfWidth)
        class(NanotubeBundle) :: this
        doubleprecision :: halfWidth

        halfWidth = this%halfWidth
    end function

    !!<summary>Return the type(nanotube) object of a given index.</summary>
    function getNanotube(this, nanotubeIndex) result(nanotube)
        class(NanotubeBundle) :: this
        integer, intent(in) :: nanotubeIndex
        type(SingleWallNanotube) :: nanotube

        nanotube = this%nanotubeArray(nanotubeIndex)
    end function

    !!<summary>Return the superlattice vector L_1^M.</summary>
    function get_lm1(this) result(lm1)
        implicit none
        class(NanotubeBundle), intent(in) :: this
        doubleprecision, dimension(2) :: lm1

        lm1 = this%lm1
    end function get_lm1

    !!<summary>Return the superlattice vector L_2^M.</summary>
    function get_lm2(this) result(lm2)
        implicit none
        class(NanotubeBundle), intent(in) :: this
        doubleprecision, dimension(2):: lm2

        lm2 = this%lm2
    end function get_lm2

    !!<summary>Return the superlattice vector G_1^M.</summary>
    function get_gm1(this) result(gm1)
        implicit none
        class(NanotubeBundle), intent(in) :: this
        doubleprecision, dimension(2) :: gm1

        gm1 = this%gm1
    end function get_gm1

    !!<summary>Return the superlattice vector G_2^M.</summary>
    function get_gm2(this) result(gm2)
        implicit none
        class(NanotubeBundle), intent(in) :: this
        doubleprecision, dimension(2) :: gm2

        gm2 = this%gm2
    end function get_gm2

    !!<summary>Return the coordinate of \bar{K}</summary>
    function get_kBar(this) result(kBar)
        implicit none
        class(NanotubeBundle), intent(in) :: this
        doubleprecision, dimension(2) :: kBar

        kBar = this%kBar
    end function get_kBar

    !!<summary>Return the coordinate of \bar{K}'</summary>
    function get_kPrimeBar(this) result(kPrimeBar)
        implicit none
        class(NanotubeBundle), intent(in) :: this
        doubleprecision, dimension(2) :: kPrimeBar

        kPrimeBar = this%kPrimeBar
    end function get_kPrimeBar

    !!<summary>Given the center of a cutoff circle, calculate its distance to the nearest cutting line in k_x direction.</summary>
    function calcShift(this, center, tubeIndex) result(shift)
        class(NanotubeBundle) :: this
        !!<parameter name="centerX">k_x component of the center of cutoff circle</parameter>
        doubleprecision, dimension(2), intent(in) :: center
        !!<paremeter name="tubeIndex">Whether the shiftX of tube1 or tube2's is to return.</parameter>
        integer, intent(in) :: tubeIndex

        doubleprecision, dimension(2) :: shift

        ! local variable
        !!<local name="cuttingLineSpaceing"></local>
        doubleprecision :: cuttingLineSpacing

        ! check if tubeIndex is valid
        if((tubeIndex==1) .or. (tubeIndex==2)) then
            cuttingLineSpacing = 8.0 * atan(1.0) / this%nanotubeArray(tubeIndex)%getChiralVectorLength()
            shift = (/-mod(center(1), cuttingLineSpacing), dble(0.0)/)
        else
            print *, __FILE__, __LINE__, "ERROR! Tube index must be 1 or 2."
        end if
    end function

    subroutine set_gm1(this, gm1)
        class(NanotubeBundle) :: this
        doubleprecision, dimension(2), intent(in) :: gm1

        this%gm1 = gm1
    end subroutine

    function get_b1(this) result(b1)
        class(NanotubeBundle) :: this
        doubleprecision, dimension(2) :: b1

        b1 = this%b1
    end function

    function get_b2(this) result(b2)
        class(NanotubeBundle) :: this
        doubleprecision, dimension(2) :: b2

        b2 = this%b2
    end function
end module

