module ClassBundleOfDoubleNanotube
    use ClassSingleWallNanotube
    use Lattice2d
    use LatticeVector
    use Class_kBase

    !!<summary>Bundle of double nanotubes class.</summary>
    type BundleOfDoubleNanotube
        !!<member name="nanotube1">One nanotube in the Bundle.</member>
        type(SingleWallNanotube), public :: nanotube1
        !!<member name="nanotube1">One nanotube in the Bundle.</member>
        type(SingleWallNanotube), public :: nanotube2
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
    contains
        procedure, public :: get_lm1
        procedure, public :: get_lm2
        procedure, public :: get_gm1
        procedure, public :: get_gm2
        procedure, public :: get_kBar
        procedure, public :: get_kPrimeBar
        procedure, public :: calc_kBase
        procedure, public :: set_gm2
    end type BundleOfDoubleNanotube

    interface BundleOfDoubleNanotube
        procedure :: constructBundleOfDoubleNanotube
    end interface BundleOfDoubleNanotube

contains
    !!<summary>Constructor of BundleOfDoubleNanotube type.</summary>
    function constructBundleOfDoubleNanotube(n_1, n_2) result(this)
        implicit none
    !!<member name="n_1">The chiral vector of tube 1 C = n_1 * a_1 + n_2 * a_2. Tube 2's is (m + n) * a_1 - n * a_2</member>
        integer, intent(in) :: n_1
    !!<member name="n_2">The chiral vector of tube 1 C = n_1 * a_1 + n_2 * a_2. Tube 2's is (m + n) * a_1 - n * a_2</member>
        integer, intent(in) :: n_2
        type(BundleOfDoubleNanotube) :: this
        !!<local name="tmp">Inverse of(R - R inverse).</local>
        doubleprecision, dimension(2,2) :: tmp
        ! G_i = (R^{-1} - R) b_i
        doubleprecision, dimension(2) :: b1, b2

        this%nanotube1 = SingleWallNanotube(n_1, n_2)
        this%nanotube2 = SingleWallNanotube(n_1 + n_2, -n_2)

        ! Rotate both of them to the zigzag direction.
        call this%nanotube1%rotate(-this%nanotube1%getTheta())
        call this%nanotube2%rotate(-this%nanotube2%getTheta())
        this%kBar = (this%nanotube1%get_k() + this%nanotube2%get_k()) / 2.0
        this%kPrimeBar = (this%nanotube1%get_kPrime() + this%nanotube2%get_kPrime()) / 2.0

        tmp = reshape((/ dble(0), dble(-1) / (dble(2) * sin(this%nanotube1%getTheta())), dble(1) / (dble(2) * sin(this%nanotube1%getTheta())), dble(0) /), shape(tmp))
        this%lm1 = matmul(tmp, a1)
        this%lm2 = matmul(tmp, a2)

        call reciprocalLatticeVector(a1, a2, b1, b2)
        this%gm1 = rotate2d(b1, -this%nanotube1%getTheta()) - rotate2d(b1, this%nanotube1%getTheta())
        this%gm2 = rotate2d(b2, -this%nanotube1%getTheta()) - rotate2d(b2, this%nanotube1%getTheta())
    end function constructBundleOfDoubleNanotube

    !!<summary>Return the superlattice vector L_1^M.</summary>
    function get_lm1(this) result(lm1)
        implicit none
        class(BundleOfDoubleNanotube), intent(in) :: this
        doubleprecision, dimension(2) :: lm1

        lm1 = this%lm1
    end function get_lm1

    !!<summary>Return the superlattice vector L_2^M.</summary>
    function get_lm2(this) result(lm2)
        implicit none
        class(BundleOfDoubleNanotube), intent(in) :: this
        doubleprecision, dimension(2):: lm2

        lm2 = this%lm2
    end function get_lm2

    !!<summary>Return the superlattice vector G_1^M.</summary>
    function get_gm1(this) result(gm1)
        implicit none
        class(BundleOfDoubleNanotube), intent(in) :: this
        doubleprecision, dimension(2) :: gm1

        gm1 = this%gm1
    end function get_gm1

    !!<summary>Return the superlattice vector G_2^M.</summary>
    function get_gm2(this) result(gm2)
        implicit none
        class(BundleOfDoubleNanotube), intent(in) :: this
        doubleprecision, dimension(2) :: gm2

        gm2 = this%gm2
    end function get_gm2

    !!<summary>Return the coordinate of \bar{K}</summary>
    function get_kBar(this) result(kBar)
        implicit none
        class(BundleOfDoubleNanotube), intent(in) :: this
        doubleprecision, dimension(2) :: kBar

        kBar = this%kBar
    end function get_kBar

    !!<summary>Return the coordinate of \bar{K}'</summary>
    function get_kPrimeBar(this) result(kPrimeBar)
        implicit none
        class(BundleOfDoubleNanotube), intent(in) :: this
        doubleprecision, dimension(2) :: kPrimeBar

        kPrimeBar = this%kPrimeBar
    end function get_kPrimeBar

    !!<summary>Give the single layer bases at descrete k points to construct Hamiltonian within the cutoff.</summary>
    function calc_kBase(this, shift, kMax) result(res)
        implicit none
        class(BundleOfDoubleNanotube) :: this
        !!<parameter name="shift">Shift from center.</parameter>
        doubleprecision, dimension(2), intent(in) :: shift
        !!<parameter name="kMax">The radius of the cutoff circle</parameter>
        real, intent(in) ::  kMax
        type(kBase) :: res
        !!<local name="g_norm">The norm of reciprocal lattice vector.<local>
        doubleprecision :: g_norm
        doubleprecision, parameter :: piOverThree = 4.0 / 3.0 * atan(1.0)
        !!<local name="d">The distance from center to the line passing center plus shift.</local>
        doubleprecision :: d
        !!<local name="chord_length">One half of the length of the chord passing center plus shift paralleled to G_M^2.</local>
        doubleprecision :: halfChordLength
        !!<local name="distance_to_chord_middle_point">The distance between center plus shift to the middle point of the chord.</local>
        doubleprecision :: distance_to_chord_middle_point

        real :: l
        doubleprecision :: delta_l
        integer :: upperbound, lowerbound
        integer :: i, j

        call res%init(shift, this%gm1, this%gm2)
        g_norm = modulus2d(this%gm1)
        l = kMax / sin(piOverThree)
        delta_l = shift(2) / sin(piOverThree)

        do i = -floor((l - delta_l) / g_norm), floor((l + delta_l) / g_norm)
            d = i * g_norm * sin(piOverThree) - shift(2)
            halfChordLength = sqrt(kMax**2 - d**2)
            lowerbound = -floor((halfChordLength + i * this%gm1(1) + shift(1)) / g_norm)
            upperbound = -ceiling((-halfChordLength + i * this%gm1(1) + shift(1)) / g_norm)
            if(upperbound > lowerbound) then
                call res%appendHashParameter(upperbound - lowerbound + 1)

                do j = lowerbound, upperbound
                    call res%appendReducedCoordinate((/i, j/))
                end do
            end if
        end do
    end function calc_kBase

    subroutine set_gm2(this, gm2)
        implicit none
        class(BundleOfDoubleNanotube) :: this
        doubleprecision, dimension(2), intent(in) :: gm2

        this%gm2 = gm2
    end function
end module ClassBundleOfDoubleNanotube
