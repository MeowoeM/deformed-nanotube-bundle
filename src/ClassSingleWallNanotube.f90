module ClassSingleWallNanotube
    use Lattice2d
    use LatticeVector

    !!<summary>Single wall nanotube class.</summary>
    type SingleWallNanotube
        !!<member name="n1">The chiral vector of the upper layer C = n1 * a_1 + n2 * a_2.</member>
        integer, private :: m
        !!<member name="n2">The chiral vector of the upper layer C = n1 * a_1 + n2 * a_2.</member>
        integer, private :: n
        !!<member name="orientation">When the nanotube is placed parallel to z direction, the direction to which the end point of the chiral vector points.</member>
        doubleprecision, private :: orientation
        !!<member name="shiftZ">Relative dislocation in the z direction.</member>
        doubleprecision, private :: shiftZ
        !!<member name="theta">The angle from the zigzag direction to the chiral vector.</member>
        doubleprecision, private :: theta
        !!<member name="a1">Lattice vector.</member>
        doubleprecision, dimension(2), private :: a1
        !!<member name="a2">Lattice vector.</member>
        doubleprecision, dimension(2), private :: a2
        !!<member name="b1">Reciprocal lattice vector.</member>
        doubleprecision, dimension(2), private :: b1
        !!<member name="b2">Reciprocal lattice vector.</member>
        doubleprecision, dimension(2), private :: b2
        !!<member name="k">K point.</member>
        doubleprecision, dimension(2), private :: k
        !!<member name="kPrime">K' point.</member>
        doubleprecision, dimension(2), private :: kPrime
    contains
        procedure, public :: initializeSingleWallNanotube
        procedure, public :: rotate
        procedure, public :: getTheta
        procedure, public :: get_k
        procedure, public :: get_kPrime
        procedure, public :: get_m
        procedure, public :: get_n
        procedure, public :: getOrientation
        procedure, public :: getChiralVectorLength
        procedure, public :: getShiftZ
    end type SingleWallNanotube
contains
    function newSingleWallNanotube(m, n, orientation) result(this)
    type(SingleWallNanotube) :: this
        !!<parameter name="n1">The chiral vector of the upper layer C = n1 * a_1 + n2 * a_2.</parameter>
        integer, intent(in) :: m
        !!<parameter name="n2">The chiral vector of the upper layer C = n1 * a_1 + n2 * a_2.</parameter>
        integer, intent(in) :: n
        !!<parameter name="orientation">When the nanotube is placed parallel to z direction, the direction to which the end point of the chiral vector points.</parameter>
        doubleprecision :: orientation

        call this%initializeSingleWallNanotube(m, n, orientation)
    end function

    !!<summary>Constructor of SingleWallNanotube type.</summary>
    subroutine initializeSingleWallNanotube(this, m, n, orientation)
        class(SingleWallNanotube) :: this
        !!<parameter name="n1">The chiral vector of the upper layer C = n1 * a_1 + n2 * a_2.</parameter>
        integer, intent(in) :: m
        !!<parameter name="n2">The chiral vector of the upper layer C = n1 * a_1 + n2 * a_2.</parameter>
        integer, intent(in) :: n
        !!<parameter name="orientation">When the nanotube is placed parallel to z direction, the direction to which the end point of the chiral vector points.</parameter>
        doubleprecision :: orientation

        this%m = m
        this%n = n
        this%orientation = orientation

        this%a1 = a1
        this%a2 = a2

        this%k = k
        this%kPrime = kPrime

        call ReciprocalLatticeVector(a1, a2, this%b1, this%b2)

        this%theta = atan((dble(m) * a1(2) + dble(n) * a2(2)) / (dble(m) * a1(1) + dble(n) * a2(1)))

        call this%rotate(-this%theta)
    end subroutine initializeSingleWallNanotube

    !!<summary>Rotate the layer by angle theta.</summary>
    subroutine rotate(this, theta)
      class(SingleWallNanotube) :: this
      !!<parameter name=$"1:para_name">Rotation angle.</parameter>
      doubleprecision, intent(in) :: theta

      this%a1 = rotate2d(this%a1, theta)
      this%a2 = rotate2d(this%a2, theta)
      this%b1 = rotate2d(this%b1, theta)
      this%b2 = rotate2d(this%b2, theta)
      this%k = rotate2d(this%k, theta)
      this%kPrime = rotate2d(this%kPrime, theta)
    end subroutine rotate

    !!<summary>Get the angle theta from zigzag direction to the chiral vector.</summary>
    function getTheta(this) result(theta)
        class(SingleWallNanotube) :: this
        doubleprecision :: theta

        theta = this%theta
    end function getTheta

    !!<summary>Get K point coordinate.</summary>
    function get_k(this) result(k)
        class(SingleWallNanotube) :: this
        doubleprecision, dimension(2) :: k

        k = this%k
    end function get_k

    !!<summary>Get K' point coordinate.</summary>
    function get_kPrime(this) result(kPrime)
        class(SingleWallNanotube) :: this
        doubleprecision, dimension(2) :: kPrime

        kPrime = this%kPrime
    end function get_kPrime

    !!<summary>Get n1 of the chiral vector C = n1 * a_1 + n2 * a_2</summary>
    function get_m(this) result(m)
        class(SingleWallNanotube) :: this
        doubleprecision :: m

        m = this%m
    end function get_m

    !!<summary>Get n2 of the chiral vector C = n1 * a_1 + n2 * a_2</summary>
    function get_n(this) result(n)
        class(SingleWallNanotube) :: this
        doubleprecision :: n2

        n = this%n
    end function get_n

    !!<summary>Get n2 of the chiral vector C = n1 * a_1 + n2 * a_2</summary>
    function getOrientation(this) result(orientation)
        class(SingleWallNanotube) :: this
        doubleprecision :: orientation

        orientation = this%orientation
    end function getOrientation

    !!<summary>Get the modulus of the chiral vector.</summary>
    function getChiralVectorLength(this) result(modulus)
        class(SingleWallNanotube), intent(in) :: this
        doubleprecision :: modulus

        modulus = modulus2d(this%m * a1 + this%n * a2)
    end function getChiralVectorLength

    function getShiftZ(this) result(shiftZ)
        class(SingleWallNanotube), intent(in) :: this
        doubleprecision :: shiftZ

        shiftZ = this%shiftZ
    end function
end module ClassSingleWallNanotube



