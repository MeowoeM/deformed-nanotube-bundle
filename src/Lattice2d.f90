module Lattice2d
    public :: crossEz2d
    public :: ezCross2d
    public :: ReciprocalLatticeVector
    public :: rotate2d
    public :: modulus2d
contains
    !!<summary>Cross product, c = a \times e_z.</summary>
    function crossEz2d(a) result(c)
        implicit none
        !!<parameter name="a">a in c = a \times e_z.</parameter>
        doubleprecision, dimension(2), intent(in)  :: a
        doubleprecision, dimension(2)              :: c

        c = (/ a(2), -a(1) /)
    end function crossEz2d


!!<member name="kPrime">K' point in reciprocal space</member>   !!<summary>Cross product, c = e_z \times a.</summary>
    function ezCross2d(a) result(c)
        implicit none
        !!<parameter name="a">a in c =  e_z \times a.</parameter>
        doubleprecision, dimension(2), intent(in) :: a
        doubleprecision, dimension(2)             :: c

        c = (/ -a(2), a(1) /)
    end function ezCross2d

    !!<summary>Reciprocal lattice vectors.</summary>
    subroutine reciprocalLatticeVector(a1, a2, b1, b2)
        !!<parameter name="a1">doubleprecision lattice vector a1.</parameter>
        doubleprecision, dimension(2), intent(in)   :: a1
        !!<parameter name="a2">doubleprecision lattice vector a2.</parameter>
        doubleprecision, dimension(2), intent(in)   :: a2
        !!<parameter name="b1">Reciprocal lattice vector b1</parameter>
        doubleprecision, dimension(2), intent(out)  :: b1
        !!<parameter name="b2">Reciprocal lattice vector b2</parameter>
        doubleprecision, dimension(2), intent(out)  :: b2

        b1 = 8.0 * atan(1.) * crossEz2d(a2) / dot_product(a1, crossEz2d(a2))
        b2 = 8.0 * atan(1.) * ezCross2d(a1) / dot_product(a2, ezCross2d(a1))
    end subroutine ReciprocalLatticeVector

    !!<summary>2D rotation about the origin.</summary>
    function rotate2d(a, theta) result(a_prime)
        !!<parameter name="a">The original vector or point.</parameter>
        doubleprecision, dimension(2), intent(in) :: a
        !!<parameter name="theta">Angle of rotation.</parameter>
        doubleprecision, intent(in)               :: theta
        doubleprecision, dimension(2)	           :: a_prime

        a_prime = (/ a(1) * cos(theta) - a(2) * sin(theta), a(2) * cos(theta) + a(1) * sin(theta) /)
    end function rotate2d

    !!<summary>Calculate the modulus of a 2d vector.</summary>
    function modulus2d(vector) result(modulus)
        !!<parameter name="vector">2d vector.</parameter>
        doubleprecision, dimension(2), intent(in) :: vector
        doubleprecision :: modulus

        modulus = sqrt(vector(1)**2 + vector(2)**2)
    end function modulus2d
end module Lattice2d
