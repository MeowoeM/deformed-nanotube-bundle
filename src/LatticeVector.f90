!!<summary>Definition of Lattice Vectors.</summary>
module LatticeVector
    implicit none
    !!<member name="a">Lattice constant.</member>
    doubleprecision, parameter :: a = 0.24612
    !!<member name="a1">Lattice vector</member>
    doubleprecision, dimension(2), parameter :: a1 = (/ 1.0, 0.0 /) * a
    !!<member name="a2">Lattice vector</member>
    doubleprecision, dimension(2), parameter :: a2 = (/ 1.0 / 2.0, sqrt(3.0) / 2.0 /) * a
    !!<member name="k">K point in reciprocal space</member>
    doubleprecision, dimension(2), parameter :: k = (/ -2.0 / 3.0, 0.0 /) * 8.0 * atan(1.0) / a
    !!<member name="kPrime">K' point in reciprocal space</member>
    doubleprecision, dimension(2), parameter :: kPrime = (/ 2.0 / 3.0, 0.0 /) * 8.0 * atan(1.0) / a
end module LatticeVector
