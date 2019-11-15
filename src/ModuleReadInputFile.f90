program TestLattice2d
    use Lattice2d
    use LatticeVector
    implicit none
    !real, dimension(2) :: a1
    !real, dimension(2) :: a2
    doubleprecision, dimension(2) :: b1
    doubleprecision, dimension(2) :: b2
    doubleprecision, dimension(2) :: c

    print *,a1
    c = crossEz2d(a1)
    print *,c

    print *,a1
    c = ezCross2d(a1)
    print *,c

    call ReciprocalLatticeVector(a1, a2, b1, b2)
    print *,b1
    print *,b2

end program TestLattice2d
