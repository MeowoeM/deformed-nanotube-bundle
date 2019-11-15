!!<summary>Calculate the spectral function for output from the eigenstates and eigenenergy of the Hamiltonian matrix.</summary>
module ClassSpectralFunction
    use Class_kBase
    use ClassVectorDouble
    use ClassVectorDoubleArray1d
    use ClassVectorInteger

    type VectorDoubleArray1dContainer
        class(VectorDoubleArray1d), allocatable :: subVector
    end type

    type SpectralFunction
        !!<parameter name="cuttingLineReducedCoordinateRange1">Reduced coordinate(x) range of cutting lines to output. tube 1</parameter>
        integer, dimension(2), private :: cuttingLineReducedCoordinateRange1
        !!<parameter name="cuttingLineReducedCoordinateRange1">Reduced coordinate(x) range of cutting lines to output. tube 2</parameter>
        integer, dimension(2), private :: cuttingLineReducedCoordinateRange2
        !!<parameter name="output1">Spectral function A of k_y and energy, to output. (/ k_y, energy, A(k_y, energy) /) tube 1<\parameter>
        type(VectorDoubleArray1dContainer), dimension(:), allocatable, private :: output1
        !!<parameter name="output2">Spectral function A of k_y and energy, to output. (/ k_y, energy, A(k_y, energy) /) tube 2<\parameter>
        type(VectorDoubleArray1dContainer), dimension(:), allocatable, private :: output2
        !!<parameter name="kPointNumberList1">A List of k point numbers on each cuttinglines. tube 1</parameter>
        type(VectorInteger), private :: kPointNumberList1
        !!<parameter name="kPointNumberList2">A List of k point numbers on each cuttinglines. tube 2</parameter>
        type(VectorInteger), private :: kPointNumberList2
        !!<parameter name="cuttingLineCatagory1">The k_x coordinates of cuttinglines to output. tube 1</parameter>
        type(VectorDouble), private :: cuttingLineCatagory1
        !!<parameter name="cuttingLineCatagory2">The k_x coordinates of cuttinglines to output. tube 2</parameter>
        type(VectorDouble), private :: cuttingLineCatagory2
    contains
        final :: deleteSpectralFunction
        procedure, public :: initSpectralFunction
        procedure, public :: extend
        procedure, public :: output
        procedure, public :: printCuttinglinePosition
        procedure, public :: updateSpectralFunction
    end type
contains
    subroutine initSpectralFunction(this, base1, base2, cuttingLineRange)
        class(SpectralFunction) :: this
        type(kBase), intent(in) :: base1, base2
        doubleprecision, intent(in) :: cuttingLineRange

        ! local variables
        integer :: tmpIndex
        doubleprecision, dimension(2) :: tmpCoordinate
        integer, dimension(2) :: tmpVectorInteger

        this%cuttingLineReducedCoordinateRange1 = base1%cuttingLineCoordinateInRange(cuttingLineRange)
        this%cuttingLineReducedCoordinateRange2 = base2%cuttingLineCoordinateInRange(cuttingLineRange)

        allocate(this%output1(this%cuttingLineReducedCoordinateRange1(2) - this%cuttingLineReducedCoordinateRange1(1) + 1))
        allocate(this%output2(this%cuttingLineReducedCoordinateRange2(2) - this%cuttingLineReducedCoordinateRange2(1) + 1))

        do i = 1, this%cuttingLineReducedCoordinateRange1(2) - this%cuttingLineReducedCoordinateRange1(1) + 1
            this%output1(i) = VectorDoubleArray1dContainer(subVector = newVectorDoubleArray1d(3))
        end do
        do i = 1, this%cuttingLineReducedCoordinateRange2(2) - this%cuttingLineReducedCoordinateRange2(1) + 1
            this%output2(i) = VectorDoubleArray1dContainer(subVector = newVectorDoubleArray1d(3))
        end do

        this%kPointNumberList1 = base1%get_kPointNumberList()
        this%kPointNumberList2 = base2%get_kPointNumberList()

        this%cuttingLineCatagory1 = VectorDouble()
        this%cuttingLineCatagory2 = VectorDouble()
        call this%updateSpectralFunction(base1, base2)
    end subroutine

    subroutine updateSpectralFunction(this, base1, base2)
        implicit none
        class(SpectralFunction) :: this
        type(kBase) :: base1, base2

        ! local variables
        integer :: tmpIndex
        doubleprecision, dimension(2) :: tmpCoordinate
        integer, dimension(2) :: tmpVectorInteger
        integer :: i

        tmpVectorInteger = this%cuttingLineReducedCoordinateRange1
        do i = 1, tmpVectorInteger(2) - tmpVectorInteger(1) + 1
            tmpIndex = this%kPointNumberList1%getSum(1, tmpVectorInteger(1) + i - 1 - base1%get_xMinReduced()) + 1
            tmpCoordinate = base1%getCoordinate(tmpIndex)
            call this%cuttingLineCatagory1%append(tmpCoordinate(1))
        end do
        tmpVectorInteger = this%cuttingLineReducedCoordinateRange2
        do i = 1, tmpVectorInteger(2) - tmpVectorInteger(1) + 1
            tmpIndex = this%kPointNumberList2%getSum(1, tmpVectorInteger(1) + i - 1 - base2%get_xMinReduced()) + 1
            tmpCoordinate = base2%getCoordinate(tmpIndex)
            call this%cuttingLineCatagory2%append(tmpCoordinate(1))
        end do
    end subroutine

    subroutine deleteSpectralFunction(this)
        type(SpectralFunction) :: this

        if(allocated(this%output1)) then
            deallocate(this%output1)
        end if

        if(allocated(this%output2)) then
            deallocate(this%output2)
        end if
    end subroutine

    !!<summary>Calculate spectral function from input eigenstates and their eigen energy. Add the result to the existing spectral function to enlarge its range.</summary>
    subroutine extend(this, base1, base2, eigenenergy, eigenvector, m)
        class(SpectralFunction) :: this
        type(kBase), intent(in) :: base1, base2
        doubleprecision, dimension(:), intent(in) :: eigenenergy
        complex(8), dimension(:, :), intent(in) :: eigenvector
        integer, intent(in) :: m

        ! local variables
        integer :: tmpInteger
        integer, dimension(2) :: tmpVectorInteger
        doubleprecision, dimension(2) :: tmpCoordinate
        doubleprecision :: tmp2, mm

        mm = 10

        do i = 1, m
            tmpVectorInteger = this%cuttingLineReducedCoordinateRange1
            do j = 1, tmpVectorInteger(2) - tmpVectorInteger(1) + 1
                tmpInteger = this%kPointNumberList1%getSum(1, tmpVectorInteger(1) + j - 1 - base1%get_xMinReduced()) + 1
                do kk = tmpInteger, tmpInteger + this%kPointNumberList1%get(tmpVectorInteger(1) - base1%get_xMinReduced() + j)
                    tmp2 = abs(eigenvector(kk                      , i)) ** 2 + &
                           abs(eigenvector(kk + base1%getQuantity(), i)) ** 2
                    tmpCoordinate = base1%getCoordinate(kk)
                    call this%output1(j)%subVector%append((/tmpCoordinate(2), eigenenergy(i), tmp2/))
                end do
            end do

            tmpVectorInteger = this%cuttingLineReducedCoordinateRange2
            do j = 1, tmpVectorInteger(2) - tmpVectorInteger(1) + 1
                tmpInteger = this%kPointNumberList2%getSum(1, tmpVectorInteger(1) + j - 1 - base2%get_xMinReduced()) + 1
                do kk = tmpInteger, tmpInteger + this%kPointNumberList2%get(tmpVectorInteger(1) - base2%get_xMinReduced() + j)
                    tmp2 = abs(eigenvector(kk + 2 * base1%getQuantity()                      , i)) ** 2 + &
                           abs(eigenvector(kk + 2 * base1%getQuantity() + base2%getQuantity(), i)) ** 2
                    tmpCoordinate = base2%getCoordinate(kk)
                    if(tmpCoordinate(2) < mm) then
                        mm = tmpCoordinate(2)
                    end if
                    call this%output2(j)%subVector%append((/tmpCoordinate(2), eigenenergy(i), tmp2/))
                end do
            end do
        end do
    end subroutine

    subroutine printCuttinglinePosition(this)
        class(SpectralFunction) :: this

        print *, __FILE__, __LINE__, 'cuttingline positions of tube 1:'
        do i = 1, this%cuttingLineCatagory1%getLength()
            print *, this%cuttingLineCatagory1%get(i)
        end do

        print *, __FILE__, __LINE__, 'cuttingline positions of tube 2:'
        do i = 1, this%cuttingLineCatagory2%getLength()
            print *, this%cuttingLineCatagory2%get(i)
        end do
    end subroutine

    subroutine output(this, outUnit)
        class(SpectralFunction) :: this
        integer, intent(in) :: outUnit

        ! local variable
        integer :: tmpInteger

        tmpInteger = this%cuttingLineCatagory1%getLength()
        write(outUnit) tmpInteger
        call this%cuttingLineCatagory1%writeVector(outUnit)
        do i = 1, tmpInteger
            write(outUnit) this%output1(i)%subVector%getLength()
            call this%output1(i)%subVector%writeVector(outUnit)
        end do
        tmpInteger = this%cuttingLineCatagory2%getLength()
        write(outUnit) tmpInteger
        call this%cuttingLineCatagory2%writeVector(outUnit)
        do i = 1, tmpInteger
            write(outUnit) this%output2(i)%subVector%getLength()
            call this%output2(i)%subVector%writeVector(outUnit)
        end do
    end subroutine
end module
