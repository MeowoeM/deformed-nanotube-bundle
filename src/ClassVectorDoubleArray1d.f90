!!<summary>Mimicing the vector template in C++ STL but only with double array type BECAUSE THERE IS NO TEMPLATE IN FORTRAN!</summary>
module ClassVectorDoubleArray1d
    !!<member name="SCALING_CONSTANT">Controls how much the capacity scales once the length of the list reach its capacity.</member>
    doubleprecision, parameter, private :: SCALING_CONSTANT = 1.5
    !!<member name="INITIAL_CAPACITY">Capacity right after initialization.</member>
    integer, parameter, private :: INITIAL_CAPACITY = 1024

    !!<summary>Mimicing the vector template in C++ STL but only with doubleprecision type BECAUSE THERE IS NO TEMPLATE IN FORTRAN!</summary>
    type VectorDoubleArray1d
        !!<member name="list">Data storage.</member>
        doubleprecision, dimension(:, :), allocatable, private :: list
        !!<member name="arraySize">The size of arrays in this vector.</member>
        integer, private :: arraySize
        !!<member name="length">Number of the element(s) in the list.</member>
        integer, private :: length
        !!<member name="capacity">Number of the elements can be held in the current list.</member>
        integer, private :: capacity
        !!<member name="isinitializeVectorDoubleArray1dd">whether the vector is initializeVectorDoubleArray1dd</member>
        logical, private :: isinitializeVectorDoubleArray1dd = .false.
    contains
        procedure, public :: append
        procedure, public :: get
        procedure, public :: getLength
        procedure, public :: writeVector
        procedure, public :: initializeVectorDoubleArray1d
        procedure, public :: constructFromArray2d
        final :: deleteVector
    end type VectorDoubleArray1d

contains
    function newVectorDoubleArray1d(arraySize) result(this)
        type(VectorDoubleArray1d) :: this
        integer, intent(in) :: arraySize

        call this%initializeVectorDoubleArray1d(arraySize)
    end function

    !!<summary>Specify the dimension of arrays in this vector.</summary>
    subroutine initializeVectorDoubleArray1d(this, arraySize)
        class(VectorDoubleArray1d) :: this
        integer, intent(in) :: arraySize

        this%arraySize = arraySize
        this%length = 0
        this% capacity = INITIAL_CAPACITY
        this%isinitializeVectorDoubleArray1dd = .true.
        if ( allocated(this%list) ) then
            deallocate(this%list)
        end if
        allocate(this%list(INITIAL_CAPACITY, this%arraySize))
    end subroutine initializeVectorDoubleArray1d

    !!<summary>Construct a 1d array vector from a 2d array.</summary>
    subroutine constructFromArray2d(this, array2d)
        doubleprecision, dimension(:, :), intent(in) :: array2d
        class(VectorDoubleArray1d) :: this

        !!<local name="sizeInfo">The size info of the 2d array</local>
        integer, dimension(2) :: sizeInfo
        sizeInfo = shape(array2d)
        this%arraySize = sizeInfo(2)
        this%length = sizeInfo(1)
        this%capacity = ceiling(sizeInfo(1) * SCALING_CONSTANT)
        this%isinitializeVectorDoubleArray1dd = .true.
        if ( allocated(this%list) ) then
            deallocate(this%list)
        end if
        allocate(this%list(this%capacity, this%arraySize))

        do i = 1, this%length
            this%list(i, :) = array2d(i, :)
        end do
    end subroutine constructFromArray2d

    !!<summary>Destructor.</summary>
    subroutine deleteVector (this)
        type(VectorDoubleArray1d) :: this

        if ( allocated(this%list) ) then
            deallocate(this%list)
        end if
    end subroutine deleteVector

    !!<summary>Append an element to the vector.</summary>delete_vector
    subroutine append(this, element)
        class(VectorDoubleArray1d) :: this
        doubleprecision, dimension(:), intent(in) :: element
        doubleprecision, dimension(:, :), allocatable :: listTmp

        if ( this%isinitializeVectorDoubleArray1dd .eqv. .false. ) then
            print *, __FILE__, __LINE__, "ERROR! Vector uninitialized."
            return
        end if

        if ( this%arraySize /= size(element) ) then
            print *, __FILE__, __LINE__, "ERROR! The size of the arry to append doesn't match this vector."
            return
        end if

        ! if length of the vector exceeds its capacity, copy all of delete_vectorits content to a new allocated memory
        if ( this%length == this%capacity ) then
            this%capacity = ceiling(this%capacity * SCALING_CONSTANT)
            allocate(listTmp(this%capacity, this%arraySize))

            do i = 1, this%length
                listTmp(i, :) = this%list(i, :)
            end do

            deallocate(this%list)
            call move_alloc(listTmp, this%list)
        end if

        this%length = this%length + 1
        this%list(this%length, :) = element
    end subroutine append

    !!<summary>Get the element.</summary>
    function get(this, index) result(element)
        class(VectorDoubleArray1d) :: this
        integer, intent(in) :: index
        doubleprecision, dimension(this%arraySize) :: element

        element = this%list(index, 1:this%arraySize)
    end function get

    !!<summary>Get the length of the vector</summary>
    function getLength(this) result(length)
        class(VectorDoubleArray1d) :: this
        integer :: length

        length = this%length
    end function getLength

    subroutine writeVector(this, outUnit)
        class(VectorDoubleArray1d) :: this
        integer :: outUnit

        doubleprecision, dimension(:, :), allocatable :: m

        allocate(m(this%length, 3))
        m = this%list(1 : this%length, :)
        write (outUnit) m
        deallocate(m)
    end subroutine
end module
