!!<summary>Mimicing the vector template in C++ STL but only with integer array type BECAUSE THERE IS NO TEMPLATE IN FORTRAN!</summary>
module ClassVectorIntegerArray1d
    !!<member name="SCALING_CONSTANT">Controls how much the capacity scales once the length of the list reach its capacity.</member>
    doubleprecision, parameter, private :: SCALING_CONSTANT = 1.5
    !!<member name="INITIAL_CAPACITY">Capacity right after initialization.</member>
    integer, parameter, private :: INITIAL_CAPACITY = 1024

    !!<summary>Mimicing the vector template in C++ STL but only with doubleprecision type BECAUSE THERE IS NO TEMPLATE IN FORTRAN!</summary>
    type VectorIntegerArray1d
        !!<member name="list">Data storage.</member>
        integer, dimension(:, :), allocatable, private :: list
        !!<member name="arraySize">The size of arrays in this vector.</member>
        integer, private :: arraySize
        !!<member name="length">Number of the element(s) in the list.</member>
        integer, private :: length
        !!<member name="capacity">Number of the elements can be held in the current list.</member>
        integer, private :: capacity
        !!<member name="isinitializeVectorIntegerArray1dd">whether the vector is initializeVectorIntegerArray1dd</member>
        logical, private :: isinitializeVectorIntegerArray1dd = .false.
    contains
        procedure, public :: append
        procedure, public :: get
        procedure, public :: getLength
        final :: deleteVector
    end type VectorIntegerArray1d

    interface VectorIntegerArray1d
        module procedure :: initializeVectorIntegerArray1d
        module procedure :: constructFromArray2d
    end interface VectorIntegerArray1d

contains
    !!<summary>Specify the dimension of arrays in this vector.</summary>
    function initializeVectorIntegerArray1d(arraySize) result(this)
        type(VectorIntegerArray1d) :: this
        integer, intent(in) :: arraySize

        this%arraySize = arraySize
        this%length = 0
        this% capacity = INITIAL_CAPACITY
        this%isinitializeVectorIntegerArray1dd = .true.
        allocate(this%list(INITIAL_CAPACITY, this%arraySize))
    end function initializeVectorIntegerArray1d

    !!<summary>Construct a 1d array vector from a 2d array.</summary>
    function constructFromArray2d(array2d) result(this)
        integer, dimension(:, :), intent(in) :: array2d
        type(VectorIntegerArray1d) :: this

        !!<local name="sizeInfo">The size info of the 2d array</local>
        integer, dimension(2) :: sizeInfo
        sizeInfo = shape(array2d)
        this%arraySize = sizeInfo(2)
        this%length = sizeInfo(1)
        this%capacity = ceiling(sizeInfo(1) * SCALING_CONSTANT)
        this%isinitializeVectorIntegerArray1dd = .true.
        allocate(this%list(this%capacity, this%arraySize))

        do i = 1, this%length
            this%list(i, :) = array2d(i, :)
        end do
    end function constructFromArray2d

    !!<summary>Destructor.</summary>
    subroutine deleteVector (this)
        type(VectorIntegerArray1d) :: this

        if ( allocated(this%list) ) then
            deallocate(this%list)
        end if
    end subroutine deleteVector

    !!<summary>Append an element to the vector.</summary>delete_vector
    subroutine append(this, element)
        class(VectorIntegerArray1d) :: this
        integer, dimension(:), intent(in) :: element
        integer, dimension(:, :), allocatable :: listTmp

        if ( this%isinitializeVectorIntegerArray1dd .eqv. .false. ) then
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
        class(VectorIntegerArray1d) :: this
        integer, intent(in) :: index
        doubleprecision, dimension(this%arraySize) :: element

        element = this%list(index, 1:this%arraySize)
    end function get

    !!<summary>Get the length of the vector</summary>
    function getLength(this) result(length)
        class(VectorIntegerArray1d) :: this
        integer :: length

        length = this%length
    end function getLength
end module ClassVectorIntegerArray1d
