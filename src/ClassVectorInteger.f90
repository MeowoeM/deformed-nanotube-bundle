module ClassVectorInteger
    !!<member name="SCALING_CONSTANT">Constrols how much the capacity scales onece the length of the list reach its capacity.</member>
    real, parameter, private :: SCALING_CONSTANT = 1.2
    !!<member name="INITIAL_CAPACITY">Capacity right after initialization.</member>
    integer, parameter, private :: INITIAL_CAPACITY = 1024

    !!<summary>Mimicing the vector template in C++ STL but only with real type BECAUSE THERE IS NO TEMPLATE IN FORTRAN!</summary>
    type VectorInteger
        !!<member name="list">Data storage.</member>
        integer, dimension(:), allocatable, private :: list
        !!<member name="length">Number of the element(s) in the list.</member>
        integer, private :: length
        !!<member name="capacity">Number of the elements can be held in the current list.</member>
        integer, private :: capacity
        !!<member name="isinitializeVectorIntegerd">whether the vector is initializeVectorIntegerd</member>
        logical, private :: isinitializeVectorIntegerd = .false.
    contains
        procedure, public :: append
        procedure, public :: get
        procedure, public :: getLength
        procedure, public :: getSum
        procedure, public :: writeVector
        procedure, public :: toArray
        final :: deleteVector
    end type VectorInteger

    interface VectorInteger
        procedure :: initializeVectorInteger
        procedure :: constructFromArray1d
    end interface VectorInteger

contains
    !!<summary>Specify the dimension of arrays in this vector.</summary>
    function initializeVectorInteger() result(this)
        type(VectorInteger) :: this

        this%length = 0
        this% capacity = INITIAL_CAPACITY
        this%isinitializeVectorIntegerd = .true.
        allocate(this%list(INITIAL_CAPACITY))
    end function initializeVectorInteger

    !!<summary>Construct a 1d array vector from a 2d array.</summary>
    function constructFromArray1d(array1d) result(this)
        integer, dimension(:), intent(in) :: array1d
        type(VectorInteger) :: this

        !!<local name="sizeInfo">The size info of the 2d array</local>
        this%length = size(array1d)
        this%capacity = ceiling(this%length * SCALING_CONSTANT)
        this%isinitializeVectorIntegerd = .true.
        allocate(this%list(this%capacity))

        do i = 1, this%length
            this%list(i) = array1d(i)
        end do
    end function constructFromArray1d

    !!<summary>Destructor.</summary>
    subroutine deleteVector (this)
        type(VectorInteger) :: this

        if ( allocated(this%list) ) then
            deallocate(this%list)
        end if
    end subroutine deleteVector

    !!<summary>Append an element to the vector.</summary>
    subroutine append(this, element)
        class(VectorInteger) :: this
        integer, intent(in) :: element
        integer, dimension(:), allocatable :: listTmp

        if ( this%isinitializeVectorIntegerd .eqv. .false. ) then
            print *, __FILE__, __LINE__, "ERROR! Vector uninitialized."
            return
        end if

        ! if length of the vector exceeds its capacity, copy all of its content to a new allocated memory
        if ( this%length == this%capacity ) then
            this%capacity = ceiling(this%capacity * SCALING_CONSTANT)
            allocate(listTmp(this%capacity))

            do i = 1, this%length
                listTmp(i) = this%list(i)
            end do

            deallocate(this%list)
            call move_alloc(listTmp, this%list)
        end if

        this%length = this%length + 1
        this%list(this%length) = element
    end subroutine append

    !!<summary>Get the element.</summary>
    function get(this, index) result(element)
        class(VectorInteger) :: this
        integer, intent(in) :: index
        integer :: element

        element = this%list(index)
    end function get

    !!<summary>Get the length of the vector</summary>
    function getLength(this) result(length)
        class(VectorInteger) :: this
        integer :: length

        length = this%length
    end function getLength

    !!<summary>Return the sum from the ith element in the vector to the jth</summary>
    function getSum(this, i, j) result(s)
        class(VectorInteger) :: this
        integer, intent(in) :: i
        integer, intent(in) :: j

        ! result
        integer :: s

        s = sum(this%list(i : j))
    end function

    subroutine writeVector(this, outUnit)
        class(VectorInteger) :: this
        integer :: outUnit

        integer, dimension(:), allocatable :: v

        allocate(v(this%length))
        v = this%list(1 : this%length)
        write (outUnit) v
        deallocate(v)
    end subroutine

    subroutine toArray(this, array)
        class(VectorInteger) :: this
        integer, allocatable, dimension(:) :: array

        allocate(array(this%length))
        array = this%list(1:this%length)
    end subroutine
end module ClassVectorInteger
