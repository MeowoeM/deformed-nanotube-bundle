module ClassVectorComplex8
    !!<member name="SCALING_CONSTANT">Constrols how much the capacity scales onece the length of the list reach its capacity.</member>
    real, parameter, private :: SCALING_CONSTANT = 1.2
    !!<member name="INITIAL_CAPACITY">Capacity right after initialization.</member>
    integer, parameter, private :: INITIAL_CAPACITY = 1024

    !!<summary>Mimicing the vector template in C++ STL but only with real type BECAUSE THERE IS NO TEMPLATE IN FORTRAN!</summary>
    type VectorComplex8
        !!<member name="list">Data storage.</member>
        complex(8), dimension(:), allocatable, private :: list
        !!<member name="length">Number of the element(s) in the list.</member>
        integer, private :: length
        !!<member name="capacity">Number of the elements can be held in the current list.</member>
        integer, private :: capacity
        !!<member name="isinitializeVectorComplex8d">whether the vector is initializeVectorComplex8d</member>
        logical, private :: isinitializeVectorComplex8d = .false.
    contains
        procedure, public :: append
        procedure, public :: get
        procedure, public :: getLength
        procedure, public :: getSum
        procedure, public :: writeVector
        procedure, public :: toArray
        final :: deleteVector
    end type VectorComplex8

    interface VectorComplex8
        procedure :: initializeVectorComplex8
        procedure :: constructFromArray1d
    end interface VectorComplex8

contains
    !!<summary>Specify the dimension of arrays in this vector.</summary>
    function initializeVectorComplex8() result(this)
        type(VectorComplex8) :: this

        this%length = 0
        this%capacity = INITIAL_CAPACITY
        this%isinitializeVectorComplex8d = .true.
        allocate(this%list(INITIAL_CAPACITY))
    end function initializeVectorComplex8

    !!<summary>Construct a 1d array vector from a 2d array.</summary>
    function constructFromArray1d(array1d) result(this)
        complex(8), dimension(:), intent(in) :: array1d
        type(VectorComplex8) :: this

        !!<local name="sizeInfo">The size info of the 2d array</local>
        this%length = size(array1d)
        this%capacity = ceiling(this%length * SCALING_CONSTANT)
        this%isinitializeVectorComplex8d = .true.
        allocate(this%list(this%capacity))

        do i = 1, this%length
            this%list(i) = array1d(i)
        end do
    end function constructFromArray1d

    !!<summary>Destructor.</summary>
    subroutine deleteVector (this)
        type(VectorComplex8) :: this

        if ( allocated(this%list) ) then
            deallocate(this%list)
        end if
    end subroutine deleteVector

    !!<summary>Append an element to the vector.</summary>
    subroutine append(this, element)
        class(VectorComplex8) :: this
        complex(8), intent(in) :: element
        complex(8), dimension(:), allocatable :: listTmp

        if ( this%isinitializeVectorComplex8d .eqv. .false. ) then
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
        class(VectorComplex8) :: this
        integer, intent(in) :: index
        real :: element

        element = this%list(index)
    end function get

    !!<summary>Get the length of the vector</summary>
    function getLength(this) result(length)
        class(VectorComplex8) :: this
        integer :: length

        length = this%length
    end function getLength

    !!<summary>Return the sum from the ith element in the vector to the jth</summary>
    function getSum(this, i, j) result(s)
        class(VectorComplex8) :: this
        integer, intent(in) :: i
        integer, intent(in) :: j

        ! result
        complex(8) :: s

        s = sum(this%list(i : j))
    end function

    subroutine writeVector(this, outUnit)
        class(VectorComplex8) :: this
        integer :: outUnit

        complex(8), dimension(:), allocatable :: v

        allocate(v(this%length))
        v = this%list(1 : this%length)
        write (outUnit) v
        deallocate(v)
    end subroutine

    subroutine toArray(this, array)
        class(VectorComplex8) :: this
        complex(8), allocatable, dimension(:) :: array

        allocate(array(this%length))
        array = this%list(1:this%length)
    end subroutine
end module ClassVectorComplex8


