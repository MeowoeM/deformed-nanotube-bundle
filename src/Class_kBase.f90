#define DIMENSIONALITY 2
!!<summary>Points in K space(relative to K or K' points) based on which the Hamiltonian matrix is to construct.</summary>
module Class_kBase
    use ClassVectorIntegerArray1d
    use ClassVectorInteger

    !!<summary>Points in K space based on which the Hamiltonian matrix is to construct.</summary>
    type kBase
        !!<member name="reducedCoordinate">A pair of integers(i, j) denoting a k point coordinate by shift + i * G_1^M - j * G_2^M</member>
        type(VectorIntegerArray1d), private :: reducedCoordinate
        !!<member name="center">Center of the cutoff circle.</member>
        doubleprecision, dimension(DIMENSIONALITY), private :: center
        !!<member name="shift">Shift from center.</member>
        doubleprecision, dimension(DIMENSIONALITY), private :: shift
        !!<member name="xStep">x step of the k base grid.</member>
        doubleprecision, private :: xStep
        !!<member name="yStep">y step of the k base grid.</member>
        doubleprecision, private :: yStep
        !!<member name="cutoffRadius">The radius of the cutoff circle.</member>
        doubleprecision, private :: cutoffRadius
        !!<member name="xMinReduced">The minimum x reduced coordinate.</member>
        integer, private :: xMinReduced
        !!<member name="xMAxReduced">The maximum x reduced coordinate.</member>
        integer, private :: xMaxReduced
        !!<member name="cuttingLineNumber">Number of cutting lines within the cutoff circle.</member>
        integer, private :: cuttingLineNumber
        !!<member name="yUpperBoundList">The list of the upper bound of the y reduced coordinates on each cutting line.</member>
        integer, dimension(:), allocatable, private :: yUpperBoundList
        !!<member name="yLowerBoundList">The list of the lower bound of the y reduced coordinates on each cutting line.</member>
        integer, dimension(:), allocatable, private :: yLowerBoundList
        !!<member name="kPointNumberList">The number of k points on each cutting line.</member>
        integer, dimension(:), allocatable, private :: kPointNumberList
    contains
        final :: delete_kBase
        procedure, public :: initialize_kBase
        procedure, public :: appendReducedCoordinate
        procedure, public :: getReducedCoordinate
        procedure, public :: getCoordinate
        procedure, public :: getCoordinateX
        procedure, public :: getCoordinateY
        procedure, public :: getIndex
        procedure, public :: getIndex2
        procedure, public :: getQuantity
        procedure, public :: kPointInRangeX
        procedure, public :: get_xStep
        procedure, public :: get_yStep
        procedure, public :: getShift
        procedure, public :: getCenter
        procedure, public :: get_kPointNumberList
        procedure, public :: get_xMaxReduced
        procedure, public :: get_xMinReduced
        procedure, public :: getCutoffRadius
        procedure, public :: cuttingLineCoordinateInRange
        procedure, public :: setShift
    end type kBase
contains
    !!<summary>Destructor.</summary>
    subroutine delete_kBase(this)
    type(kBase) :: this

    if(allocated(this%yUpperBoundList)) then
        deallocate(this%yUpperBoundList)
    end if

    if(allocated(this%yLowerBoundList)) then
        deallocate(this%yLowerBoundList)
    end if

    if(allocated(this%kPointNumberList)) then
        deallocate(this%kPointNumberList)
    end if
    end subroutine

    function new_kBase(center, shift, xStep, yStep, cutoffRadius) result(this)
        type(kBase) :: this
        !!<parameter name="center">Center of the cutoff circle.</center>
        doubleprecision, dimension(DIMENSIONALITY), intent(in) :: center
        !!<parameter name="shift">Shift from the center.</parameter>
        doubleprecision, dimension(DIMENSIONALITY), intent(in) :: shift
        !!<parameter name="gm1">x step of the k base grid.</parameter>
        doubleprecision, intent(in) :: xStep
        !!<parameter name="gm2">y step of the k base grid.</parameter>
        doubleprecision, intent(in) :: yStep
        !!<parameter name="cutoffRadius">The radius of the cutoff circle.</parameter>
        doubleprecision, intent(in) :: cutoffRadius

        call this%initialize_kBase(center, shift, xStep, yStep, cutoffRadius)
    end function

    !!<summary>Initialize, reducedCoordinate,the integer array vector, and the center of the cutoff circle as well as the shift from it.</summary>
    subroutine initialize_kBase(this, center, shift, xStep, yStep, cutoffRadius)
        class(kBase) :: this
        !!<parameter name="center">Center of the cutoff circle.</center>
        doubleprecision, dimension(DIMENSIONALITY), intent(in) :: center
        !!<parameter name="shift">Shift from the center.</parameter>
        doubleprecision, dimension(DIMENSIONALITY), intent(in) :: shift
        !!<parameter name="gm1">x step of the k base grid.</parameter>
        doubleprecision, intent(in) :: xStep
        !!<parameter name="gm2">y step of the k base grid.</parameter>
        doubleprecision, intent(in) :: yStep
        !!<parameter name="cutoffRadius">The radius of the cutoff circle.</parameter>
        doubleprecision, intent(in) :: cutoffRadius

        ! tmp variable
        !!<local name="halfChordLength"><\local>
        doubleprecision :: halfChordLength
        integer :: counter

        this%reducedCoordinate = VectorIntegerArray1d(DIMENSIONALITY)
        this%center = center
        this%shift = shift
        this%xStep = xStep
        this%yStep = yStep
        this%cutoffRadius = cutoffRadius
        this%xMinReduced = -floor((cutoffRadius + shift(1)) / xStep)
        this%xMaxReduced = floor((cutoffRadius - shift(1)) / xStep)
        this%cuttingLineNumber = this%xMaxReduced - this%xMinReduced + 1

        allocate(this%yUpperBoundList(this%cuttingLineNumber))
        allocate(this%yLowerBoundList(this%cuttingLineNumber))
        allocate(this%kPointNumberList(this%cuttingLineNumber))

        counter = 0
        do i = this%xMinReduced, this%xMaxReduced
            counter = counter + 1
            halfChordLength = sqrt(cutoffRadius ** 2 - (shift(1) + i * xStep) ** 2)
            this%yLowerBoundList(counter) = -floor((halfChordLength) / yStep)
            this%yUpperBoundList(counter) = floor((halfChordLength) / yStep)
            if(this%yUpperBoundList(counter) >= this%yLowerBoundList(counter)) then
                do j = this%yLowerBoundList(counter), this%yUpperBoundList(counter)
                    call this%appendReducedCoordinate((/i, j/))
                end do
            end if
        end do
        this%kPointNumberList = this%yUpperBoundList - this%yLowerBoundList + 1
    end subroutine

    !!<summary>Add a k base.</summary>
    subroutine appendReducedCoordinate(this, reducedCoordinate)
        class(kBase) :: this
        integer, dimension(DIMENSIONALITY), intent(in) :: reducedCoordinate

        call this%reducedCoordinate%append(reducedCoordinate)
    end subroutine appendReducedCoordinate

    !!<summary>Return the reduced coordinate of one k base.</summary>
    function getReducedCoordinate(this, i) result(reducedCoordinate)
        class(kBase) :: this
        !!<parameter name="i">Index.<parameter>
        integer, intent(in) :: i

        ! result
        integer, dimension(DIMENSIONALITY) :: reducedCoordinate

        reducedCoordinate = this%reducedCoordinate%get(i)
    end function getReducedCoordinate

    !!<summary>Return the coordinate of one k base.<summary>
    function getCoordinate(this, i) result(coordinate)
        class(kBase) :: this
        !!<parameter name="i">Index.<parameter>
        integer, intent(in) :: i

        ! result
        doubleprecision, dimension(DIMENSIONALITY) :: coordinate

        !!<local></local>
        integer, dimension(DIMENSIONALITY) :: reducedCoordinate

        reducedCoordinate = this%reducedCoordinate%get(i)
        coordinate = reducedCoordinate * (/this%xStep, this%yStep/) + this%shift + this%center
    end function

    function getCoordinateX(this, xReducedCoordinate) result(coordinateX)
        class(kBase) :: this
        !!<parameter name="xReducedCoordinate"><parameter>
        integer, intent(in) :: xReducedCoordinate

        ! result
        doubleprecision :: coordinateX

        coordinateX = xReducedCoordinate * this%xStep + this%shift(1) + this%center(1)
     end function

     function getCoordinateY(this, yReducedCoordinate) result(coordinateY)
        class(kBase) :: this
        !!<parameter name="yReducedCoordinate"><parameter>
        integer, intent(in) :: yReducedCoordinate

        ! result
        doubleprecision :: coordinateY

        coordinateY = yReducedCoordinate * this%yStep + this%shift(2) + this%center(2)
     end function

    !!<summary>Given the reduced coordinate, return its index in O(1) time.</summary>
    function getIndex(this, reducedCoordinate) result(i)
        class(kBase) :: this
        integer, dimension(DIMENSIONALITY), intent(in) :: reducedCoordinate

        ! result
        integer :: i

        ! local
        integer, dimension(DIMENSIONALITY) :: tmpReducedCoordinate
        integer :: tmpIndex
        integer :: jMax
        integer :: jMin

        ! Return zero when out of the cutoff circle
        !tmpIndex = this%reducedCoordinate%getLength()
        !tmpReducedCoordinate = this%reducedCoordinate%get(tmpIndex)
        !iMax = tmpReducedCoordinate(1)
        !tmpIndex = 1
        !tmpReducedCoordinate = this%reducedCoordinate%get(tmpIndex)
        !iMin = tmpReducedCoordinate(1)

        if((reducedCoordinate(1) < this%xMinReduced) .or. (reducedCoordinate(1) > this%xMaxReduced)) then
            i = 0
            return
        end if

        tmpIndex = 0
        do i = 1, (reducedCoordinate(1) - this%xMinReduced + 1)
            tmpIndex = tmpIndex + this%kPointNumberList(i)
        end do
        tmpReducedCoordinate = this%reducedCoordinate%get(tmpIndex)
        jMax = tmpReducedCoordinate(2)
        jMin = jMax - this%kPointNumberList(reducedCoordinate(1) - this%xMinReduced + 1) + 1

        if ((reducedCoordinate(2) < jMin) .or. (reducedCoordinate(2) > jMax)) then
            i = 0
            return
        end if
        i = tmpIndex - (jMax - reducedCoordinate(2))
    end function

    function getIndex2(this, coordinate) result(i)
        class(kBase) :: this
        doubleprecision, dimension(DIMENSIONALITY), intent(in) :: coordinate

        ! result
        integer :: i

        ! local
        integer, dimension(DIMENSIONALITY) :: tmpReducedCoordinate, reducedCoordinate
        integer :: tmpIndex
        integer :: jMax
        integer :: jMin

        ! Return zero when out of the cutoff circle
        !tmpIndex = this%reducedCoordinate%getLength()
        !tmpReducedCoordinate = this%reducedCoordinate%get(tmpIndex)
        !iMax = tmpReducedCoordinate(1)
        !tmpIndex = 1
        !tmpReducedCoordinate = this%reducedCoordinate%get(tmpIndex)
        !iMin = tmpReducedCoordinate(1)

        reducedCoordinate(1) = floor(dot_product((/1.0, 0.0/), coordinate - this%shift - this%center) / this%xStep + 0.5)
        reducedCoordinate(2) = floor(dot_product((/0.0, 1.0/), coordinate - this%shift - this%center) / this%yStep + 0.5)
!        print *, __FILE__, __LINE__, reducedCoordinate

        if((reducedCoordinate(1) < this%xMinReduced) .or. (reducedCoordinate(1) > this%xMaxReduced)) then
            i = 0
            return
        end if

        tmpIndex = 0
        do i = 1, (reducedCoordinate(1) - this%xMinReduced + 1)
            tmpIndex = tmpIndex + this%kPointNumberList(i)
        end do
        tmpReducedCoordinate = this%reducedCoordinate%get(tmpIndex)
        jMax = tmpReducedCoordinate(2)
        jMin = jMax - this%kPointNumberList(reducedCoordinate(1) - this%xMinReduced + 1) + 1
        if ((reducedCoordinate(2) < jMin) .or. (reducedCoordinate(2) > jMax)) then
            i = 0
            return
        end if
        i = tmpIndex - (jMax - reducedCoordinate(2))
    end function

    !!<summary>Return the quantity of bases</summary>
    function getQuantity(this) result(quantity)
        class(kBase) :: this

        ! result
        integer :: quantity

        quantity = this%reducedCoordinate%getLength()
    end function

    !!<summary>Return the x points within given range. The range is described by a center and a x range.</summary>
    function kPointInRangeX(this, xCenter, xRange) result(xReducedCoordinate)
        class(kBase) :: this
        !!<parameter name="center">Center of the range<\parameter>
        doubleprecision, intent(in) :: xCenter
        !!<parameter name="xRange">Choosing k points in the range of center plus minus rangeX.</parameter>
        doubleprecision, intent(in) :: xRange

        type(VectorInteger) :: xReducedCoordinate

        ! local variables
        integer :: xMin, xMax

        xMin = ceiling((xCenter - xRange - this%shift(1) - this%center(1)) / this%xStep)
        xMin = max(xMin, this%xMinReduced)
        xMax = floor((xCenter + xRange - this%shift(1) - this%center(1)) / this%xStep)
        xMax = min(xMax, this%xMaxReduced)
        xReducedCoordinate = VectorInteger()
        do i = xMin, xMax
            call xReducedCoordinate%append(i)
        end do
    end function

    function get_xStep(this) result(xStep)
        class(kBase) :: this
        doubleprecision :: xStep3

        xStep = this%xStep
    end function

    function get_yStep(this) result(yStep)
        class(kBase) :: this
        doubleprecision :: yStep

        yStep = this%yStep
    end function

    function getShift(this) result(shift)
        class(kBase) :: this
        doubleprecision, dimension(DIMENSIONALITY) :: shift

        shift = this%shift
    end function

    function getCenter(this) result(center)
        class(kBase) :: this
        doubleprecision, dimension(DIMENSIONALITY) :: center

        center = this%center
    end function

    !!<summary>Returning a list of k point number on each cutting line.</summary>
    function get_kPointNumberList(this) result(kPointNumberList)
        class(kBase) :: this
        type(VectorInteger) :: kPointNumberList

        kPointNumberList = VectorInteger(this%kPointNumberList(1 : this%cuttingLineNumber))
    end function

    function get_xMinReduced(this) result(xMinReduced)
        class(kBase) :: this
        integer :: xMinReduced

        xMinReduced = this%xMinReduced
    end function

    function getCutoffRadius(this) result(cutoffRadius)
        class(kBase) :: this
        doubleprecision :: cutoffRadius

        cutoffRadius = this%cutoffRadius
    end function

    function get_xMaxReduced(this) result(xMaxReduced)
        class(kBase) :: this
        integer :: xMaxReduced

        xMaxReduced = this%xMaxReduced
    end function

    function cuttingLineCoordinateInRange(this, xBound) result(i)
        class(kBase) :: this
        doubleprecision, intent(in) :: xBound
        integer, dimension(2) :: i

        i(1) = max(ceiling((-xBound - this%shift(1)) / this%xStep), this%xMinReduced)
        i(2) = min(floor((xBound - this%shift(1)) / this%xStep), this%xMaxReduced)
    end function

    subroutine setShift(this, shift)
        class(kBase) :: this
        doubleprecision, dimension(DIMENSIONALITY), intent(in) :: shift

        this%shift = shift
    end subroutine
end module Class_kBase






