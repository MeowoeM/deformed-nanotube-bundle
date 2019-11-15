!!<summary>Coupling strength decaying between two tubes.</summary>
module ClassInterlayerAttenuation
implicit none
!!<summary>Coupling strength decaying between two tubes.</summary>
    type InterlayerAttenuation
        !!<member name="n">Number of subintervals on x axis. </member>
        integer, private :: n
        !!<member name="w">The width unit for the attenuation function, preferrably the width of Gaussian for the original undeformed tubes</member>
        doubleprecision, private :: w
        !!<member name="nw">Length of flat region in the middle of attenuation fuction in unit of w</member>
        doubleprecision, private :: nw
        !!<member name="xMin">Range to carry out the Fourier transformation of attenuation function, -(nw + 3) by default.</member>
        doubleprecision, private :: xMin
        !!<member name="xMax">Range to carry out the Fourier transformation of attenuation function, (nw + 3) by default.</member>
        doubleprecision, private :: xMax
        !!<member name="kMin">Cutoff for Fourier transformed attenuation function.</member>
        doubleprecision, private :: kMin
        !!<member name="kMax">Cutoff for Fourier transformed attenuation function.</member>
        doubleprecision, private :: kMax
        !!<member name="nk">Number of intervals in k space. Dimension of the table for Fourier transformed attenuation function is (k + 1, 2).</member>
        integer, private :: nk
        !!<member name="kTable">Table for Fourier transformed attenuation function.</member>
        doubleprecision, dimension(:, :), allocatable, private :: kTable
    contains
        procedure, public :: init => initializeInterlayerAttenuation
        procedure, public :: attenuationFunction
        procedure, public :: attenuationFunctionFT
        procedure, public :: genterate_kTable
        final :: deleteInterlayerAttenuation
    end type

    interface InterlayerAttenuation
        procedure :: newInterlayerAttenuation
    end interface
contains
    !!<summary>Constructor.</summary>
    function newInterlayerAttenuation(n, nw, w) result(this)
        type(InterlayerAttenuation) :: this
        integer, intent(in) :: n
        doubleprecision, intent(in) :: nw
        doubleprecision, intent(in) :: w

        call this%init(n, nw, w)
    end function

    !!<summary>Initialize the class.</summary>
    subroutine initializeInterlayerAttenuation(this, n, nw, w)
        class(InterlayerAttenuation) :: this
        integer, intent(in) :: n
        doubleprecision, intent(in) :: nw
        doubleprecision, intent(in) :: w

        this%n = n
        this%nw = nw
        this%w = w
        this%xMin = -dble(nw + 3) * w
        this%xMax = -this%xMin
    end subroutine

    !!<summary>Destructor.</summary>
    subroutine deleteInterlayerAttenuation(this)
        type(InterlayerAttenuation) :: this

        if(allocated(this%kTable)) then
            deallocate(this%kTable)
        end if
    end subroutine

    !!<summary>Attenuation function.</summary>
    function attenuationFunction(this, x) result(y)
        class(InterlayerAttenuation) :: this
        doubleprecision, intent(in) :: x
        doubleprecision :: y

        if(x < -dble(this%nw) * this%w) then
            y = exp(-(x + dble(this%nw) * this%w) ** 2 / 2 / this%w ** 2)
        elseif(x <= dble(this%nw) * this%w) then
            y = 1.0
        else
            y = exp(-(x - dble(this%nw) * this%w) ** 2 / 2 / this%w ** 2)
        end if
    end function

    !!<summary>Generate kTable.</summary>
    subroutine genterate_kTable(this, kMin, kMax, nk)
        class(InterlayerAttenuation) :: this
        doubleprecision, intent(in) :: kMin, kMax
        integer, intent(in) :: nk

        ! local vars
        !!<local name="integrand">exp[i k x] attenuationFunction</local>
        complex(8), dimension(:), allocatable :: integrand
        !!<local name="simpsonCoefficient">Coefficients in Simpson's rule, 1 4 2 4 ... 2 4 1</local>
        doubleprecision, dimension(:), allocatable :: simpsonCoefficient
        ! iterator(s)
        integer :: i, j
        doubleprecision :: deltaK, deltaX, x

        this%kMin = kMin
        this%kMax = kMax
        this%nk = nk

        i = this%n + 1
        allocate(integrand(i))
        allocate(simpsonCoefficient(i))
        allocate(this%kTable(nk + 1, 2))

        ! Coefficients for Simpson's Rule
        simpsonCoefficient(1) = 1.0
        do i = 2, this%n
            if(mod(i, 2) == 0) then
                simpsonCoefficient(i) = 4.0
            else
                simpsonCoefficient(i) = 2.0
            end if
        end do
        simpsonCoefficient(this%n + 1) = 1.0

        deltaK = (kMax - kMin) / dble(nk)
        deltaX = (this%xMax - this%xMin) / dble(this%n)
        do i = 0, nk
            this%kTable(i + 1, 1) = kMin + dble(i) * deltaK
            ! calculate integrand list for a specific k
            do j = 0, this%n
                x = this%xMin + dble(j) * deltaX
                integrand(j + 1) = exp(cmplx(0.0, this%kTable(i + 1, 1) * x, kind=kind(1.0d0))) * this%attenuationFunction(x)
            end do
            ! Integrate via Simpson's Rule
            this%kTable(i + 1, 2) = real(sum(simpsonCoefficient * integrand)) * deltaX / 3.0
        end do

        deallocate(integrand)
        deallocate(simpsonCoefficient)
    end subroutine

    !!<summary>Fourier transformed attenuation function.</summary>
    function attenuationFunctionFT(this, k) result(y)
        class(InterlayerAttenuation) :: this
        doubleprecision, intent(in) :: k
        doubleprecision :: y

        ! local var(s)
        integer :: upperBoundIndex, lowerBoundIndex, tmpIndex

        if(k < this%kMin .or. k > this%kMax) then
            print *, __FILE__, __LINE__, "ERROR! k exceeds the range of Fourier transformed attenuation function. Modifying kMin and kMax might solve this problem."
        end if

        upperBoundIndex = this%nk + 1
        lowerBoundIndex = 1

        do while(upperBoundIndex /= lowerBoundIndex + 1)
            tmpIndex = (lowerBoundIndex + upperBoundIndex) / 2
            if(k > this%kTable(tmpIndex, 1)) then
                lowerBoundIndex = tmpIndex
            else
                upperBoundIndex = tmpIndex
            end if
        end do

        y = this%kTable(lowerBoundIndex, 2) + &
            (k                               - this%kTable(lowerBoundIndex, 1)) * &
            (this%kTable(upperBoundIndex, 2) - this%kTable(lowerBoundIndex, 2)) / &
            (this%kTable(upperBoundIndex, 1) - this%kTable(lowerBoundIndex, 1))
    end function
end module
