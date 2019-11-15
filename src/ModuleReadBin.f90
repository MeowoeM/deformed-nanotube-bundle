!!1<summary>For readding the binary output cache file.</summary>
module ModuleReadBin
    integer :: outUnit = 16
contains
    subroutine openFile(fileName)
        !!<parameter name="filename">Input cache file.</parameter>
        character(len = *), intent(in) :: fileName

        open (unit=outUnit, form='unformatted', file=fileName, action="read")
    end subroutine

    function readInteger() result(i)
        integer :: i
        read (outUnit) i
    end function

    subroutine readIntegerVector(n, v)
        integer, intent(in) :: n
        integer, dimension(n) :: v

        read (outUnit) v
    end subroutine

    subroutine readDoubleVector(n, v)
        integer, intent(in) :: n
        doubleprecision, dimension(n) :: v

        read (outUnit) v
    end subroutine

    subroutine readComplex8Vector(n, v)
        integer, intent(in) :: n
        complex(8), dimension(n) :: v

        read (outUnit) v
    end subroutine

    subroutine readDoubleMatrix(m, n, A)
        integer, intent(in) :: m
        integer, intent(in) :: n
        doubleprecision, dimension(m, n) :: A

        read (outUnit) A
        !print *,A
    end subroutine

    subroutine readComplex8Matrix(m, n, A)
        integer, intent(in) :: m
        integer, intent(in) :: n
        complex(8), dimension(m, n) :: A

        read (outUnit) A
    end subroutine

    subroutine closeFile()
        close(outUnit)
    end subroutine
end module
