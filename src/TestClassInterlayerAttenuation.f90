program TestClassInterlayerAttenuation
    use ClassInterlayerAttenuation

    integer :: n
    doubleprecision :: w
    doubleprecision :: nw
    doubleprecision :: xMin
    doubleprecision :: xMax
    doubleprecision :: kMin
    doubleprecision :: kMax
    integer :: nk
    character(len = 64) :: arg, functionName
    type(InterlayerAttenuation) :: f
    doubleprecision :: x, k
    integer :: outUnit = 16

    call get_command_argument(0, functionName)
    do while(i <= command_argument_count())
        call get_command_argument(i, arg)

        select case(arg)
        case('-n')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) n
        case('-w')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) w
        case('--nw')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) nw
        case('--nk')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) nk
        case('--xMin')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) xMin
        case('--xMax')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) xMax
        case('--kMin')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) kMin
        case('--kMax')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) kMax
        end select

        i = i + 1
    end do

    call f%init(n, nw, w)
    call f%genterate_kTable(kMin, kMax, nk)

    open (unit=outUnit, file='test.txt', action="write", status="replace")
    do i = 0, n
        x = xMin + i * (xMax - xMin) / dble(n)
        write(outUnit, *) '{', x, ',', f%attenuationFunction(x), '},'
    end do

    write(outUnit, *) ' '
    write(outUnit, *) '==============================='
    write(outUnit, *) ' '

    do i = 0, nk
        k = kMin + i * (kMax - kMin) / dble(nk)
        write(outUnit, *) '{', k, ',', f%attenuationFunctionFT(k), '},'
    end do
end program
