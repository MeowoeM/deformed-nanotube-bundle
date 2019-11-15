program main
    use Lattice2d
    use ClassSingleWallNanotube
    use ClassNanotubeBundle
    use class_kBase
    use ClassHamiltonianConstructor
    use ModuleHermitianMatrixEigenvalueSolver
    use ClassSpectralFunction
!    use ModuleReadInputFile
    use ClassVectorInteger
    use ClassVectorIntegerArray1d
    use ClassVectorDouble
    use ClassVectorDoubleArray1d
    use ClassVectorComplex8
    use Fox_DOM, only: &
        Node, &
        NodeList, &
        parseFile, &
        getDocumentElement, &
        extractDataAttribute, &
        getElementsByTagName, &
        extractDataAttribute, &
        item, &
        getNodeListLength => getLength

    implicit none
    type(NanotubeBundle) :: bundle
    type(SingleWallNanotube) :: tmpTube, tmpTube2
    type(kBase) :: base1, base2
    type(HamiltonianConstructor) :: constructor
    type(VectorComplex8) :: valuesVector
    complex(8), dimension(:, :), allocatable :: H
    doubleprecision, allocatable, dimension(:) :: w
    complex(8), dimension(:), allocatable :: values
    type(VectorInteger) :: columnsVector
    integer, dimension(:), allocatable :: columns
    type(VectorInteger) :: rowIndexVector
    integer, dimension(:), allocatable :: rowIndex
    type(SpectralFunction) :: spectrum

    doubleprecision :: energy0, energy1, energy2, deltaE
    integer :: totalDimension, flag1, eigenvalueFoundNumber
    integer, dimension(2) :: fittingCoefficient

    !!<local>Chiral Vector.</local>
    integer :: m1, n1, m2, n2
    !!<local name="xi">\xi indicates whether K point(+1) or K' point(-1).</local>
    integer :: xi
    !!<local name="shift1">Shift from K or K'.</local>
    doubleprecision, dimension(2) :: shift1
    !!<local name="shift2">Shift from K or K'.</local>
    doubleprecision, dimension(2) :: shift2
    !!<local name="kMAx"> The radius of cutoff circle.</local>
    doubleprecision :: kMax
    !!<local name="shiftCenter">The center of a series of cutoff circle because \bar{K} or \bar{K}' do not necessarily falls on the cutting line.</local>
    doubleprecision, dimension(2) :: shiftCenter
    !!<local name="shiftStep">The step that two neighbouring cutoff circle centers differ.</local>
    !doubleprecision :: shiftStep
    ! Controls the density of k points in k_y direction
    doubleprecision :: kStep
    !!<local name="interlayerInteraction">Interlayer interaction intensity.</local>
    doubleprecision :: interlayerInteraction
    !!<local name="sampleNumberRS">Number of samples for interlayer attenuation function in real space.</local>
    integer :: sampleNumberRS
    !!<local name="sampleNumberFT">Number of samples for Fourier transformed interlayer attenuation function.</local>
    integer :: sampleNumberFT
    !!<local name="flatRegionWidthHalf">A half of the width of flat intertube interacting region in unit of Gaussian half width(see InterlayerAttenuation class for details).</local>
    doubleprecision :: flatRegionWidthHalf
    !!<local name="isInfinite">Whether it is an infinite bundle or one consisting of only 2 tubes.</local>
    logical :: isInfinite
    !!<local name="diameterRatio">Ratio of the diameter of tube 2 against that of tube 1.</local>
    integer :: diameterRatio


    character(len = 64) :: arg, functionName, inputFileName, outputfileName

    integer :: i, j, kk, ll, mode, tmpInteger
    integer, dimension(2) :: cuttingLineCoordinate1, cuttingLineCoordinate2
    integer :: outUnit = 16

    !!<local name="kBaseIndex">Storing the index of k points on the cutting lines of interest.</local>
    !type(VectorInteger) :: kBaseIndex
    !!<local name="output">Output for mode 1, the first column stores the ky coordinates of k points, the second column energy and the third A</local>
    type(VectorDoubleArray1d), dimension(:), allocatable :: outputBase1, outputBase2
    !!<local name="cuttingLineRange">Select cutting lines for sepctral function with the bound, only available in mode 1.</local>
    doubleprecision :: cuttingLineRange
    !!<local name="cuttingLineCategory">Storing which cutting line k points belong to.</local>
    type(VectorDouble) :: cuttingLineCategory1, cuttingLineCategory2
    type(VectorInteger) :: tmpIntVector1, tmpIntVector2

    doubleprecision, dimension(2) :: tmpVector
    doubleprecision :: tmp, tmp2
    doubleprecision :: accuracy
    doubleprecision :: halfWidth
    integer :: shiftN

    doubleprecision, allocatable, dimension(:, :) :: baseCoordinates

    type(VectorDoubleArray1d) :: x, y

    type(VectorIntegerArray1d) :: couplingComponentVector
    integer, dimension(2) :: tmpIntVector
    type(Node), pointer :: inputFile, root, chirality, couplingComponent
    type(NodeList), pointer :: chiralityList, couplingComponentList

    outputfileName = "result.txt"
    interlayerInteraction = 1.0

    !------------- Reading Input Arguments -------------!
    i = 1
    cuttingLineRange = 9
    halfWidth = 0
    mode = 1
    accuracy = 1e-4
    call get_command_argument(0, functionName)
    do while(i <= command_argument_count())
        call get_command_argument(i, arg)

        select case(arg)
        case('-h', '--help')
            call help(functionName)
            stop
        case('--m1')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) m1
        case('-i', '--input')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) inputFileName
        case('--n1')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) n1
        case('-r')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) kMax
        case('--xi')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) xi
        case('-n')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) shiftN
        case('-o')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) outputfileName
        case('-c')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) cuttingLineRange
        case('-w')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) halfWidth
        case('-m')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) mode

            ! check whether mode is legitimate
            if((mode /= 0) .and. (mode /= 1)) then
                print *, __FILE__, __LINE__, " ERROR! Mode(-m) must be 0 or 1."
            end if
        case('--interaction')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, *) interlayerInteraction
        case default
            print '(a,a,/)', 'Unrecognized command-line option: ', arg
            call help(functionName)
            stop
        end select

        i = i + 1
    end do

    inputFile => parseFile(inputFileName)
    root => getDocumentElement(inputFile)
    chiralityList => getElementsByTagName(root, "chirality")
    couplingComponentList => getElementsByTagName(root, "couplingComponent")
    print *, __FILE__, __LINE__, inputFileName, 'read.'

    call extractDataAttribute(root, "cutoffCircleRadius", kMax)
    print *, __FILE__, __LINE__, "cutoffCircleRadius", kMax
    call extractDataAttribute(root, "xi", xi)
    print *, __FILE__, __LINE__, "xi", xi
    call extractDataAttribute(root, 'shiftN', shiftN)
    print *, __FILE__, __LINE__, "shiftN", shiftN
    call extractDataAttribute(root, 'outputFileName', outputfileName)
    print *, __FILE__, __LINE__, "outputFileName", outputFileName
    call extractDataAttribute(root, 'cuttingLineRange', cuttingLineRange)
    print *, __FILE__, __LINE__, "cuttingLineRange", cuttingLineRange
    call extractDataAttribute(root, 'couplingAmplifier', interlayerInteraction)
    print *, __FILE__, __LINE__, "couplingAmplifier", interlayerInteraction
    call extractDataAttribute(root, 'sampleNumberRS', sampleNumberRS)
    print *, __FILE__, __LINE__, "sampleNumberRS", sampleNumberRS
    call extractDataAttribute(root, 'sampleNumberFT', sampleNumberFT)
    print *, __FILE__, __LINE__, "sampleNumberFT", sampleNumberFT
    call extractDataAttribute(root, 'flatRegionWidthHalf', flatRegionWidthHalf)
    print *, __FILE__, __LINE__, "flatRegionWidthHalf", flatRegionWidthHalf, "nw"
    call extractDataAttribute(root, 'isInfinite', isInfinite)
    print *, __FILE__, __LINE__, "isInfinite", isInfinite
    call extractDataAttribute(root, 'diameterRatio', diameterRatio)
    print *, __FILE__, __LINE__, "diameterRatio", diameterRatio
    call extractDataAttribute(root, 'halfWidth', halfWidth)
    if(halfWidth < 1e-6) then
        print *, __FILE__, __LINE__, "halfWidth auto"
    else
        print *, __FILE__, __LINE__, "halfWidth", halfWidth
    end if
    call extractDataAttribute(root, 'mode', mode)
    print *, __FILE__, __LINE__, 'mode', mode


    chirality => item(chiralityList, 0)
    call extractDataAttribute(chirality, 'm', m1)
    call extractDataAttribute(chirality, 'n', n1)

    couplingComponentVector = VectorIntegerArray1d(2)
    do i = 0, getNodeListLength(couplingComponentList) - 1
        couplingComponent => item(couplingComponentList, i)
        call extractDataAttribute(couplingComponent, 'Gm1', tmpIntVector(1))
        call extractDataAttribute(couplingComponent, 'Gm2', tmpIntVector(2))
        call couplingComponentVector%append(tmpIntVector)
    end do

    !print *,"1"
    !print *, n1, n2, shiftCenter, shiftN, kMax
    ! zigzag
    if((dble(n1) / dble(m1) < 0.2) .and. (dble(n1) / dble(m1) > -0.2)) then
        bundle = newNanotubeBundle(newSingleWallNanotube(m1, n1, dble(0.0)), &
                                   newSingleWallNanotube((m1 + n1) * diameterRatio, -n1 * diameterRatio, &
                                   dble(0.0)), dble(0.0), isInfinite)
        print *, __FILE__, __LINE__, 'Chiral vector is close to zigzag direction.'
    ! armchair
    else if((dble(n1) / dble(m1) > 0.8) .and. (dble(n1) / dble(m1) < 1.8)) then
        bundle = newNanotubeBundle(newSingleWallNanotube(m1, n1, dble(0.0)), &
                                   newSingleWallNanotube(n1 * diameterRatio, m1 * diameterRatio, &
                                   dble(0.0)), dble(0.0), isInfinite)
        print *, __FILE__, __LINE__, 'Chiral vector is close to armchair direction.'
    else
        print *, __FILE__, __LINE__, 'Chiral vector should be close to zigzag/armchair direction.'
    end if
    tmpVector = bundle%get_gm1()
    ! print *, __FILE__, __LINE__, bundle%get_gm1(), bundle%get_gm2()
    ! print *, __FILE__, __LINE__, bundle%get_gm1(), bundle%get_gm2()
    tmpTube = bundle%getNanotube(1)
    tmpTube2 = bundle%getNanotube(2)

    xi = 1
    open (unit=outUnit, form='unformatted', file=outputfileName, action="write", status="replace")

    print *, 'mode:',  mode
    if(mode == 1) then

        !print *, "continue?"
        !read (*, *) tmpInteger
        print *, __FILE__, __LINE__, 'gm1', bundle%get_gm1()
        print *, __FILE__, __LINE__, 'gm2', bundle%get_gm2()
        print *, __FILE__, __LINE__, 'kbar', bundle%get_kBar()
        print *, __FILE__, __LINE__, 'deltaky', dot_product(bundle%get_gm1(), (/0.0, 1.0/)) / dble(shiftN)
        print *, __FILE__, __LINE__, "Initializing..."
        shift1 = bundle%calcShift(bundle%get_kBar(), 1)
        tmpTube = bundle%getNanotube(1)
        base1 = new_kBase(bundle%get_kBar(), shift1, 8.0 * atan(1.0) / tmpTube%getChiralVectorLength(), abs(dot_product((/0.0, 1.0/), bundle%get_gm1())), kMax)
        shift2 = bundle%calcShift(bundle%get_kBar(), 2)
        tmpTube = bundle%getNanotube(2)
        base2 = new_kBase(bundle%get_kBar(), shift2, 8.0 * atan(1.0) / tmpTube%getChiralVectorLength(), abs(dot_product((/0.0, 1.0/), bundle%get_gm1())), kMax)
        totalDimension = 2.0 * (base1%getQuantity() + base2%getQuantity())
        print *, __FILE__, __LINE__, totalDimension
        call spectrum%initSpectralFunction(base1, base2, cuttingLineRange)
        call spectrum%printCuttinglinePosition()
        open(unit=6, form='FORMATTED', carriagecontrol='FORTRAN')

        if(abs(halfWidth - 0.0) < accuracy) then
            call bundle%setHalfWidthAuto((tmpTube%getChiralVectorLength() + tmpTube2%getChiralVectorLength()) / 2.0)
            print *, __FILE__, __LINE__, "halfWidth", bundle%getHalfWidth()
        else
            call bundle%setHalfWidth(halfWidth)
        end if

        call constructor%set_xCutoff(min(4.0 / bundle%getHalfWidth(), base1%getCutoffRadius()))
        call bundle%initAttenuation(sampleNumberRS, flatRegionWidthHalf)
        call bundle%genterate_kTable(-constructor%get_xCutoff(), constructor%get_xCutoff(), sampleNumberFT)
        do i = 0, shiftN - 1
            write(6, 100) 'Progress: ', i + 1, '/', shiftN
            100 format('+', A, I0, A, I0)
            shift1 = shift1 + (/dble(0), dot_product((/0.0, 1.0/), bundle%get_gm1()) / dble(shiftN)/)
            call base1%setShift(shift1)
            shift2 = shift2 + (/dble(0), dot_product((/0.0, 1.0/), bundle%get_gm1()) / dble(shiftN)/)
            call base2%setShift(shift2)
!            print *, __FILE__, __LINE__, "Constructing Hamiltonian matrix..."
            call constructor%constructHamiltonianMatrix(bundle, base1, base2, xi, H, interlayerInteraction, couplingComponentVector)

!            tmpIntVector1 = base1%get_kPointNumberList()
!            tmpIntVector2 = base2%get_kPointNumberList()
!
!            do j = 1, tmpIntVector1%getLength()
!                print *, base1%getCoordinate(tmpIntVector1%getSum(1, j - 1) + 1)
!            end do
!
!            do j = 1, tmpIntVector2%getLength()
!                print *, base2%getCoordinate(tmpIntVector2%getSum(1, j - 1) + 1)
!            end do
!
!            do kk = 1, 2 * base1%getQuantity()
!                do ll = 2 * base1%getQuantity() + 1, totalDimension
!                    if(&
!                        .not.(((kk > tmpIntVector1%getSum(1, 1)) .and. (kk <= tmpIntVector1%getSum(1, 2))) .or. &
!                        ((kk > base1%getQuantity() + tmpIntVector1%getSum(1, 1)) .and. (kk <= base1%getQuantity() + tmpIntVector1%getSum(1, 2)))) .or. &
!                        .not.(((ll > 2 * base1%getQuantity() + tmpIntVector2%getSum(1, 1)) .and. (ll <= 2 * base1%getQuantity() + tmpIntVector2%getSum(1, 2))) .or. &
!                        ((ll > 2 * base1%getQuantity() + base2%getQuantity() + tmpIntVector2%getSum(1, 1)) .and. (ll <= 2 * base1%getQuantity() + base2%getQuantity() + tmpIntVector2%getSum(1, 2)))) &
!                    ) then
!                        H(kk, ll) = 0
!                    end if
!                end do
!            end do

!            do kk = 1, 2 * base1%getQuantity()
!                do ll = 2 * base1%getQuantity() + 1, totalDimension
!                    if(&
!                        .not.(((kk > 0) .and. (kk <= tmpIntVector1%get(1))) .or. &
!                        ((kk > base1%getQuantity()) .and. (kk <= base1%getQuantity() + tmpIntVector1%get(1)))) .or. &
!                        .not.(((ll > 2 * base1%getQuantity()) .and. (ll <= 2 * base1%getQuantity() + tmpIntVector2%get(1))) .or. &
!                        ((ll > 2 * base1%getQuantity() + base2%getQuantity()) .and. (ll <= 2 * base1%getQuantity() + base2%getQuantity() + tmpIntVector2%get(1)))) &
!                    ) then
!                        H(kk, ll) = 0
!                    end if
!                end do
!            end do

!            x = newVectorDoubleArray1d(2)
!            do j = 1, base1%getQuantity()
!                call x%append(base1%getCoordinate(j))
!            end do
!            write (outUnit) base1%getQuantity()
!            call x%writeVector(outUnit)
!
!            y = newVectorDoubleArray1d(2)
!            do j = 1, base2%getQuantity()
!                call y%append(base2%getCoordinate(j))
!            end do
!            write (outUnit) base2%getQuantity()
!            call y%writeVector(outUnit)
!
!            write (outUnit) totalDimension
!            write (outUnit) H
!
!
!        close(unit=outUnit)

    !        write (outUnit) totalDimension
    !        write (outUnit) valuesVector%getLength()
    !        call valuesVector%writeVector(outUnit)
    !        write (outUnit) columnsVector%getLength()
    !        call columnsVector%writeVector(outUnit)
    !        write (outUnit) rowIndexVector%getLength()
    !        call rowIndexVector%writeVector(outUnit)
    !        close(unit=outUnit)
            !write(outUnit) size(H, 1)
            !write(outUnit) H
!            print *, __FILE__, __LINE__, "Size of H:", size(H, 1)
!            print *, __FILE__, __LINE__, "Solving it..."
    !        tmpInteger = base%getQuantity()
            !print *, H(1:(2*tmpInteger),1)
            !print *, '------------'
            !print *, H((2*tmpInteger+1):(4*tmpInteger),1)

            call hermitianSolver(H, w)
!            kk = base1%getQuantity()
!            ll = base2%getQuantity()
!            call hamiltonianPerturbationSolver(H, w, kk, ll)
!            print *, __FILE__, __LINE__, "writing results..."
            call spectrum%extend(base1,base2, w, H, size(w))


!-------------- output base coordinates, eigenvalues and eigenvectors --------------
!            kk = base1%getQuantity()
!            allocate(baseCoordinates(kk, 2))
!            do j = 1, kk
!                tmpVector = base1%getCoordinate(j)
!                baseCoordinates(j, 1) = tmpVector(1)
!                baseCoordinates(j, 2) = tmpVector(2)
!            end do
!            write (outUnit) kk
!            write (outUnit) baseCoordinates
!            deallocate(baseCoordinates)
!            ll = base2%getQuantity()
!            allocate(baseCoordinates(ll, 2))
!            do j = 1, ll
!                tmpVector = base2%getCoordinate(j)
!                baseCoordinates(j, 1) = tmpVector(1)
!                baseCoordinates(j, 2) = tmpVector(2)
!            end do
!            write (outUnit) ll
!            write (outUnit) baseCoordinates
!            deallocate(baseCoordinates)
!            write (outUnit) totalDimension
!            write (outUnit) w
!            write (outUnit) H
!-------------- output base coordinates, eigenvalues and eigenvectors --------------

            deallocate(H)
            deallocate(w)

    !        i = output%getLength()
    !        write (outUnit) i
    !        call output%writeVector(outUnit)
    !        call cuttingLineCategory%writeVector(outUnit)
        end do


        print *, __FILE__, __LINE__, "Done!"
        call spectrum%output(outUnit)
        close(unit=outUnit)
    end if

    close(unit=outUnit)
contains
    subroutine help(functionName)
        character(len = 32), intent(in) :: functionName

        print *,"NAME"
        print *,"   ",functionName,"    - calculate the band structure of double nanotube bundle"
        print *,""
        print *,"SYNOPSIS"
        print *,"   ",functionName," --help --n1 --n2 -x -y -n -r -o"
        print *,""
        print *,"OPTIONS"
        print *,"   -s              Controls the density of k points in k_y direction."
        print *,"   -r              The radius of cutoff circles"
        print *,"   --xi            K point(+1) or K' point(-1)"
        print *,"   --interaction   interlayer interaction (on/off)"
        print *,"   -i              input file name(len < 32 chars)"
        print *,"   -o              output file name(len < 32 chars)"
        print *,"   -m              output mode, 0 or 1"
        print *,"                       0: write kBase, Hamiltionians and their eigen energies and vectors"
        print *,"                       1: write the A on the cutting lines in the middle"
        print *,"   -c              number of cutting lines for spectral function, only available in mode 1"
        print *,"   -w              specify halfwidth of the interlayer interaction Gaussian function factor manually"
        print *,"   --m0            specify the subspace dimension maximum for solving sparse matrix eigenvalues, default 400."
        print *,"   --eRange        energy range for searching eigenenergy, (-eRange, eRange)."
        print *,""
        print *,"OUTPUT FORMAT:"
        print *,"   number of cutoff circles, say n1"
        print *,"   repeating n1 times"
        print *,"       k bases"
        print *,"       number of row/column"
        print *,"       eigenvalues"
        print *,"       eigenvectors"
        print *,"* All of them are seperated by \n"

    end subroutine help
end program main
