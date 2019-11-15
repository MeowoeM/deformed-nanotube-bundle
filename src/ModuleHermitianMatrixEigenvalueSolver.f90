module ModuleHermitianMatrixEigenvalueSolver

contains
    !!<summary>A wrapping of zfeast_hcsrev, a sparse hermitian eigenvalue solver
    function sparseHermitianSolver(uplo, n, values, rowIndex, columns, eMin, eMax, m0, e, x, m) result(info)
        !!<\parameter name="uplo">Must be 'U' or 'L' or 'F'. If uplo = 'U', values stores the upper triangular parts of A. If uplo = 'L', values stores the lower triangular parts of A. If uplo= 'F' , values stores the full matrix A.<\parameter>
        CHARACTER, intent(in) :: uplo
        !!<parameter name="n">Sets the size of the problem. n > 0.<\parameter>
        integer, intent(in) :: n
        ! the following three parameters describe the input sparse matrix in CSR format
        !!<parameter name="values">Array containing the nonzero elements of either the full matrix A or the upper or lower triangular part of the matrix A, as specified by uplo.<\parameter>
        complex(8), dimension(:), intent(in) :: values
        !!<parameter name="rowIndex">Array of length n + 1, containing indices of elements in the array values , such that values(i) is the index in the array values of the first non-zero element from the row i . The value of the last element values(n + 1) is equal to the number of non-zeros plus one.<\parameter>
        integer, dimension(n + 1), intent(in) :: rowIndex
        !!<parameter name="columns">Array containing the column indices for each non-zero element of the matrix A being represented in the array values . Its length is equal to the length of the array values.<\parameter>
        integer, dimension(:), intent(in) :: columns
        !!<parameter name="eMin">The lower bound of the interval to be searched for eigenvalues</parameter>
        doubleprecision, intent(in) :: eMin
        !!<parameter name="eMax">The upper bound of the interval to be searched for eigenvalues</parameter>
        doubleprecision, intent(in) :: eMax
        !!<parameter name="m0">On entry, specifies the initial guess for subspace dimension to be used, 0 < m0 <= n. Set m0 >= m where m is the total number of eigenvalues located in the interval [emin, emax]. If the initial guess is wrong, Extended Eigensolver routines return info=3.</parameter>
        integer, intent(inout) :: m0
        !!<parameter name="e">Array of length m0. On output, the first m entries of e are eigenvalues found in the interval.<\parameter>
        doubleprecision, dimension(m0), intent(inout) :: e
        !!<parameter name="x">On output, the first m columns of x contain the orthonormal eigenvectors corresponding to the computed eigenvalues e, with the i-th column of x holding the eigenvector associated with e(i).<\parameter>
        complex(8), dimension(n, m0), intent(inout) :: x
        !!<parameter name="m">The total number of eigenvalues found in the interval [emin, emax]: 0 < m <= m0.<\parameter>
        integer, intent(inout) :: m

        !!<parameter name="info">If info=0, the execution is successful.<\parameter>
        integer :: info

        !!<local name="fpm">FEAST parameter.</local>
        integer, dimension(128) :: fpm
        !!<local name="epsout">Outputting is supressed. Contains the relative error on the trace: |tracei - tracei-1| /max(|emin|, |emax|)</local>
        doubleprecision :: epsout
        !!<local name="loop">Outputting is supressed. Contains the number of refinement loop executed. Ignored on input.<\local>
        integer :: loop
        !!<local name="res">Outputting is supressed. Array of length m0. On exit, the first m components contain the relative residual vector<\local>
        doubleprecision, dimension(m0) :: res
        ! Check if the sparse matrix in CSR format is valid
        integer :: m0tmp
        if(.not.((size(values) == size(columns)) .and. (rowIndex(n + 1) == size(values) + 1))) then
            print *, __FILE__, __LINE__, "ERROR! The dimension of values array, columns array and rowIndex(n + 1) do not match!"
        end if

        call feastinit(fpm)
!print *, __FILE__, __LINE__, size(rowIndex), n
        m0tmp = m0
        call zfeast_hcsrev(uplo, n, values, rowIndex, columns, fpm, epsout, loop, eMin, eMax, m0tmp, e, x, m, res, info)
        call outputInfo(info)
    end function

    subroutine outputInfo(info)
        integer, intent(in) :: info

        select case(info)
            case(202)
                print *, __FILE__, __LINE__, "ERROR! Problem with size of the system n (n≤0)."
            case(201)
                print *, __FILE__, __LINE__, "ERROR! Problem with size of initial subspace m0 (m0≤0 or m0>n)."
            case(200)
                print *, __FILE__, __LINE__, "ERROR! Problem with emin,emax (emin≥emax)."
            case(4)
                print *, __FILE__, __LINE__, "WARNING! Successful return of only the computed subspace after call with fpm[14] = 1."
            case(3)
                print *, __FILE__, __LINE__, "WARNING! Size of the subspace m0 is too small (m0<m)."
            case(2)
                print *, __FILE__, __LINE__, "WARNING! No Convergence (number of iteration loops >fpm[4])."
            case(1)
                print *, __FILE__, __LINE__, "WARNING! No eigenvalue found in the search interval. See remark below for further details."
            case(-1)
                print *, __FILE__, __LINE__, "ERROR! Internal error for allocation memory."
            case(-2)
                print *, __FILE__, __LINE__, "ERROR! Internal error of the inner system solver. Possible reasons: not enough memory for inner linear system solver or inconsistent input."
            case(-3)
                print *, __FILE__, __LINE__, "ERROR! Internal error of the reduced eigenvalue solver. \n Possible cause: matrix B may not be positive definite. It can be checked by setting fpm(28) = 1 before calling an Extended Eigensolver routine, or by using LAPACK routines."
            case(-4)
                print *, __FILE__, __LINE__, "ERROR! Matrix B is not positive definite."
        end select
    end subroutine

    !!<summary>A wrapping of zheevr. The input hermitian matrix is a lower triangular one.</summary>
    subroutine hermitianSolver(a, w)
        !!<parameter name="a">The input matrix.</parameter>
        complex(8), dimension(:, :), intent(inout) :: a
        !!<parameter name="w">Array, size at least max(1, n), contains the selected eigenvalues in ascending order, stored in w(1) to w(m).</parameter>
        doubleprecision, allocatable, dimension(:), intent(inout) :: w

        !-----------------!
        ! Local Constants !
        !-----------------!

        !!<local name="n">The order of the matrix a.</local>
        integer :: n
        !!<local name="lda">The leading dimension of the array a.</local>
        integer :: lda
        !!<local name="LWMAX">The maximum size of local workspace</local>
        integer, parameter :: lwmax = 1000

        !---------------!
        ! Local Scalars !
        !---------------!

        !!<local name="info">If info = 0, the execution is successful. If info = -i, the i-th parameter had an illegal value.  If info = i, an internal error has occurred.</local>
        integer :: info
        !!<local name="lwork">The dimension of the array work.</local>
        integer :: lwork

        !--------------!
        ! Local Arrays !
        !--------------!

        !!<local name="rwork">Workspace array, size max(1, lwork).</local>
        doubleprecision, allocatable, dimension(:) :: rwork
        !!<local name="work">A workspace array, its dimension max(1, lwork).</local>
        complex(8), dimension(lwmax) :: work

        !----------------- Exacution Phase -----------------!
        n = size(a, 1)
        lda = n
        allocate(rwork(3 * n - 2))
        allocate(w(n))

        ! Query the optimal workspace.
        lwork = -1

        !print *, n, lda, vl, vu, il, iu, abstol, m, ldz, lrwork
        !print *, n
        CALL zheev('v', 'U', n, a, lda, w, work, lwork, rwork, info)
        lwork = max( min( lwmax, int( work( 1 ) ) ), 2 * n - 1 )

        ! Solve eigenproblem.
        !print *, n
        CALL zheev('v', 'U', n, a, lda, w, work, lwork, rwork, info)
                ! Check for convergence.
        IF( info.GT.0 ) THEN
            WRITE(*,*)'The algorithm failed to compute eigenvalues.'
            STOP
        END IF

        deallocate(rwork)
    end subroutine hermitianSolver

    subroutine hamiltonianPerturbationSolver(H, w, base1quantity, base2quantity)
        integer, intent(in) :: base1quantity
        integer, intent(in) :: base2quantity
        complex(8), dimension(:, :), intent(inout) :: H
        doubleprecision, allocatable, dimension(:), intent(inout) :: w

        ! local variables
        complex(8), dimension(:, :), allocatable :: delta_h
        doubleprecision, dimension(:), allocatable :: eigenvalue
        complex(8), dimension(:, :), allocatable :: eigenvector

        allocate(delta_h(2 *(base1quantity + base2quantity), 2 *(base1quantity + base2quantity)))
        allocate(eigenvalue(2 *(base1quantity + base2quantity)))
        allocate(eigenvector(2 *(base1quantity + base2quantity), 2 *(base1quantity + base2quantity)))

!        print *, __FILE__, __LINE__, 'a'
        allocate(w(size(H, 1)))

        do i = 1, size(H, 1)
            do j = 1, size(H, 2)
                delta_h(j, i) = conjg(H(i, j))
            end do
        end do
        H = H + delta_h

        do i = 1, base1quantity
            do j = 1, base1quantity
                delta_h(i, j) = 0
            end do
        end do
!        print *, __FILE__, __LINE__, 'b'
        do i = 2 * base1quantity + 1, 2 * base1quantity + base2quantity
            do j = 2 * base1quantity + 1, 2 * base1quantity + base2quantity
                delta_h(i, j) = 0
            end do
        end do
        do i = 1, base1quantity
            do j = 2 * base1quantity + 1, 2 * base1quantity + base2quantity
                delta_h(i, j) = H(i, j)
                delta_h(j, i) = conjg(H(i, j))
            end do
        end do

!        do i = 1, size(delta_h, 1)
!            do j = 1, size(delta_h, 2)
!                if(abs(delta_h(i, j)) /= 0) then
!                    print *, i, j, delta_h(i, j)
!                end if
!            end do
!        end do

!        print *, __FILE__, __LINE__, 'c'
        call noncouplingHamiltonianSolver(H, base1quantity, base2quantity, eigenvalue, eigenvector)
!        print *, __FILE__, __LINE__, 'd'
        call hamiltonianPerturbation(eigenvalue, eigenvector, delta_h, w, H, base1quantity, base2quantity)

        deallocate(delta_h)
        deallocate(eigenvalue)
        deallocate(eigenvector)
    end subroutine

    !!<summary>Given the solution to original eigenvalue problem: H_0, v_0 = e_0 v_0 and a small perturbation \delta H to H_0, solving H v = e v where H = H_0 + \delta H, e = e_0 + \delta e and v = v_0 + \delta v</summary>
    subroutine hamiltonianPerturbation(e0, v0, delta_h, e, v, base1quantity, base2quantity)
        integer, intent(in) :: base1quantity
        integer, intent(in) :: base2quantity
        doubleprecision, dimension(:), intent(in) :: e0
        complex(8), dimension(:, :), intent(in) :: v0
        complex(8), dimension(:, :), intent(in) :: delta_h
        doubleprecision, dimension(:), intent(out) :: e
        complex(8), dimension(:, :), intent(out) :: v

        ! local variables
        !!<local>Dimension of the Hamiltonian.</local>
        integer :: n
        !!<local>Just for fill in the parameter for zhemv.</local>
        complex(8), allocatable, dimension(:) :: y
        complex(8) :: tmp

        n = size(delta_h, 1)
        allocate(y(n))
        do i = 1, n
            if(i <= base1quantity) then
                do j = 1, 2 * base1quantity
                    y(j) = 0
                end do
                do j = 2 * base1quantity + 1, 2 * (base1quantity + base2quantity)
                    y(j) = delta_h(j, i) * v0(i, i) + delta_h(j, base1quantity + i) * v0(base1quantity + i, i)
                end do
            elseif(i <= 2 * base1quantity) then
                do j = 1, 2 * base1quantity
                    y(j) = 0
                end do
                do j = 2 * base1quantity + 1, 2 * (base1quantity + base2quantity)
                    y(j) = delta_h(j, i) * v0(i, i) + delta_h(j, i - base1quantity) * v0(i - base1quantity, i)
                end do
            elseif(i <= 2 * base1quantity + base2quantity) then
                do j = 1, 2 * base1quantity
                    y(j) = delta_h(j, i) * v0(i, i) + &
                           delta_h(j, i + base2quantity) * v0(i + base2quantity, i)
                end do
                do j = 2 * base1quantity + 1, 2 * (base1quantity + base2quantity)
                    y(j) = 0
                end do
            elseif(i <= 2 * (base1quantity + base2quantity)) then
                do j = 1, 2 * base1quantity
                    y(j) = delta_h(j, i) * v0(i, i) + &
                           delta_h(j, i - base2quantity) * v0(i - base2quantity, i)
                end do
                do j = 2 * base1quantity + 1, 2 * (base1quantity + base2quantity)
                    y(j) = 0
                end do
            end if
            e(i) = e0(i) + dot_product(v0(:, i), y)
            v(:, i) = v0(:, i)
            do j = 1, n
                if(i /= j) then
                    tmp = dot_product(v0(:, j), y)
                    if(abs(tmp) > 1e-8) then
!                        print *, __FILE__, __LINE__, i, j, tmp
                        e(i) = e(i) + tmp ** 2 / (e0(i) - e0(j))
                        v(:, i) = v(:, i) + tmp / (e0(i) - e0(j)) * v0(:, j)
!                        if(abs(e0(i) - e0(j)) < 1e-8) then
!                            print *, __FILE__, __LINE__, i, j, e0(i), e0(j)
!                        end if
                    end if
                end if
            end do
!            do j = 1, n
!                if(v(j, i) /= 0) then
!                    print *, i, j, v(j, i)
!                end if
!            end do
        end do

!        print *, __FILE__, __LINE__, e - e0
        deallocate(y)
    end subroutine

    !!<summary>Sovle the Hamiltonian, H, of two non-interacting nanotube. H = H_1 \oplus H_2</summary>
    subroutine noncouplingHamiltonianSolver(H, base1quantity, base2quantity, eigenvalue, eigenvector)
        integer, intent(in) :: base1quantity
        integer, intent(in) :: base2quantity
        complex(8), dimension(2 *(base1quantity + base2quantity), 2 *(base1quantity + base2quantity)), intent(in) :: H
        doubleprecision, dimension(2 *(base1quantity + base2quantity)), intent(out) :: eigenvalue
        complex(8), dimension(2 *(base1quantity + base2quantity), 2 *(base1quantity + base2quantity)), intent(out) :: eigenvector

        ! local variables
        complex(8), dimension(2) :: tmpEigenvalue
        complex(8), dimension(2, 2) :: tmpEigenvector

        ! TODO: If ram runs out, varible 'eigenvector' can be combined to 'H' which turns into intent(inout).

        do i = 1, base1quantity
            call eigensystem2by2(H(i, i + base1quantity), H(i + base1quantity, i), tmpEigenvalue, tmpEigenvector)
            eigenvalue(i                ) = dble(tmpEigenvalue(1))
            eigenvalue(i + base1quantity) = dble(tmpEigenvalue(2))
            eigenvector(i                , i                ) = tmpEigenvector(1, 1)
            eigenvector(i                , i + base1quantity) = tmpEigenvector(1, 2)
            eigenvector(i + base1quantity, i                ) = tmpEigenvector(2, 1)
            eigenvector(i + base1quantity, i + base1quantity) = tmpEigenvector(2, 2)
        end do

        do i = 1, base2quantity
            call eigensystem2by2(H(2 * base1quantity + i, 2 * base1quantity + i + base2quantity), &
                                 H(2 * base1quantity + i + base2quantity, 2 * base1quantity + i), tmpEigenvalue, tmpEigenvector)
            eigenvalue(2 * base1quantity + i                ) = dble(tmpEigenvalue(1))
            eigenvalue(2 * base1quantity + i + base2quantity) = dble(tmpEigenvalue(2))
            eigenvector(2 * base1quantity + i                , 2 * base1quantity + i                ) = tmpEigenvector(1, 1)
            eigenvector(2 * base1quantity + i                , 2 * base1quantity + i + base2quantity) = tmpEigenvector(1, 2)
            eigenvector(2 * base1quantity + i + base2quantity, 2 * base1quantity + i                ) = tmpEigenvector(2, 1)
            eigenvector(2 * base1quantity + i + base2quantity, 2 * base1quantity + i + base2quantity) = tmpEigenvector(2, 2)
        end do
    end subroutine

    !!<summary>((0, a), (b, 0))</summary>
    subroutine eigensystem2by2(a, b, eigenvalue, eigenvector)
        complex(8), intent(in) :: a
        complex(8), intent(in) :: b
        complex(8), dimension(2), intent(out) :: eigenvalue
        complex(8), dimension(2, 2), intent(out) :: eigenvector

        eigenvalue = (/ sqrt(a * b), -sqrt(a * b) /)
        eigenvector(1, :) = (/ -sqrt(a / (a + b)), sqrt(a / (a + b)) /)
        eigenvector(2, :) = (/ sqrt(b / (a + b)), sqrt(b / (a + b)) /)
    end subroutine
end module
