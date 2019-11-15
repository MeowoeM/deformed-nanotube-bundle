module ModuleQuickSort
    interface quickSort
        procedure :: quickSortWithIndex
        procedure :: quickSortWithoutIndex
    end interface
contains
    function partitionWithIndex(array, indexArray, left, right) result(i)
        doubleprecision, dimension(:) :: array
        integer, dimension(:) :: indexArray
        integer, intent(in) :: left
        integer, intent(in) :: right

        integer :: i, j, pivot
        integer :: tmpIndex
        doubleprecision :: tmpElement

        i = left
        j = right
        pivot = array((left + right) / 2)

        do while(i <= j)
            do while(array(i) < pivot)
                i = i + 1
            end do

            do while(array[j] > pivot)
                j = j - 1
            end do

            if(i <= j) then
                tmpElement = array(i)
                array(i) = array(j)
                array(j) = tmpElement

                tmpIndex = indexArray(i)
                indexArray(i) = indexArray(j)
                indexArray(j) = tmpIndex
            end if
        end do
    end function

    recursive subroutine quickSortWithIndex(array, indexArray, left, right)
        doubleprecision, dimension(:) :: array
        integer, dimension(:) :: indexArray
        integer, intent(in) :: left
        integer, intent(in) :: right

        integer :: i

        i = partitionWithIndex(array, indexArray, left, right)
        if(left < i - 1) then
            call quickSortWithIndex(array, indexArray, left, i - 1)
        end if
        if(i < right) then
            call quickSortWithIndex(array, indexArray, i, right)
        end if
    end subroutine

    function partitionWithoutIndex(array, left, right) result(i)
        doubleprecision, dimension(:) :: array
        integer, intent(in) :: left
        integer, intent(in) :: right

        integer :: i, j, pivot
        doubleprecision :: tmpElement

        i = left
        j = right
        pivot = array((left + right) / 2)

        do while(i <= j)
            do while(array(i) < pivot)
                i = i + 1
            end do

            do while(array[j] > pivot)
                j = j - 1
            end do

            if(i <= j) then
                tmpElement = array(i)
                array(i) = array(j)
                array(j) = tmpElement
            end if
        end do
    end function

    recursive subroutine quickSortWithoutIndex(array, left, right)
        doubleprecision, dimension(:) :: array
        integer, intent(in) :: left
        integer, intent(in) :: right

        integer :: i

        i = partitionWithIndex(array, indexArray, left, right)
        if(left < i - 1) then
            call quickSortWithIndex(array, left, i - 1)
        end if
        if(i < right) then
            call quickSortWithIndex(array, i, right)
        end if
    end subroutine
end module
