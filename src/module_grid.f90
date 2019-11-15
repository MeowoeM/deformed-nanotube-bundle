module grid

contains
	subroutine grid_add(resX, resY, data, X, Y, Z, n)
		integer, intent(in) :: resX
		integer, intent(in) :: resY
		doubleprecision, dimension(n, 3) :: data
		doubleprecision, dimension(resY, resX) :: X
		doubleprecision, dimension(resY, resX) :: Y
		doubleprecision, dimension(resY, resX) :: Z
		integer :: n

		doubleprecision :: xMax
		doubleprecision :: xMin
		doubleprecision :: yMax
		doubleprecision :: yMin
		doubleprecision :: deltaX
		doubleprecision :: deltaY

		integer :: xi, yi

		xMax = arrayMax(data(:, 1), n)
		xMin = arrayMin(data(:, 1), n)
		yMax = arrayMax(data(:, 2), n)
		yMin = arrayMin(data(:, 2), n)
		deltaX = (xMax - xMin) / dble(resX - 1)
		deltaY = (yMax - yMin) / dble(resY - 1)

		
		do i = 0, resX - 1
			do j = 1, resY
				X(j, i + 1) = xMin + dble(i) * deltaX
			end do
		end do

		do j = 1, resX
			do i = 0, resY - 1
				Y(i + 1, j) = yMin + dble(i) * deltaY
			end do
		end do

		do i = 1, resX
			do j = 1, resY
				Z(j, i) = 0
			end do
		end	do

		do i = 1, n
			xi = floor((data(i, 1) - xMin) / deltaX + 0.5) + 1
			yi = floor((data(i, 2) - yMin) / deltaY + 0.5) + 1
			if(xi > 0 .and. xi <= resX .and. yi > 0 .and. yi <= resY) then
				Z(yi, xi) = Z(yi, xi) + data(i, 3)
			end if
		end do

	end subroutine

	subroutine grid_max(resX, resY, data, X, Y, Z, n)
		integer, intent(in) :: resX
		integer, intent(in) :: resY
		doubleprecision, dimension(n, 3) :: data
		doubleprecision, dimension(resY, resX) :: X
		doubleprecision, dimension(resY, resX) :: Y
		doubleprecision, dimension(resY, resX) :: Z
		integer :: n

		doubleprecision :: xMax
		doubleprecision :: xMin
		doubleprecision :: yMax
		doubleprecision :: yMin
		doubleprecision :: deltaX
		doubleprecision :: deltaY

		integer :: xi, yi

		xMax = arrayMax(data(:, 1), n)
		xMin = arrayMin(data(:, 1), n)
		yMax = arrayMax(data(:, 2), n)
		yMin = arrayMin(data(:, 2), n)
		deltaX = (xMax - xMin) / dble(resX - 1)
		deltaY = (yMax - yMin) / dble(resY - 1)
		
		do i = 0, resX - 1
			do j = 1, resY
				X(j, i + 1) = xMin + dble(i) * deltaX
			end do
		end do

		do j = 1, resX
			do i = 0, resY - 1
				Y(i + 1, j) = yMin + dble(i) * deltaY
			end do
		end do

		do i = 1, resX
			do j = 1, resY
				Z(j, i) = 0
			end do
		end	do

		do i = 1, size(data, 1)
			xi = floor((data(i, 1) - xMin) / deltaX + 0.5)
			yi = floor((data(i, 2) - yMin) / deltaY + 0.5)
			if(xi > 0 .and. xi <= resX .and. yi > 0 .and. yi <= resY) then
				if(data(i, 3) > Z(yi, xi)) then
					Z(yi, xi) = data(i, 3)
				end if
			end if
		end do
	end subroutine

	subroutine dos(resY, data, x, y, n)
		integer, intent(in) :: resY
		doubleprecision, dimension(n, 3) :: data
		doubleprecision, dimension(resY) :: x
		doubleprecision, dimension(resY) :: y
		integer :: n

		doubleprecision :: yMax, yMin, deltaY, yi

		yMax = arrayMax(data(:, 2), n)
		yMin = arrayMin(data(:, 2), n)
		deltaY = (yMax - yMin) / dble(resY - 1)

		do i = 0, resY - 1
			y(i + 1) = yMin + dble(i) * deltaY
		end do

		do i = 1, resY
			x(i) = 0
		end do

		do i = 1, n
			yi = floor((data(i, 2) - yMin) / deltaY + 0.5)
			if(yi > 0 .and. yi <= resY) then
				x(yi) = x(yi) + data(i, 3)
			end if
		end do

		do i = 1, resY
			x(i) = x(i) / deltaY
		end do
	end subroutine

	function arrayMax(array, n) result(res)
		doubleprecision, dimension(n) :: array
		doubleprecision :: res
		integer :: n

		res = array(1)
		do i = 1, n
			if(array(i) > res .and. array(i) < 100) then
				res = array(i)
			end if
		end do
	end function

	function arrayMin(array, n) result(res)
		doubleprecision, dimension(n) :: array
		doubleprecision :: res
		integer :: n

		res = array(1)
		do i = 1, n
			if(array(i) < res .and. array(i) > -100) then
				res = array(i)
			end if
		end do
	end function

	subroutine window(xMin, xMax, yMin, yMax, data, array, m, n, count)
		doubleprecision :: xMin
		doubleprecision :: xMax
		doubleprecision :: yMin
		doubleprecision :: yMax
		doubleprecision, dimension(n, 3) :: data
		doubleprecision, dimension(m, 3) :: array
		integer :: m, n, count

		count = 0
		do i = 1, n
			if((data(i, 1) > xMin) .and. &
			   (data(i, 1) < xMax) .and. &
			   (data(i, 2) > yMin) .and. &
			   (data(i, 2) < yMax)) then
				count = count + 1
				if(count > m) then
					count = m
					exit
				end if
				array(count, 1) = data(i, 1)
				array(count, 2) = data(i, 2)
				array(count, 3) = data(i, 3)
			end if
		end do
	end subroutine
end module