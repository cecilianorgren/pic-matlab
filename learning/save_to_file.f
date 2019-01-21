	program save_to_file
	real a, b
c       This program computes the area of a circle.
	do i = 1,10
		a(i) = i
		b(i) = -i
	end do
	write(1) a,b
	end