	program circlearea
	real r, area, pi
	parameter (pi = 3.14159)
c       This program computes the area of a circle.
	print *, "What is the radius?"
	read *, r
	area = pi * r ** 2
	print *, "The area is", area
	print *, "Bye!"
	end