	subroutine inv(A,QT,R,NA)
	implicit none
	integer*4 NA
	real*8 A(NA*NA), QT(NA*NA), R(NA*NA)
	integer*4 n,m,k
	real*8 tmp
	call qr(A,QT,R,NA)
	
	do n= NA-1,0,-1
		tmp = R(n*NA + n + 1)
		do m = 0, NA-1
			R (n*NA + m + 1) = R (n*NA + m + 1)/ tmp
			QT(n*NA + m + 1) = QT(n*NA + m + 1)/ tmp
		enddo
		do m=n+1,NA-1
			tmp = R(n*NA + m + 1)
			R( n*NA + m + 1) = 0
			do k=0,NA-1
				QT(n*NA + k + 1) = QT(n*NA + k + 1) - QT( m*NA + k + 1) * tmp
			enddo
		enddo
	enddo
	
	end
	
