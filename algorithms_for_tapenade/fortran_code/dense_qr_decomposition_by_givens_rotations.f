	SUBROUTINE QR(A,QT,R,NA)
	implicit none
	integer*4 NA
	real*8 A(NA*NA), QT(NA*NA), R(NA*NA)
	
	
	real*8 tmp, at,bt,rt,c,s,Rnk, QTnk
	integer*4 n,m,k 
c     PREPARE Q AND R
	do n=0,NA-1
		do m=0,NA-1
			if ( n .EQ. m) then
				tmp = 1.0
			else 
				tmp = 0.0
			endif
			QT(n*NA+m + 1)=tmp
		enddo
	enddo
	 
	do n=1,NA*NA
		R(n)=A(n)
	enddo
	
c	MAIN ALGORITHM
	do n=0,NA-1
		do m=n+1,NA-1
			at = R(n*NA+n+1)
			bt = R(m*NA+n+1)
			rt = sqrt( at*at + bt*bt)
			c = at/rt
			s = bt/rt
			do k=1,NA
c				UPDATE R
				Rnk = R(n*NA + k )
				R(n*NA + k ) = c*Rnk + s*R(m*NA + k )
				R(m*NA + k ) =-s*Rnk + c*R(m*NA + k )
c				UPDATE Q
				QTnk = QT(n*NA + k )
				QT(n*NA + k ) = c * QTnk + s*QT(m*NA + k )
				QT(m*NA + k ) =-s * QTnk + c*QT(m*NA + k )
			enddo
		enddo
	enddo
	
	end

