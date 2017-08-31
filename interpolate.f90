SUBROUTINE interpolate(xa,ya,n,x,y,dy,myid,myvertex)

  INTEGER:: n,myid,myvertex
  DOUBLE PRECISION:: dy,x,y
  INTEGER,PARAMETER:: nmax=10
  DOUBLE PRECISION, DIMENSION(1:n):: xa,ya
  INTEGER:: i,m,ns
  DOUBLE PRECISION:: den,dif,dift,ho,hp,w
  DOUBLE PRECISION, DIMENSION(1:nmax):: c,d

  IF(n.GT.nmax) THEN
     WRITE(*,*)'  !! interpolation polynomial order too high!'
     STOP
  ENDIF

  ns=1
  dif=ABS(x-xa(1))

  DO i=1,n
     dift=ABS(x-xa(i))
     IF(dift.LT.dif) THEN
        ns=i
        dif=dift
     ENDIF
     c(i)=ya(i)
     d(i)=ya(i)
  ENDDO

  y=ya(ns)
  ns=ns-1

  DO m=1,n-1
     DO i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        IF(den.EQ.0.0d0) THEN
           !  PRINT*, 'intern',xa(i),xa(i+m)
           PRINT*, 'intern',xa(1),xa(2),xa(3),myid,myvertex
           WRITE(*,*)'failure in interpolation!'
           STOP
        ENDIF
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
     ENDDO
     IF(2*ns.LT.n-m) THEN
        dy=c(ns+1)
     ELSE
        dy=d(ns)
        ns=ns-1
     ENDIF
     y=y+dy
  ENDDO

END SUBROUTINE interpolate
