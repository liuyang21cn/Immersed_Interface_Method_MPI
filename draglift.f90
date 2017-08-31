SUBROUTINE  draglift

  USE MPI
  USE para
  USE Euler
  USE stretchxy
  USE Lagrange

  INTEGER:: ic,jc,id,jd
  INTEGER,PARAMETER:: many=3,nany=3
  DOUBLE PRECISION:: ds,xn,yn,gacobi,foo,pplus,dist
  DOUBLE PRECISION,DIMENSION(0:many):: xx,yy
  DOUBLE PRECISION,DIMENSION(1:many):: xa,ya,xb,yb,pp

  DOUBLE PRECISION,DIMENSION(nvertex4proc(nobj4proc)):: fx,fy
  INTEGER:: myvertex,myobj,myobj_local,myobj_global
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE:: cd_temp,cl_temp
  DOUBLE PRECISION, DIMENSION(1:2):: A,corner0,corner1

  ALLOCATE(cd_temp(1:nobj))
  ALLOCATE(cl_temp(1:nobj))

  ! initialize cd cl to zero
  IF(myid==0) THEN
     cd =zero
     cl =zero
  ENDIF
  cd_temp=zero
  cl_temp=zero
  fx=zero
  fy=zero

  CALL mpi_barrier(comm2d,ierr)
  ds=1.01d0*SQRT(dxmax*dxmax+dymax*dymax)
  DO myobj=1,nobj4proc
     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)
        A(1) = xs(myvertex)
        A(2) = ys(myvertex)
        corner0=(/xf(iubeg-1),yf(jvbeg-1)/)
        corner1=(/xf(iuend+1),yf(jvend+1)/)
        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN
           xn= taoyn(myvertex)
           yn=-taoxn(myvertex)
           gacobi=SQRT(xn*xn+yn*yn)
           xn=xn/gacobi
           yn=yn/gacobi
           xx(0)=xs(myvertex)
           yy(0)=ys(myvertex)

           DO n=1,nany
              xx(n)= xx(0)+DBLE(n)*ds*xn
              yy(n)= yy(0)+DBLE(n)*ds*yn

              i=INT((x2xi(xx(n))-xi0)/dxi*2)
              j=INT((y2eta(yy(n))-eta0)/deta*2)
              IF(MOD(i,2).EQ.0) THEN
                 ic=i/2
              ELSE
                 ic=(i+1)/2
              ENDIF
              IF(MOD(j,2).EQ.0) THEN
                 jc=j/2
              ELSE
                 jc=(j+1)/2
              ENDIF

              id=INT(SIGN(1.0d0,xn))
              jd=INT(SIGN(1.0d0,yn))
              IF(id.LT.0) THEN
                 ic=ic+1
              ENDIF
              IF(jd.LT.0) THEN
                 jc=jc+1
              ENDIF

              ! interpolate pp(n)
              DO i=1,many
                 DO j=1,many
                    xa(j)=yc(jc+jd*(j-1))
                    ya(j)=p(ic+id*(i-1),jc+jd*(j-1))
                 ENDDO
                 xb(i)=xc(ic+id*(i-1))
                 CALL interpolate(xa,ya,many,yy(n),yb(i),foo,myid,myvertex)
              ENDDO
              CALL interpolate(xb,yb,many,xx(n),pp(n),foo,myid,myvertex)
           ENDDO !enddo n=1,nany
           pplus=3.0d0*pp(1)-3.0d0*pp(2)+pp(3)
           fx(myvertex)=-pplus*xn+Re1*ujcn(myvertex)
           fy(myvertex)=-pplus*yn+Re1*vjcn(myvertex)
        ENDIF
     ENDDO

     DO myvertex=nvertex4proc(myobj),nvertex4proc(myobj-1)+1,-1
        A(1) = xs(myvertex)
        A(2) = ys(myvertex)
        corner0=(/xf(iubeg-1),yf(jvbeg-1)/)
        corner1=(/xf(iuend+1),yf(jvend+1)/)
        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN
           IF ( myvertex ==  nvertex4proc(myobj-1)+1) THEN
              tx = -taoxn(nvertex4proc(myobj)-1)
              ty = -taoyn(nvertex4proc(myobj)-1)
           ELSE
              tx = -taoxn(myvertex-1)
              ty = -taoyn(myvertex-1)
           END IF
           gacobi=SQRT(tx*tx+ty*ty)
           gacobi2=tx*tx+ty*ty

           ! normal direction
           xn=-ty
           yn= tx
           xn=xn/gacobi
           yn=yn/gacobi

           xx(0)=xs(myvertex)
           yy(0)=ys(myvertex)
           DO n=1,nany
              xx(n)= xx(0)+DBLE(n)*ds*xn
              yy(n)= yy(0)+DBLE(n)*ds*yn
              i=INT((x2xi(xx(n))-xi0)/dxi*2)
              j=INT((y2eta(yy(n))-eta0)/deta*2)
              IF(MOD(i,2).EQ.0) THEN
                 ic=i/2
              ELSE
                 ic=(i+1)/2
              ENDIF
              IF(MOD(j,2).EQ.0) THEN
                 jc=j/2
              ELSE
                 jc=(j+1)/2
              ENDIF

              id=INT(SIGN(1.0d0,xn))
              jd=INT(SIGN(1.0d0,yn))
              IF(id.LT.0) THEN
                 ic=ic+1
              ENDIF
              IF(jd.LT.0) THEN
                 jc=jc+1
              ENDIF

              ! interpolate pp(n)
              DO i=1,many
                 DO j=1,many
                    xa(j)=yc(jc+jd*(j-1))
                    ya(j)=p(ic+id*(i-1),jc+jd*(j-1))
                 ENDDO
                 xb(i)=xc(ic+id*(i-1))
                 CALL interpolate(xa,ya,many,yy(n),yb(i),foo,myid,myvertex)
              ENDDO
              CALL interpolate(xb,yb,many,xx(n),pp(n),foo,myid,myvertex)
           ENDDO               !enddo n=1,nany
           pplus=3.0d0*pp(1)-3.0d0*pp(2)+pp(3)
           fx(myvertex)=0.5d0*(fx(myvertex)+(-pplus*xn+Re1*ujcn(myvertex)))
           fy(myvertex)=0.5d0*(fy(myvertex)+(-pplus*yn+Re1*vjcn(myvertex)))
        ENDIF
     ENDDO

     ! calculate local cd & cl
     DO myobj_global = 1, nobj
        IF ( obj4proc(myobj) == object_list(myobj_global) ) THEN
           DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
              IF(xs(myvertex).LE.xc(ipend+1) .AND. xs(myvertex).GT.xc(ipbeg)) THEN
                 IF(ys(myvertex).LE.yc(jpend+1) .AND. ys(myvertex).GT.yc(jpbeg)) THEN
                    dist=SQRT((xs(myvertex+1)-xs(myvertex))*(xs(myvertex+1)-xs(myvertex))+ &
                         (ys(myvertex+1)-ys(myvertex))*(ys(myvertex+1)-ys(myvertex)))
                    cd_temp(myobj_global)=cd_temp(myobj_global)+(fx(myvertex)+fx(myvertex+1))*dist
                    cl_temp(myobj_global)=cl_temp(myobj_global)+(fy(myvertex)+fy(myvertex+1))*dist
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO

  ! match local object to global object and reduce cd,cl
  DO myobj = 1, nobj
     CALL MPI_reduce(cd_temp(myobj),cd(myobj),1, MPI_DOUBLE_PRECISION,MPI_SUM,0,comm2d,ierr)
     CALL MPI_reduce(cl_temp(myobj),cl(myobj),1, MPI_DOUBLE_PRECISION,MPI_SUM,0,comm2d,ierr)
  ENDDO

  DEALLOCATE(cd_temp)
  DEALLOCATE(cl_temp)
END SUBROUTINE draglift
