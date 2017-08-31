SUBROUTINE jc_firstsecond(krk)

  USE mpi
  USE para
  USE stretchxy
  USE Euler
  USE Lagrange

  INTEGER:: myvertex,myobj,krk

  DOUBLE PRECISION:: gacobi2,gacobi,tx,ty,xn,yn
  DOUBLE PRECISION:: duxjc,duyjc,dduxjc,dduyjc,dduxyjc,ddpxjc0,ddpxjc1
  DOUBLE PRECISION:: dvxjc,dvyjc,ddvxjc,ddvyjc,ddvxyjc,ddpyjc0,ddpyjc1
  DOUBLE PRECISION:: dpxjc,dpyjc
  DOUBLE PRECISION:: dudnjc,dudnjcn,dudnjcp,dvdnjc,dvdnjcn,dvdnjcp
  DOUBLE PRECISION:: ddudnp,ddvdnp,qx,qy,xs0,ys0,xs1,ys1

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: deltaujc1,deltavjc1,deltaujc2,deltavjc2
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE:: ujc2,vjc2,pjc2,dunjcp,dvnjcp

  ! dudxdt: derivative of dudx over tao
  DOUBLE PRECISION:: dudxdt,dvdxdt,dpdxdt,dudydt,dvdydt,dpdydt

  ! variables used in interpolation of dudnjc ddudnjc
  DOUBLE PRECISION:: ds,dist,foo,r3
  DOUBLE PRECISION,DIMENSION(0:3):: uu,vv,xx,yy,distance
  INTEGER:: ie,ic,je,jc,id,jd,iu,ju,iv,jv,io,jo
  INTEGER, PARAMETER:: many=3,nany=3
  DOUBLE PRECISION,DIMENSION(1:many):: xa,ya,xb,yb
  DOUBLE PRECISION, DIMENSION(1:2):: A,corner0,corner1

  ALLOCATE(deltaujc1(1:nvertex4proc(nobj4proc)))
  ALLOCATE(deltavjc1(1:nvertex4proc(nobj4proc)))
  ALLOCATE(deltaujc2(1:nvertex4proc(nobj4proc)))
  ALLOCATE(deltavjc2(1:nvertex4proc(nobj4proc)))

  ALLOCATE(ujc2(1:6, 1:nvertex4proc(nobj4proc)))
  ALLOCATE(vjc2(1:6, 1:nvertex4proc(nobj4proc)))
  ALLOCATE(pjc2(1:8, 1:nvertex4proc(nobj4proc)))

  ujc = zero
  vjc = zero
  pjc = zero
  ujc2= zero
  vjc2= zero
  pjc2= zero

  uu = zero
  vv = zero
  xx = zero
  yy = zero
  distance = zero
  IF(krk.EQ.3 .AND. JCsaved.EQV..TRUE.) THEN
     ujcn=zero
     vjcn=zero
     taoxn=zero
     taoyn=zero
  ENDIF
  ds=1.01d0*SQRT(dxmax*dxmax+dymax*dymax)
  ! counter clockwise direction
  DO myobj=1,nobj4proc
     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)
        A(1) = xs(myvertex)
        A(2) = ys(myvertex)
        corner0=(/xf(iubeg-1),yf(jvbeg-1)/)
        corner1=(/xf(iuend+1),yf(jvend+1)/)
        corner0=(/xf(iubeg-3),yf(jvbeg-3)/)
        corner1=(/xf(iuend+3),yf(jvend+3)/)
        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN
           ! tao direction
           tx=taox(myvertex)
           ty=taoy(myvertex)

           gacobi=SQRT(tx*tx+ty*ty)
           gacobi2=tx*tx+ty*ty

           ! normal direction
           xn= ty
           yn=-tx
           xn=xn/gacobi
           yn=yn/gacobi

           ! velocity and position of vertices
           uu(0)=us(myvertex)
           vv(0)=vs(myvertex)
           xx(0)=xs(myvertex)
           yy(0)=ys(myvertex)
           ! dudn jump conditions on negative side
           ! dudnjcn= thetat X n
           dudnjcn = -thetat(myobj)*yn
           dvdnjcn =  thetat(myobj)*xn

           ! dudn jump conditions on positive side
           ! ============ put this part into a subroutine ============
           ! a. find 3 points on normal direction s1, s2, s3
           DO n=1,nany
              ! coordinates of interpolation point
              xx(n)= xx(0)+DBLE(n)*ds*xn
              yy(n)= yy(0)+DBLE(n)*ds*yn

              ! indices of interpolation point
              i=INT((x2xi(xx(n))-xi0)/dxi*2)
              j=INT((y2eta(yy(n))-eta0)/deta*2)
              IF(MOD(i,2)==0) THEN
                 ie=i/2
                 ic=ie
              ELSE
                 ic=(i+1)/2
                 ie=ic-1
              ENDIF
              IF(MOD(j,2)==0) THEN
                 je=j/2
                 jc=je
              ELSE
                 jc=(j+1)/2
                 je=jc-1
              ENDIF
              iu=ie
              ju=jc
              iv=ic
              jv=je

              id=INT(SIGN(1.d0,xn))
              jd=INT(SIGN(1.d0,yn))
              IF(id<0.0d0) THEN
                 iu=iu+1
                 iv=iv+1
              ENDIF
              IF(jd<0.0d0) THEN
                 ju=ju+1
                 jv=jv+1
              ENDIF

              ! interpolation of uu(n)
              DO j=1,many
                 DO i=1,many
                    xa(i)=xf(iu+id*(i-1))
                    ya(i)=u(iu+id*(i-1),ju+jd*(j-1))
                 ENDDO
                 xb(j)=yc(ju+jd*(j-1))
                 CALL interpolate(xa,ya,many,xx(n),yb(j),foo,myid,myvertex)
              ENDDO
              CALL interpolate(xb,yb,many,yy(n),uu(n),foo,myid,myvertex)

              ! interpolation of vv(n)
              DO i=1,many
                 DO j=1,many
                    xa(j)=yf(jv+jd*(j-1))
                    ya(j)=v(iv+id*(i-1),jv+jd*(j-1))
                 ENDDO
                 xb(i)=xc(iv+id*(i-1))
                 CALL interpolate(xa,ya,many,yy(n),yb(i),foo,myid,myvertex)
              ENDDO
              CALL interpolate(xb,yb,many,xx(n),vv(n),foo,myid,myvertex)
           ENDDO

           ! b. dudn jump conditions on positive side
           dudnjcp=(4.0d0*uu(1)-uu(2)-3.0d0*uu(0))/ds/2.0d0
           dvdnjcp=(4.0d0*vv(1)-vv(2)-3.0d0*vv(0))/ds/2.0d0

           IF(krk.EQ.3 .AND. JCsaved.EQV..TRUE.) THEN
              ujcn(myvertex)=dudnjcp
              vjcn(myvertex)=dvdnjcp
              taoxn(myvertex)=taox(myvertex)
              taoyn(myvertex)=taoy(myvertex)
           ENDIF

           ddudnp=(-5.0d0*uu(1)+4.0d0*uu(2)-1.0d0*uu(3)+2.0d0*uu(0))/ds/ds
           ddvdnp=(-5.0d0*vv(1)+4.0d0*vv(2)-1.0d0*vv(3)+2.0d0*vv(0))/ds/ds

           ! c. dudn jump conditions
           dudnjc = dudnjcp - dudnjcn
           dvdnjc = dvdnjcp - dvdnjcn

           ! laplacian u jump conditions
           ! ddudn jump conditions
           deltaujc1(myvertex) = ddudnp + vertex(3,myvertex)*dudnjc
           deltavjc1(myvertex) = ddvdnp + vertex(3,myvertex)*dvdnjc

           ! dudx dvdx dudy dvdy
           duxjc = -ty/(tx*yn-ty*xn)*dudnjc
           dvxjc = -ty/(tx*yn-ty*xn)*dvdnjc
           duyjc =  tx/(tx*yn-ty*xn)*dudnjc
           dvyjc =  tx/(tx*yn-ty*xn)*dvdnjc

           ! dpdx dpdy jump conditions
           ! body force qx, qy
           qx =  thetatt(myobj)*(ys(myvertex)-ysc(myobj))
           qy = -thetatt(myobj)*(xs(myvertex)-xsc(myobj))

           ! dpdx, dpdy jump conditions
           dpxjc = 1.0d0/Re*deltaujc1(myvertex) + qx
           dpyjc = 1.0d0/Re*deltavjc1(myvertex) + qy

           ujc(1,myvertex)=duxjc
           vjc(1,myvertex)=dvxjc
           pjc(1,myvertex)=dpxjc
           ujc(2,myvertex)=duyjc
           vjc(2,myvertex)=dvyjc
           pjc(2,myvertex)=dpyjc

        ENDIF
     ENDDO ! end myvertex

     ! second order jump conditions
     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
        A(1)=xs(myvertex)
        A(2)=ys(myvertex)
        corner0=allcorner0(:,myid)
        corner1=allcorner1(:,myid)
        corner0=(/xf(iubeg-1),yf(jvbeg-1)/)
        corner1=(/xf(iuend+1),yf(jvend+1)/)

        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN

           ! tao direction
           tx=taox(myvertex)
           ty=taoy(myvertex)
           gacobi2=tx*tx+ty*ty

           ! dudxdt, dvdxdt, dudydt, dxdydt
           ! dudxdt = (dudxB-dudxA)/|AB|
           ! position of neighbor vertex
           xs0=xs(myvertex)
           xs1=xs(myvertex+1)
           ys0=ys(myvertex)
           ys1=ys(myvertex+1)
           dist=SQRT((ys1-ys0)*(ys1-ys0)+(xs1-xs0)*(xs1-xs0))

           dudxdt = (ujc(1,myvertex+1) - ujc(1,myvertex))/dist
           dvdxdt = (vjc(1,myvertex+1) - vjc(1,myvertex))/dist
           dpdxdt = (pjc(1,myvertex+1) - pjc(1,myvertex))/dist

           dudydt = (ujc(2,myvertex+1) - ujc(2,myvertex))/dist
           dvdydt = (vjc(2,myvertex+1) - vjc(2,myvertex))/dist
           dpdydt = (pjc(2,myvertex+1) - pjc(2,myvertex))/dist

           ! dduxjc dduyjc dduxyjc
           dduxjc=(tx*dudxdt-ty*dudydt+ty*ty*deltaujc1(myvertex))/gacobi2
           ddvxjc=(tx*dvdxdt-ty*dvdydt+ty*ty*deltavjc1(myvertex))/gacobi2

           dduyjc=(ty*dudydt-tx*dudxdt+tx*tx*deltaujc1(myvertex))/gacobi2
           ddvyjc=(ty*dvdydt-tx*dvdxdt+tx*tx*deltavjc1(myvertex))/gacobi2

           dduxyjc=(ty*dudxdt+tx*dudydt-tx*ty*deltaujc1(myvertex))/gacobi2
           ddvxyjc=(ty*dvdxdt+tx*dvdydt-tx*ty*deltavjc1(myvertex))/gacobi2

           ! ddpdx ddpdy jump conditions
           r3 = 0.0d0
           ddpxjc0  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc0  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

           r3 = 1.0d0
           ddpxjc1  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc1  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

           ujc(3,myvertex)=dduxjc
           vjc(3,myvertex)=ddvxjc
           pjc(3,myvertex)=ddpxjc1
           ujc(5,myvertex)=dduxyjc
           vjc(5,myvertex)=ddvxyjc
           pjc(7,myvertex)=ddpxjc0

           ujc(4,myvertex)=dduyjc
           vjc(4,myvertex)=ddvyjc
           pjc(4,myvertex)=ddpyjc1
           ujc(6,myvertex)=dduxyjc
           vjc(6,myvertex)=ddvxyjc
           pjc(8,myvertex)=ddpyjc0
        ENDIF
     ENDDO !end my vertex
     ujc(3,nvertex4proc(myobj)) = ujc(3,nvertex4proc(myobj-1)+1)
     vjc(3,nvertex4proc(myobj)) = vjc(3,nvertex4proc(myobj-1)+1)
     pjc(3,nvertex4proc(myobj)) = pjc(3,nvertex4proc(myobj-1)+1)
     ujc(5,nvertex4proc(myobj)) = ujc(5,nvertex4proc(myobj-1)+1)
     vjc(5,nvertex4proc(myobj)) = vjc(5,nvertex4proc(myobj-1)+1)
     pjc(7,nvertex4proc(myobj)) = pjc(7,nvertex4proc(myobj-1)+1)

     ujc(4,nvertex4proc(myobj)) = ujc(4,nvertex4proc(myobj-1)+1)
     vjc(4,nvertex4proc(myobj)) = vjc(4,nvertex4proc(myobj-1)+1)
     pjc(4,nvertex4proc(myobj)) = pjc(4,nvertex4proc(myobj-1)+1)
     ujc(6,nvertex4proc(myobj)) = ujc(6,nvertex4proc(myobj-1)+1)
     vjc(6,nvertex4proc(myobj)) = vjc(6,nvertex4proc(myobj-1)+1)
     pjc(8,nvertex4proc(myobj)) = pjc(8,nvertex4proc(myobj-1)+1)

  ENDDO ! end myobj

  ! clockwise direction
  DO myobj=1,nobj4proc
     DO myvertex=nvertex4proc(myobj),nvertex4proc(myobj-1)+1,-1
        A(1) = xs(myvertex)
        A(2) = ys(myvertex)
        corner0=(/xf(iubeg-1),yf(jvbeg-1)/)
        corner1=(/xf(iuend+1),yf(jvend+1)/)
        corner0=(/xf(iubeg-3),yf(jvbeg-3)/)
        corner1=(/xf(iuend+3),yf(jvend+3)/)
        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN
           ! tao direction
           IF ( myvertex ==  nvertex4proc(myobj-1)+1) THEN
              tx = -taox(nvertex4proc(myobj)-1)
              ty = -taoy(nvertex4proc(myobj)-1)
           ELSE
              tx = -taox(myvertex-1)
              ty = -taoy(myvertex-1)
           END IF

           gacobi=SQRT(tx*tx+ty*ty)
           gacobi2=tx*tx+ty*ty

           ! normal direction
           xn=-ty
           yn= tx
           xn=xn/gacobi
           yn=yn/gacobi

           ! velocity and position of vertices
           uu(0)=us(myvertex)
           vv(0)=vs(myvertex)
           xx(0)=xs(myvertex)
           yy(0)=ys(myvertex)
           distance(0)=zero

           ! dudn jump conditions on negative side
           ! dudnjcn= thetat X n

           dudnjcn = -thetat(myobj)*yn
           dvdnjcn =  thetat(myobj)*xn

           ! dudn jump conditions on positive side
           ! ============ put this part into a subroutine ============
           ! a. find 3 points on normal direction s1, s2, s3

           DO n=1,nany
              ! coordinates of interpolation point
              xx(n)= xx(0)+DBLE(n)*ds*xn
              yy(n)= yy(0)+DBLE(n)*ds*yn

              ! indices of interpolation point
              i=INT((x2xi(xx(n))-xi0)/dxi*2)
              j=INT((y2eta(yy(n))-eta0)/deta*2)
              IF(MOD(i,2)==0) THEN
                 ie=i/2
                 ic=ie
              ELSE
                 ic=(i+1)/2
                 ie=ic-1
              ENDIF
              IF(MOD(j,2)==0) THEN
                 je=j/2
                 jc=je
              ELSE
                 jc=(j+1)/2
                 je=jc-1
              ENDIF
              iu=ie
              ju=jc
              iv=ic
              jv=je

              id=INT(SIGN(1.d0,xn))
              jd=INT(SIGN(1.d0,yn))
              IF(id<0.0d0) THEN
                 iu=iu+1
                 iv=iv+1
              ENDIF
              IF(jd<0.0d0) THEN
                 ju=ju+1
                 jv=jv+1
              ENDIF

              ! interpolation of uu(n)
              DO j=1,many
                 DO i=1,many
                    xa(i)=xf(iu+id*(i-1))
                    ya(i)=u(iu+id*(i-1),ju+jd*(j-1))
                 ENDDO
                 xb(j)=yc(ju+jd*(j-1))
                 CALL interpolate(xa,ya,many,xx(n),yb(j),foo,myid,myvertex)
              ENDDO
              CALL interpolate(xb,yb,many,yy(n),uu(n),foo,myid,myvertex)

              ! interpolation of vv(n)
              DO i=1,many
                 DO j=1,many
                    xa(j)=yf(jv+jd*(j-1))
                    ya(j)=v(iv+id*(i-1),jv+jd*(j-1))
                 ENDDO
                 xb(i)=xc(iv+id*(i-1))
                 CALL interpolate(xa,ya,many,yy(n),yb(i),foo,myid,myvertex)
              ENDDO
              CALL interpolate(xb,yb,many,xx(n),vv(n),foo,myid,myvertex)
           ENDDO

           ! b. dudn jump conditions on positive side
           dudnjcp=(4.0d0*uu(1)-uu(2)-3.0d0*uu(0))/ds/2.0d0
           dvdnjcp=(4.0d0*vv(1)-vv(2)-3.0d0*vv(0))/ds/2.0d0

           ddudnp=(-5.0d0*uu(1)+4.0d0*uu(2)-1.0d0*uu(3)+2.0d0*uu(0))/ds/ds
           ddvdnp=(-5.0d0*vv(1)+4.0d0*vv(2)-1.0d0*vv(3)+2.0d0*vv(0))/ds/ds

           ! save dudn|+ for draglift
           IF(krk.EQ.3 .AND. JCsaved.EQV..TRUE.) THEN
              ujcn(myvertex)=0.5d0*(dudnjcp+ujcn(myvertex))
              vjcn(myvertex)=0.5d0*(dvdnjcp+vjcn(myvertex))
           ENDIF

           ! c. dudn jump conditions
           dudnjc = dudnjcp - dudnjcn
           dvdnjc = dvdnjcp - dvdnjcn

           ! laplacian u jump conditions
           ! ddudn jump conditions
           deltaujc2(myvertex) = ddudnp + vertex(3,myvertex)*dudnjc
           deltavjc2(myvertex) = ddvdnp + vertex(3,myvertex)*dvdnjc

           ! dudx dvdx dudy dvdy
           duxjc = -ty/(tx*yn-ty*xn)*dudnjc
           dvxjc = -ty/(tx*yn-ty*xn)*dvdnjc
           duyjc =  tx/(tx*yn-ty*xn)*dudnjc
           dvyjc =  tx/(tx*yn-ty*xn)*dvdnjc

           ! dpdx dpdy jump conditions
           ! body force qx, qy
           qx =  thetatt(myobj)*(ys(myvertex)-ysc(myobj))
           qy = -thetatt(myobj)*(xs(myvertex)-xsc(myobj))

           ! dpdx, dpdy jump conditions
           dpxjc = 1.0d0/Re*deltaujc2(myvertex) + qx
           dpyjc = 1.0d0/Re*deltavjc2(myvertex) + qy

           ujc2(1,myvertex)=duxjc
           vjc2(1,myvertex)=dvxjc
           pjc2(1,myvertex)=dpxjc

           ujc2(2,myvertex)=duyjc
           vjc2(2,myvertex)=dvyjc
           pjc2(2,myvertex)=dpyjc

        ENDIF
     ENDDO

     ! second order jump conditions
     DO myvertex=nvertex4proc(myobj),nvertex4proc(myobj-1)+2,-1
        A(1)=xs(myvertex)
        A(2)=ys(myvertex)
        corner0=allcorner0(:,myid)
        corner1=allcorner1(:,myid)
        corner0=(/xf(iubeg-1),yf(jvbeg-1)/)
        corner1=(/xf(iuend+1),yf(jvend+1)/)

        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN

           ! tao direction
           tx = -taox(myvertex-1)
           ty = -taoy(myvertex-1)
           gacobi2=tx*tx+ty*ty

           ! dudxdt, dvdxdt, dudydt, dxdydt
           ! dudxdt = (dudxB-dudxA)/|AB|
           ! position of neighbor vertex
           xs0=xs(myvertex)
           xs1=xs(myvertex-1)
           ys0=ys(myvertex)
           ys1=ys(myvertex-1)
           dist=SQRT((ys1-ys0)*(ys1-ys0)+(xs1-xs0)*(xs1-xs0))

           ! if not the last vertex, just next minus current
           dudxdt = (ujc2(1,myvertex-1) - ujc2(1,myvertex))/dist
           dvdxdt = (vjc2(1,myvertex-1) - vjc2(1,myvertex))/dist
           dpdxdt = (pjc2(1,myvertex-1) - pjc2(1,myvertex))/dist

           dudydt = (ujc2(2,myvertex-1) - ujc2(2,myvertex))/dist
           dvdydt = (vjc2(2,myvertex-1) - vjc2(2,myvertex))/dist
           dpdydt = (pjc2(2,myvertex-1) - pjc2(2,myvertex))/dist

           ! dduxjc dduyjc dduxyjc
           dduxjc=(tx*dudxdt-ty*dudydt+ty*ty*deltaujc2(myvertex))/gacobi2
           ddvxjc=(tx*dvdxdt-ty*dvdydt+ty*ty*deltavjc2(myvertex))/gacobi2

           dduyjc=(ty*dudydt-tx*dudxdt+tx*tx*deltaujc2(myvertex))/gacobi2
           ddvyjc=(ty*dvdydt-tx*dvdxdt+tx*tx*deltavjc2(myvertex))/gacobi2

           dduxyjc=(ty*dudxdt+tx*dudydt-tx*ty*deltaujc2(myvertex))/gacobi2
           ddvxyjc=(ty*dvdxdt+tx*dvdydt-tx*ty*deltavjc2(myvertex))/gacobi2

           ! ddpdx ddpdy jump conditions
           r3 = 0.0d0
           ddpxjc0  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc0  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

           r3 = 1.0d0
           ddpxjc1  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc1  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

           ujc2(3,myvertex)=dduxjc
           vjc2(3,myvertex)=ddvxjc
           pjc2(3,myvertex)=ddpxjc1
           ujc2(5,myvertex)=dduxyjc
           vjc2(5,myvertex)=ddvxyjc
           pjc2(7,myvertex)=ddpxjc0

           ujc2(4,myvertex)=dduyjc
           vjc2(4,myvertex)=ddvyjc
           pjc2(4,myvertex)=ddpyjc1
           ujc2(6,myvertex)=dduxyjc
           vjc2(6,myvertex)=ddvxyjc
           pjc2(8,myvertex)=ddpyjc0

        ENDIF
     ENDDO ! end myvertex
     ujc2(3,nvertex4proc(myobj-1)+1) = ujc2(3,nvertex4proc(myobj))
     vjc2(3,nvertex4proc(myobj-1)+1) = vjc2(3,nvertex4proc(myobj))
     pjc2(3,nvertex4proc(myobj-1)+1) = pjc2(3,nvertex4proc(myobj))
     ujc2(5,nvertex4proc(myobj-1)+1) = ujc2(5,nvertex4proc(myobj))
     vjc2(5,nvertex4proc(myobj-1)+1) = vjc2(5,nvertex4proc(myobj))
     pjc2(7,nvertex4proc(myobj-1)+1) = pjc2(7,nvertex4proc(myobj))

     ujc2(4,nvertex4proc(myobj-1)+1) = ujc2(4,nvertex4proc(myobj))
     vjc2(4,nvertex4proc(myobj-1)+1) = vjc2(4,nvertex4proc(myobj))
     pjc2(4,nvertex4proc(myobj-1)+1) = pjc2(4,nvertex4proc(myobj))
     ujc2(6,nvertex4proc(myobj-1)+1) = ujc2(6,nvertex4proc(myobj))
     vjc2(6,nvertex4proc(myobj-1)+1) = vjc2(6,nvertex4proc(myobj))
     pjc2(8,nvertex4proc(myobj-1)+1) = pjc2(8,nvertex4proc(myobj))

  ENDDO ! end myobj

  DO myobj=1,nobj4proc
     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)

        ujc(1,myvertex)=(ujc(1,myvertex)+ujc2(1,myvertex))/2.0d0
        ujc(2,myvertex)=(ujc(2,myvertex)+ujc2(2,myvertex))/2.0d0
        ujc(3,myvertex)=(ujc(3,myvertex)+ujc2(3,myvertex))/2.0d0
        ujc(4,myvertex)=(ujc(4,myvertex)+ujc2(4,myvertex))/2.0d0
        ujc(5,myvertex)=(ujc(5,myvertex)+ujc2(5,myvertex))/2.0d0
        ujc(6,myvertex)=(ujc(6,myvertex)+ujc2(6,myvertex))/2.0d0

        vjc(1,myvertex)=(vjc(1,myvertex)+vjc2(1,myvertex))/2.0d0
        vjc(2,myvertex)=(vjc(2,myvertex)+vjc2(2,myvertex))/2.0d0
        vjc(3,myvertex)=(vjc(3,myvertex)+vjc2(3,myvertex))/2.0d0
        vjc(4,myvertex)=(vjc(4,myvertex)+vjc2(4,myvertex))/2.0d0
        vjc(5,myvertex)=(vjc(5,myvertex)+vjc2(5,myvertex))/2.0d0
        vjc(6,myvertex)=(vjc(6,myvertex)+vjc2(6,myvertex))/2.0d0

        pjc(1,myvertex)=(pjc(1,myvertex)+pjc2(1,myvertex))/2.0d0
        pjc(2,myvertex)=(pjc(2,myvertex)+pjc2(2,myvertex))/2.0d0
        pjc(3,myvertex)=(pjc(3,myvertex)+pjc2(3,myvertex))/2.0d0
        pjc(4,myvertex)=(pjc(4,myvertex)+pjc2(4,myvertex))/2.0d0
        pjc(5,myvertex)=(pjc(5,myvertex)+pjc2(5,myvertex))/2.0d0
        pjc(6,myvertex)=(pjc(6,myvertex)+pjc2(6,myvertex))/2.0d0
        pjc(7,myvertex)=(pjc(7,myvertex)+pjc2(7,myvertex))/2.0d0
        pjc(8,myvertex)=(pjc(8,myvertex)+pjc2(8,myvertex))/2.0d0
     ENDDO
  ENDDO

  DEALLOCATE(deltaujc1)
  DEALLOCATE(deltavjc1)
  DEALLOCATE(deltaujc2)
  DEALLOCATE(deltavjc2)

  DEALLOCATE(ujc2)
  DEALLOCATE(vjc2)
  DEALLOCATE(pjc2)



END SUBROUTINE jc_firstsecond
