SUBROUTINE principle_jc_p

  USE mpi
  USE para
  USE stretchxy
  USE Euler
  USE Lagrange

  ! create group structure for later mpi_group use
  TYPE cg_group
     INTEGER:: group, comm
     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: r_temp
     INTEGER, DIMENSION(:), ALLOCATABLE:: p_rank, req
     INTEGER, DIMENSION(:,:), ALLOCATABLE:: STATUS
  END TYPE cg_group

  TYPE(cg_group):: group_mpi(nobj)
  INTEGER:: myvertex_left,myvertex_right,myobj,group0,myobj_local
  INTEGER:: rank_count, rank_total(1:nobj), rank_temp(0:nprocs-1),rank_c(1:nobj)
  DOUBLE PRECISION, DIMENSION(1:2):: A,B,corner0,corner1

  ! variables for princilple jump conditions
  DOUBLE PRECISION:: tx,ty,tx_left,ty_left,tx_right,ty_right,xn,yn,gacobi
  DOUBLE PRECISION:: dudnjcp,dudnjcn,dudnjc,dvdnjcp,dvdnjcn,dvdnjc
  DOUBLE PRECISION:: ddudnp,ddvdnp,deltaujc,deltavjc,qx,qy,temp1,temp2,temp

  ! variables for interpolation
  DOUBLE PRECISION:: xx(0:3),yy(0:3),uu(0:3),vv(0:3),dist,foo
  INTEGER:: ie,ic,je,jc,id,jd,iu,ju,iv,jv,io,jo
  INTEGER, PARAMETER:: many=3,nany=3
  DOUBLE PRECISION:: xa(many),ya(many),xb(many),yb(many)

  ! Gaussian Elimination solver
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: aa, bb, cc
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: x_sol

  double PRECISION:: time0,time1

  ALLOCATE(dpxjcM(1:nvertex4proc(nobj4proc)))
  ALLOCATE(dpyjcM(1:nvertex4proc(nobj4proc)))
  ALLOCATE(jcp_rhs(1:nvertex4proc(nobj4proc)))

  dpxjcM = zero
  dpyjcM = zero
  jcp_rhs = zero

  ds=1.01d0*SQRT(dxmax*dxmax+dymax*dymax)

  DO myobj=1,nobj4proc
     !  1. compute jump conditions at middle points
     DO myvertex = nvertex4proc(myobj-1)+1, nvertex4proc(myobj)-1
        A(1) = xs(myvertex)
        A(2) = ys(myvertex)
        corner0=(/xf(iubeg-1),yf(jvbeg-1)/)
        corner1=(/xf(iuend+1),yf(jvend+1)/)
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
           uu(0)=(us(myvertex)+us(myvertex+1))/2.0d0
           vv(0)=(vs(myvertex)+vs(myvertex+1))/2.0d0
           xx(0)=(xs(myvertex)+xs(myvertex+1))/2.0d0
           yy(0)=(ys(myvertex)+ys(myvertex+1))/2.0d0

           ! dudn jump conditions on negative side
           dudnjcn = -thetat(myobj)*yn
           dvdnjcn =  thetat(myobj)*xn

           ! find 3 points on normal direction s1, s2, s3
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
              END IF
              IF(MOD(j,2)==0) THEN
                 je=j/2
                 jc=je
              ELSE
                 jc=(j+1)/2
                 je=jc-1
              END IF
              iu=ie
              ju=jc
              iv=ic
              jv=je

              id=INT(SIGN(1.d0,xn))
              jd=INT(SIGN(1.d0,yn))
              IF(id<0.0d0) THEN
                 iu=iu+1
                 iv=iv+1
              END IF
              IF(jd<0.0d0) THEN
                 ju=ju+1
                 jv=jv+1
              END IF

              ! interpolation of uu(n)
              DO j=1,many
                 DO i=1,many
                    xa(i)=xf(iu+id*(i-1))
                    ya(i)=u(iu+id*(i-1),ju+jd*(j-1))
                 END DO
                 xb(j)=yc(ju+jd*(j-1))
                 CALL interpolate(xa,ya,many,xx(n),yb(j),foo,myid,myvertex)
              END DO
              CALL interpolate(xb,yb,many,yy(n),uu(n),foo,myid,myvertex)

              ! interpolation of vv(n)
              DO i=1,many
                 DO j=1,many
                    xa(j)=yf(jv+jd*(j-1))
                    ya(j)=v(iv+id*(i-1),jv+jd*(j-1))
                 END DO
                 xb(i)=xc(iv+id*(i-1))
                 CALL interpolate(xa,ya,many,yy(n),yb(i),foo,myid,myvertex)
              END DO
              CALL interpolate(xb,yb,many,xx(n),vv(n),foo,myid,myvertex)
           END DO

           ! dudn jump conditions on positive side
           dudnjcp=(4.0d0*uu(1)-uu(2)-3.0d0*uu(0))/ds/2.0d0
           dvdnjcp=(4.0d0*vv(1)-vv(2)-3.0d0*vv(0))/ds/2.0d0

           ddudnp=(-5.0d0*uu(1)+4.0d0*uu(2)-1.0d0*uu(3)+2.0d0*uu(0))/ds/ds
           ddvdnp=(-5.0d0*vv(1)+4.0d0*vv(2)-1.0d0*vv(3)+2.0d0*vv(0))/ds/ds

           ! dudn jump conditions
           dudnjc = dudnjcp - dudnjcn
           dvdnjc = dvdnjcp - dvdnjcn

           ! laplacian u jump conditions
           ! ddudn jump conditions
           deltaujc = ddudnp + vertex(3,myvertex)*dudnjc
           deltavjc = ddvdnp + vertex(3,myvertex)*dvdnjc

           ! body force qx, qy
           qx =  thetatt(myobj)*(xx(0)-ysc(myobj))
           qy = -thetatt(myobj)*(yy(0)-xsc(myobj))

           ! dpdx, dpdy jump conditions
           dpxjcM(myvertex) = 1.0d0/Re*deltaujc + qx
           dpyjcM(myvertex) = 1.0d0/re*deltavjc + qy
        ENDIF
     END DO
     dpxjcM(nvertex4proc(myobj)) = dpxjcM(nvertex4proc(myobj-1)+1)
     dpyjcM(nvertex4proc(myobj)) = dpyjcM(nvertex4proc(myobj-1)+1)

     !  2. Assemble pressure matrix
     DO myvertex = nvertex4proc(myobj-1)+1, nvertex4proc(myobj)-1
        A(1)=xs(myvertex)
        A(2)=ys(myvertex)
        corner0=allcorner0(:,myid)
        corner1=allcorner1(:,myid)
        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN

           ! find left and right index of vertex
           myvertex_left = myvertex + 1
           myvertex_right = myvertex - 1
           IF(myvertex == nvertex4proc(myobj)) THEN
              myvertex_left = nvertex4proc(myobj-1)+2
              myvertex_right = myvertex-1
           ELSEIF(myvertex == nvertex4proc(myobj-1)+1) THEN
              myvertex_left = myvertex+1
              myvertex_right = nvertex4proc(myobj)-1
           END IF

           ! tao direction
           tx_left  =  taox(myvertex)
           ty_left  =  taoy(myvertex)
           tx_right = -taox(myvertex_right)
           ty_right = -taoy(myvertex_right)

           ! a. compute distance of two panels
           ds1 = SQRT((xs(myvertex_left)-xs(myvertex))* &
                (xs(myvertex_left)-xs(myvertex))  &
                +(ys(myvertex_left)-ys(myvertex))* &
                (ys(myvertex_left)-ys(myvertex)))

           ds2 = SQRT((xs(myvertex_right)-xs(myvertex))* &
                (xs(myvertex_right)-xs(myvertex))  &
                +(ys(myvertex_right)-ys(myvertex))* &
                (ys(myvertex_right)-ys(myvertex)))

           ! left
           temp1=(pjc(1,myvertex)+4.0d0*dpxjcM(myvertex)+pjc(1,myvertex+1))*tx_left+ &
                (pjc(2,myvertex)+4.0d0*dpyjcM(myvertex)+pjc(2,myvertex+1))*ty_left
           temp1=temp1*ds1/6.0d0

           ! right
           IF (myvertex == nvertex4proc(myobj-1)+1) THEN
              myvertex_temp = nvertex4proc(myobj)-1
              temp2=(pjc(1,myvertex)+4.0d0*dpxjcM(myvertex_temp)+pjc(1,myvertex_temp))*tx_right+&
                   (pjc(2,myvertex)+4.0d0*dpyjcM(myvertex_temp)+pjc(2,myvertex_temp))*ty_right
           ELSE
              temp2=(pjc(1,myvertex)+4.0d0*dpxjcM(myvertex-1)+pjc(1,myvertex-1))*tx_right+&
                   (pjc(2,myvertex)+4.0d0*dpyjcM(myvertex-1)+pjc(2,myvertex-1))*ty_right
           ENDIF
           temp2=temp2*ds2/6.0d0

           ! update right hand side for conjugate gradient
           jcp_rhs(myvertex)=temp1+temp2
        END IF
     END DO
     jcp_rhs(nvertex4proc(myobj)) = jcp_rhs(nvertex4proc(myobj-1)+1)
  END DO
  DEALLOCATE(dpxjcM)
  DEALLOCATE(dpyjcM)

  time0 = MPI_WTIME()
  ALLOCATE(x_sol(1:nvertex4proc(nobj4proc)))
  x_sol = zero
  rank_c = 0
  ! find ranks of each group
  CALL MPI_COMM_GROUP(comm2d,group0,ierr)
  DO myobj = 1,nobj
     ! if object at local processor, then = 1, else 0
     DO myobj_local = 1, nobj4proc
        IF ( obj4proc(myobj_local) == object_list(myobj) ) THEN
           rank_c(myobj) = 1
           EXIT  ! use exit, otherwise it may be updated to zero
        ELSE
           rank_c(myobj) = 0
        ENDIF
     ENDDO
     rank_temp = 0
     CALL MPI_ALLGATHER(rank_c(myobj),1,MPI_INTEGER,rank_temp,1,MPI_INTEGER,comm2d,ierr)
     ! computer for each obj, how many procs have it
     rank_total(myobj) = SUM(rank_temp)
     ALLOCATE(group_mpi(myobj)%p_rank(rank_total(myobj)))
     ! know for current global object, know which processor has it, get rank
     rank_count = 0
     DO i = 0, nprocs-1
        IF(rank_temp(i) == 1) THEN
           rank_count = rank_count +1
           group_mpi(myobj)%p_rank(rank_count) = i
        ENDIF
     ENDDO
     ! create group for current global object
     CALL MPI_GROUP_INCL(group0,rank_total(myobj),group_mpi(myobj)%p_rank,group_mpi(myobj)%group,ierr)
     CALL MPI_COMM_CREATE(comm2d,group_mpi(myobj)%group,group_mpi(myobj)%comm,ierr)
  ENDDO

  ! gather r_new
  DO myobj = 1, nobj
     ! loop over local object
     IF(group_mpi(myobj)%comm /= MPI_COMM_NULL) THEN
        DO myobj_local = 1, nobj4proc
           IF ( obj4proc(myobj_local) == object_list(myobj) ) THEN
              ! If global and local object matches, then create temporary r vector to store jcp_rhs
              ALLOCATE(group_mpi(myobj)%r_temp(nvertex4proc(myobj_local)-nvertex4proc(myobj_local-1)))
              ! update x_sol & r_new
              DO myvertex = nvertex4proc(myobj_local-1)+1, nvertex4proc(myobj_local)
                 A = vertex(:,myvertex)
                 IF( (A(1)>xc(ipbeg)).AND.(A(2)>yc(jpbeg)).AND. &
                      (A(1)<xc(ipend+1)).AND.(A(2)<yc(jpend+1))) THEN
                    CALL MPI_ALLreduce(jcp_rhs(myvertex), group_mpi(myobj)%r_temp(myvertex-nvertex4proc(myobj_local-1)), &
                         1, MPI_DOUBLE_PRECISION,MPI_SUM,group_mpi(myobj)%comm,ierr)
                 ELSE
                    CALL MPI_ALLreduce(0.0d0, group_mpi(myobj)%r_temp(myvertex-nvertex4proc(myobj_local-1)), &
                         1, MPI_DOUBLE_PRECISION,MPI_SUM,group_mpi(myobj)%comm,ierr)
                 ENDIF
              ENDDO ! myvertex
           END IF ! object match test
        END DO ! myobj_local
     ENDIF
  END DO ! myobj

  DO myobj = 1, nobj
     IF(group_mpi(myobj)%comm /= MPI_COMM_NULL) THEN
        ! loop over local object
        DO myobj_local = 1, nobj4proc
           IF ( obj4proc(myobj_local) == object_list(myobj) ) THEN
              DO myvertex = nvertex4proc(myobj_local-1)+1, nvertex4proc(myobj_local)
                ! get jcp_rhs from allreduced results
                 jcp_rhs(myvertex) = group_mpi(myobj)%r_temp(myvertex-nvertex4proc(myobj_local-1))
              ENDDO ! myvertex
              DEALLOCATE(group_mpi(myobj)%r_temp)
           ENDIF
        END DO ! myobj_local
     END IF !endif(norm)
  END DO ! myobj

  time1 = MPI_WTIME()
  pTime = pTime + time1-time0

  ! solve Topliz matrix
  ! set up matrix A
  ALLOCATE(aa(1:nvertex4proc(nobj4proc)))
  ALLOCATE(bb(1:nvertex4proc(nobj4proc)))
  ALLOCATE(cc(1:nvertex4proc(nobj4proc)))

  aa = 1.0d0
  bb = -2.0d0
  cc = 1.0d0

  DO myobj = 1, nobj4proc

     ! set last x is 0
     x_sol(nvertex4proc(myobj)-1)=0.0d0

     ! change matrix to upper diagonal
     cc(nvertex4proc(myobj-1)+1) = cc(nvertex4proc(myobj-1)+1)/bb(nvertex4proc(myobj-1)+1)
     jcp_rhs(nvertex4proc(myobj-1)+1) = jcp_rhs(nvertex4proc(myobj-1)+1)/bb(nvertex4proc(myobj-1)+1)

     DO myvertex = nvertex4proc(myobj-1)+2, nvertex4proc(myobj)-3
        temp = bb(myvertex) - aa(myvertex)*cc(myvertex-1)
        cc(myvertex) = cc(myvertex)/temp
        jcp_rhs(myvertex) = (jcp_rhs(myvertex)-aa(myvertex)*jcp_rhs(myvertex-1))/temp
     ENDDO
     jcp_rhs(nvertex4proc(myobj)-2) = (jcp_rhs(nvertex4proc(myobj)-2) - &
          aa(nvertex4proc(myobj)-2)*jcp_rhs(nvertex4proc(myobj)-3)) / &
          (bb(nvertex4proc(myobj)-2)-aa(nvertex4proc(myobj)-2)*cc(nvertex4proc(myobj)-3))

     ! back substitute
     x_sol(nvertex4proc(myobj)-2)=jcp_rhs(nvertex4proc(myobj)-2)
     DO myvertex = nvertex4proc(myobj)-3, nvertex4proc(myobj-1)+1, -1
        x_sol(myvertex) = jcp_rhs(myvertex) - cc(myvertex)*x_sol(myvertex+1)
     ENDDO

  ENDDO

  DO myobj = 1, nobj4proc
     DO myvertex = nvertex4proc(myobj-1)+1, nvertex4proc(myobj)-1
        pjc(5, myvertex) = x_sol(myvertex)
        pjc(6, myvertex) = x_sol(myvertex)
     ENDDO
     pjc(5, nvertex4proc(myobj)) = pjc(5, nvertex4proc(myobj-1)+1)
     pjc(6, nvertex4proc(myobj)) = pjc(6, nvertex4proc(myobj-1)+1)
  ENDDO

  DO myobj = 1,nobj
     IF(group_mpi(myobj)%comm /= MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(group_mpi(myobj)%comm, ierr)
     ENDIF
     CALL MPI_GROUP_FREE(group_mpi(myobj)%group, ierr)
  ENDDO

  DEALLOCATE(jcp_rhs)
  DEALLOCATE(x_sol)
  DEALLOCATE(aa)
  DEALLOCATE(bb)
  DEALLOCATE(cc)

END SUBROUTINE principle_jc_p
