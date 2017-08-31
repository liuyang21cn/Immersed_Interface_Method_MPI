SUBROUTINE surface_property

  USE mpi
  USE Lagrange
  USE para
  DOUBLE PRECISION:: gacobi,taoxx,taoyy
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE::coord
  DOUBLE PRECISION, DIMENSION(1:9):: objMove
  INTEGER:: myobj_global

  IF(object_move.EQV..TRUE.) THEN
     IF(objAllocated.EQV..FALSE.) THEN
        ! objects properties
        ALLOCATE(xsc(1:nobj4proc))
        ALLOCATE(ysc(1:nobj4proc))
        ALLOCATE(theta(1:nobj4proc))
        ALLOCATE(xsct(1:nobj4proc))
        ALLOCATE(ysct(1:nobj4proc))
        ALLOCATE(thetat(1:nobj4proc))
        ALLOCATE(xsctt(1:nobj4proc))
        ALLOCATE(ysctt(1:nobj4proc))
        ALLOCATE(thetatt(1:nobj4proc))
        ALLOCATE(us(1:nvertex4proc(nobj4proc)))
        ALLOCATE(vs(1:nvertex4proc(nobj4proc)))
        ALLOCATE(taox(1:nvertex4proc(nobj4proc)))
        ALLOCATE(taoy(1:nvertex4proc(nobj4proc)))
        objAllocated=.TRUE.
     ENDIF
  ENDIF

  DO myobj=1,nobj4proc
     ! update xsc,ysc
     IF(object_move.EQV..TRUE.) THEN
        DO myobj_global = 1, nobj
           IF ( obj4proc(myobj) == object_list(myobj_global) ) THEN
              CALL objUpdate(myobj_global,0,objMove)
              xsc(myobj)     = objMove(1)
              ysc(myobj)     = objMove(2)
              theta(myobj)   = objMove(3)
              xsct(myobj)    = objMove(4)
              ysct(myobj)    = objMove(5)
              thetat(myobj)  = objMove(6)
              xsctt(myobj)   = objMove(7)
              ysctt(myobj)   = objMove(8)
              thetatt(myobj) = objMove(9)
           ENDIF
        ENDDO
        ! compute coordinates and update velocity of vertices
        DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)
           ! update velocity of vertices
           xs(myvertex) = xsc(myobj)+vertex(1,myvertex)*COS(theta(myobj))-vertex(2,myvertex)*SIN(theta(myobj))
           ys(myvertex) = ysc(myobj)+vertex(1,myvertex)*SIN(theta(myobj))+vertex(2,myvertex)*COS(theta(myobj))
           us(myvertex) = xsct(myobj)-thetat(myobj)*(ys(myvertex)-ysc(myobj))
           vs(myvertex) = ysct(myobj)+thetat(myobj)*(xs(myvertex)-xsc(myobj))
        ENDDO
     ENDIF

     ! compute tao_x and tao_y
     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
        taoxx = xs(myvertex+1)-xs(myvertex)
        taoyy = ys(myvertex+1)-ys(myvertex)
        gacobi=SQRT(taoxx*taoxx+taoyy*taoyy)
        taox(myvertex) = taoxx/gacobi
        taoy(myvertex) = taoyy/gacobi
     ENDDO
     taox(nvertex4proc(myobj)) = taox(nvertex4proc(myobj-1)+1)
     taoy(nvertex4proc(myobj)) = taoy(nvertex4proc(myobj-1)+1)
  ENDDO

END SUBROUTINE surface_property
