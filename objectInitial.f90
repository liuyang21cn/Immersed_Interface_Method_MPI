SUBROUTINE objectInitial

  USE para
  USE Lagrange

  INTEGER:: myobj,myobj_global
  DOUBLE PRECISION,DIMENSION(1:9)::objMove

  ! allocate object properties
  IF(isingular.EQV..TRUE. .AND. myid==0) THEN
     ! if moving, allocate drag lift
     ALLOCATE(cd(1:nobj))
     ALLOCATE(cl(1:nobj))
  END IF
  JCsaved=.FALSE.

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
  ALLOCATE(xs(1:nvertex4proc(nobj4proc)))
  ALLOCATE(ys(1:nvertex4proc(nobj4proc)))
  ALLOCATE(taox(1:nvertex4proc(nobj4proc)))
  ALLOCATE(taoy(1:nvertex4proc(nobj4proc)))
  objAllocated=.TRUE.

  t=ZERO
  t0=ZERO
  tmove=0.0d0

  ! pull local move pattern from gloabl ObjPosition
  ! loop over local object
  DO myobj = 1, nobj4proc
     ! loop over global object
     DO myobj_global = 1, nobj
        IF ( obj4proc(myobj) == object_list(myobj_global) ) THEN
           xsc(myobj)     = ObjPosition(1,myobj_global)
           ysc(myobj)     = ObjPosition(2,myobj_global)
           theta(myobj)   = ObjPosition(3,myobj_global)
           IF(object_move.EQV..TRUE.) THEN
              xsct(myobj)    = ObjPosition(4,myobj_global)
              ysct(myobj)    = ObjPosition(5,myobj_global)
              thetat(myobj)  = ObjPosition(6,myobj_global)
              xsctt(myobj)   = ObjPosition(7,myobj_global)
              ysctt(myobj)   = ObjPosition(8,myobj_global)
              thetatt(myobj) = ObjPosition(9,myobj_global)
           ENDIF
        END IF
     END DO
     ! initial poistion and velocity of vertices
     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)
        ! if vertex is at pre-decided position, then just update xs,ys
        xs(myvertex) = xsc(myobj)+vertex(1,myvertex)*COS(theta(myobj))-vertex(2,myvertex)*SIN(theta(myobj))
        ys(myvertex) = ysc(myobj)+vertex(1,myvertex)*SIN(theta(myobj))+vertex(2,myvertex)*COS(theta(myobj))
        us(myvertex) = ZERO
        vs(myvertex) = ZERO
     END DO
  END DO
  DEALLOCATE(ObjPosition)

END SUBROUTINE objectInitial
