SUBROUTINE surface_move(krk,fac)

  USE para
  USE Euler
  USE Lagrange

  INTEGER,INTENT(in):: krk
  DOUBLE PRECISION:: fac
  INTEGER:: myobj,myobj_global
  DOUBLE PRECISION, DIMENSION(1:9):: objMove

  IF(krk.EQ.1) tmove=tmove+0.5d0*dt
  IF(krk.EQ.3) tmove=tmove+0.5d0*dt
  IF(krk.EQ.5) tmove=tmove+dt

  DO myobj=1,nobj4proc
     ! compute new xsc,ysc
     DO myobj_global = 1, nobj
        IF ( obj4proc(myobj) == object_list(myobj_global) ) THEN
           ! calculate new object movement
           CALL objUpdate(myobj_global,krk,objMove)
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

     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)
        xs(myvertex)=xsc(myobj)+vertex(1,myvertex)*COS(theta(myobj))-vertex(2,myvertex)*SIN(theta(myobj))
        ys(myvertex)=ysc(myobj)+vertex(1,myvertex)*SIN(theta(myobj))+vertex(2,myvertex)*COS(theta(myobj))
     ENDDO
  ENDDO

  CALL exchangeobject

  IF(objAllocated.EQV..TRUE.) THEN
     DEALLOCATE(xsc)
     DEALLOCATE(ysc)
     DEALLOCATE(theta)
     DEALLOCATE(xsct)
     DEALLOCATE(ysct)
     DEALLOCATE(thetat)
     DEALLOCATE(xsctt)
     DEALLOCATE(ysctt)
     DEALLOCATE(thetatt)
     DEALLOCATE(us)
     DEALLOCATE(vs)
     DEALLOCATE(taox)
     DEALLOCATE(taoy)
     objAllocated=.FALSE.
  ENDIF

END SUBROUTINE surface_move
