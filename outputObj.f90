SUBROUTINE outputObj

  USE para
  USE Lagrange

  CHARACTER(5) filesuffix
  INTEGER::myobj,nvertex,ndimension
  DOUBLE PRECISION::xsc0,ysc0,theta0,xs0,ys0
  DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE::vv,vertices

  NAMELIST /vertexinfo/nvertex,ndimension

  IF(myid==0) THEN
     IF(object_move.EQV..FALSE.) THEN
        ALLOCATE(ObjPosition(1:3,1:nobj))
        OPEN(300,file='INPUT/OBJECT/ObjPosition')
        READ(300,*)((ObjPosition(i,myobj),i=1,3),myobj=1,nobj)
        CLOSE(300)
     ELSE
        ALLOCATE(ObjPosition(1:9,1:nobj))
        DO myobj=1,nobj
           CALL objUpdate(myobj,0,ObjPosition(:,myobj))
        ENDDO
     ENDIF

     DO myobj=1,nobj
        ! udpate to current xs,ys
        xsc0  =ObjPosition(1,myobj)
        ysc0  =ObjPosition(2,myobj)
        theta0=ObjPosition(3,myobj)
        ! load original object vertices
        WRITE(unit=filesuffix,fmt='(i0)') object_list(myobj)
        OPEN(200,file='INPUT/OBJECT/vobj'//filesuffix)
        READ(200,vertexinfo)
        ALLOCATE(vertices(3,0:nvertex))
        ALLOCATE(vv(2,0:nvertex))
        READ(200,*)((vertices(i,myvertex),i=1,3),myvertex=0,nvertex)
        CLOSE(200)

        WRITE(unit=filesuffix,fmt='(i0)') myobj ! (i0): using minimum field width appropriate
        OPEN(unit=12,file='OUTPUT/object'//filesuffix,status='unknown')
        DO myvertex=0,nvertex
           vv(1,myvertex)=xsc0+vertices(1,myvertex)*COS(theta0)-vertices(2,myvertex)*SIN(theta0)
           vv(2,myvertex)=ysc0+vertices(1,myvertex)*SIN(theta0)+vertices(2,myvertex)*COS(theta0)
           WRITE(12,100)vv(1,myvertex),vv(2,myvertex)
        ENDDO
        CLOSE(12)
        DEALLOCATE(vertices)
        DEALLOCATE(vv)
     ENDDO
     DEALLOCATE(ObjPosition)
  ENDIF
100 FORMAT(1x,20e16.6)


END SUBROUTINE outputObj
