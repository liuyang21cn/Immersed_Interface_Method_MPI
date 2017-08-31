SUBROUTINE loadobject
  USE para
  USE Euler
  USE Lagrange
  INTEGER:: nlocvertex
  INTEGER:: countvertex
  INTEGER:: myobj,myvertex,myproc,neighbor,myobj_global
  INTEGER:: object,nvertex,ndimension

  DOUBLE PRECISION, DIMENSION(1:2):: A,B,corner0,corner1
  DOUBLE PRECISION, DIMENSION(1:3):: vinfo
  DOUBLE PRECISION:: xx,yy,xsc0,ysc0,theta0

  CHARACTER(5) filesuffix
  NAMELIST /vertexinfo/nvertex,ndimension

  INTERFACE
     FUNCTION panelinPE(A,B,corner0,corner1)
       DOUBLE PRECISION, DIMENSION(1:2), INTENT(in):: A,B,corner0,corner1
       LOGICAL panelinPE
     END FUNCTION panelinPE
  END INTERFACE

  ! load vertices
  ALLOCATE(nvertex4proc(0:nobj4proc))
  nvertex4proc=0
  nlocvertex=0
  DO myobj=1,nobj4proc
     object=obj4proc(myobj)
     WRITE(filesuffix,'(i0)') object
     OPEN(100,file='INPUT/OBJECT/vobj'//filesuffix)
     READ(100,vertexinfo)
     nlocvertex=nlocvertex+(nvertex+1)
     nvertex4proc(myobj)=nlocvertex
     CLOSE(100)
  END DO
  ALLOCATE(vertex(3,nlocvertex))

  DO myobj=1,nobj4proc
     object=obj4proc(myobj)
     WRITE(filesuffix,'(i0)') object
     OPEN(100,file='INPUT/OBJECT/vobj'//filesuffix)
     READ(100,vertexinfo)
     countvertex=nvertex4proc(myobj-1)
     DO myvertex=1,nvertex+1
        READ(100,*)(vertex(i,countvertex+myvertex),i=1,3)
     ENDDO
     CLOSE(100)
  END DO

  ALLOCATE(proc4obj(0:8,nobj4proc))
  proc4obj=0

  DO myobj=1,nobj4proc
     DO myobj_global = 1, nobj
        IF ( obj4proc(myobj) == object_list(myobj_global) ) THEN
           xsc0  =ObjPosition(1,myobj_global)
           ysc0  =ObjPosition(2,myobj_global)
           theta0=ObjPosition(3,myobj_global)
        ENDIF
     ENDDO
     DO neighbor=0,8
        myproc=nearproc(neighbor)
        DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
           A(1)=xsc0+vertex(1,myvertex)*COS(theta0)-vertex(2,myvertex)*SIN(theta0)
           A(2)=ysc0+vertex(1,myvertex)*SIN(theta0)+vertex(2,myvertex)*COS(theta0)
           B(1)=xsc0+vertex(1,myvertex+1)*COS(theta0)-vertex(2,myvertex+1)*SIN(theta0)
           B(2)=ysc0+vertex(1,myvertex+1)*SIN(theta0)+vertex(2,myvertex+1)*COS(theta0)
           corner0=allcorner0(:,myproc)
           corner1=allcorner1(:,myproc)
           IF(panelinPE(A,B,corner0,corner1)) THEN
              proc4obj(neighbor,myobj)=1
              EXIT ! regard the object on myproc if one of its panels on myproc
           ENDIF
        END DO
     END DO
  END DO

END SUBROUTINE loadobject
