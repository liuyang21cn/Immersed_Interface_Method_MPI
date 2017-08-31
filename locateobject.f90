SUBROUTINE locateobject

  USE MPI
  USE para
  USE Euler
  USE Lagrange
  INTEGER:: status(MPI_STATUS_SIZE)
  INTEGER:: countproc,countobj
  INTEGER:: myproc,myobj,myvertex
  INTEGER:: nvertex,ndimension
  INTEGER:: nlocobj4myproc

  INTEGER, DIMENSION(:), ALLOCATABLE:: object
  INTEGER, DIMENSION(:,:), ALLOCATABLE:: processor4object
  INTEGER, DIMENSION(:,:), ALLOCATABLE:: object4processor

  INTEGER, DIMENSION(:), ALLOCATABLE:: nobjs4procs
  INTEGER, DIMENSION(:), ALLOCATABLE:: objsend2proc

  DOUBLE PRECISION::xx,yy,xsc0,ysc0,theta0

  DOUBLE PRECISION, DIMENSION(1:2):: A,B,corner0,corner1
  DOUBLE PRECISION, DIMENSION(1:3):: vinfo
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: vertices

  CHARACTER(5) filesuffix
  NAMELIST /objectinfo/nobj,object_move
  NAMELIST /vertexinfo/nvertex,ndimension

  INTERFACE
     FUNCTION panelinPE(A,B,corner0,corner1)
       DOUBLE PRECISION, DIMENSION(1:2), INTENT(in):: A,B,corner0,corner1
       LOGICAL panelinPE
     END FUNCTION panelinPE
  END INTERFACE

  ! load object list
  OPEN(100,file='INPUT/OBJECT/objlist')
  READ(100,objectinfo)
  ALLOCATE(object(nobj))
  READ(100,*)(object(myobj),myobj=1,nobj)
  CLOSE(100)

  ! store the object name
  ALLOCATE(object_list(nobj))
  object_list = object

  ! load object moving pattern
  IF(object_move.EQV..FALSE.) THEN
     ALLOCATE(ObjPosition(1:3,1:nobj))
     OPEN(300,file='INPUT/OBJECT/ObjPosition')
     READ(300,*)((ObjPosition(i,myobj),i=1,3),myobj=1,nobj)
     CLOSE(300)
  ELSE
     tmove=zero
     ALLOCATE(ObjPosition(1:9,1:nobj))
     DO myobj=1,nobj
        CALL objUpdate(myobj,0,ObjPosition(:,myobj))
     ENDDO
  ENDIF

  ! object na
  ALLOCATE(object4processor(1:nobj,0:nprocs-1))
  object4processor=0
  IF(myid==0) THEN
     ALLOCATE(processor4object(0:nprocs-1,1:nobj))
     processor4object=-1
     DO myobj=1,nobj
        xsc0  =ObjPosition(1,myobj)
        ysc0  =ObjPosition(2,myobj)
        theta0=ObjPosition(3,myobj)
        WRITE(unit=filesuffix,fmt='(i0)') object(myobj)
        OPEN(200,file='INPUT/OBJECT/vobj'//filesuffix)
        READ(200,vertexinfo)
        ALLOCATE(vertices(3,0:nvertex))
        DO myvertex=0,nvertex
           READ(200,*)(vinfo(i),i=1,3)
           vertices(1,myvertex)=xsc0+vinfo(1)*COS(theta0)-vinfo(2)*SIN(theta0)
           vertices(2,myvertex)=ysc0+vinfo(1)*SIN(theta0)+vinfo(2)*COS(theta0)
           vertices(3,myvertex)=vinfo(3)
        ENDDO
        CLOSE(200)

        countproc=0 ! count number of processors for object
        DO myproc=0,nprocs-1
           DO myvertex=0,nvertex-1
              A=vertices(1:2,myvertex)
              B=vertices(1:2,myvertex+1)
              corner0=allcorner0(:,myproc)
              corner1=allcorner1(:,myproc)
              IF(panelinPE(A,B,corner0,corner1)) THEN
                 processor4object(countproc,myobj)=myproc
                 countproc=countproc+1
                 EXIT ! regard object on myproc if one of its panels on myproc
              ENDIF
           END DO
        END DO
        DEALLOCATE(vertices)
     END DO ! end do myobj=1,nobj

     DO myobj=1,nobj
        DO countproc=0,nprocs-1
           IF(processor4object(countproc,myobj)/=-1) THEN
              myproc=processor4object(countproc,myobj)
              DO countobj=1,nobj
                 IF(object4processor(countobj,myproc)==0) THEN
                    object4processor(countobj,myproc)=object(myobj)
                    EXIT ! push object into stack
                 END IF
              END DO
           END IF
        END DO
     END DO
     DEALLOCATE(object)
     DEALLOCATE(processor4object)
  END IF

  ALLOCATE(nobjs4procs(0:nprocs-1))
  IF(myid==0) THEN
     DO myproc=0,nprocs-1
        nobjs4procs(myproc)=0
        DO myobj=1,nobj
           IF(object4processor(myobj,myproc)/=0) THEN
              nobjs4procs(myproc)=nobjs4procs(myproc)+1
           END IF
        END DO
     END DO
  END IF
  CALL MPI_SCATTER(nobjs4procs,1,MPI_INTEGER,nobj4proc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  IF(myid==0) THEN
     DO myproc=1,nprocs-1
        nlocobj4myproc=nobjs4procs(myproc)
        ALLOCATE(objsend2proc(nlocobj4myproc))
        objsend2proc=object4processor(1:nlocobj4myproc,myproc)
        CALL MPI_SEND(objsend2proc,nlocobj4myproc,MPI_INTEGER,myproc,0,MPI_COMM_WORLD,ierr)
        DEALLOCATE(objsend2proc)
     END DO
     ALLOCATE(obj4proc(nobj4proc))
     obj4proc=object4processor(1:nobj4proc,0)
     DEALLOCATE(object4processor)
  ELSE
     ALLOCATE(obj4proc(nobj4proc))
     CALL MPI_RECV(obj4proc,nobj4proc,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
  END IF
  DEALLOCATE(nobjs4procs)


END SUBROUTINE locateobject
