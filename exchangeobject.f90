SUBROUTINE exchangeobject

  USE para
  USE Euler
  USE Lagrange
  USE MPI

  INTEGER:: status(MPI_STATUS_SIZE)
  INTEGER:: countobj,countvertex,countsend,countrecv
  INTEGER:: ntotalobj4send,ntotalobj4recv,ntotalvertex4recv
  INTEGER:: oldnobj4proc,nobj4delete,nobj4repeat,nlocvertex
  INTEGER:: myobj,myvertex,myproc,neighbor,mylocalvertex,myobj_global

  DOUBLE PRECISION:: xsc0,ysc0,theta0
  INTEGER, DIMENSION(0:8):: nobj4send,nobj4recv
  INTEGER, DIMENSION(:), ALLOCATABLE:: locindex4sendobj,nvertex4recv
  INTEGER, DIMENSION(:), ALLOCATABLE:: oldobj4proc,oldnvertex4proc
  INTEGER, DIMENSION(:,:), ALLOCATABLE:: oldproc4obj,obj4send,obj4recv
  DOUBLE PRECISION, DIMENSION(1:9):: objMove

  DOUBLE PRECISION, DIMENSION(1:2):: A,B,corner0,corner1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: oldvertex,newvertex
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: oldxs,oldys,newxs,newys
  DOUBLE PRECISION:: time0,time1

  INTERFACE
     FUNCTION panelinPE(A,B,corner0,corner1)
       DOUBLE PRECISION, DIMENSION(1:2), INTENT(in):: A,B,corner0,corner1
       LOGICAL panelinPE
     END FUNCTION panelinPE
  END INTERFACE

  time0 = MPI_WTIME()

  ! determine neighboring processors for each object
  ALLOCATE(oldproc4obj(0:8,nobj4proc))
  oldproc4obj=proc4obj
  proc4obj=0
  DO myobj=1,nobj4proc
     DO neighbor=0,8
        myproc=nearproc(neighbor)
        DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
           A(1)=xs(myvertex)
           A(2)=ys(myvertex)
           B(1)=xs(myvertex+1)
           B(2)=ys(myvertex+1)
           corner0=allcorner0(:,myproc)
           corner1=allcorner1(:,myproc)
           IF(panelinPE(A,B,corner0,corner1)) THEN
              proc4obj(neighbor,myobj)=1
              EXIT ! regard the object on myproc if one of its panels on myproc
           ENDIF
        END DO
     END DO
  END DO
  proc4obj=proc4obj-oldproc4obj ! delte(@ 0): -1=0-1, send (@1-26): 1=1-0, no action: 0=1-1=0-0
  DEALLOCATE(oldproc4obj)

  ! determine number of objects to be sent to each neighbor processor
  nobj4send=0
  DO neighbor=1,8
     nobj4send(neighbor)=0
     DO myobj=1,nobj4proc
        IF(proc4obj(neighbor,myobj)==1) THEN
           nobj4send(neighbor)=nobj4send(neighbor)+1
        END IF
     END DO
  END DO

  ! send and receive number of objects for exchange
  nobj4recv=0
  DO neighbor=2,8,2 !send from one direction and receive from its opposite direction, then alternate
     CALL MPI_SENDRECV(nobj4send(neighbor),1,MPI_INTEGER,nearproc(neighbor),0,&
          nobj4recv(neighbor-1),1,MPI_INTEGER,nearproc(neighbor-1),0,comm2d,status,ierr)
     CALL MPI_SENDRECV(nobj4send(neighbor-1),1,MPI_INTEGER,nearproc(neighbor-1),0,&
          nobj4recv(neighbor),1,MPI_INTEGER,nearproc(neighbor),0,comm2d,status,ierr)
  END DO

  ! determine ending index of each object in object stack
  DO neighbor=1,8
     nobj4send(neighbor)=nobj4send(neighbor-1)+nobj4send(neighbor)
  END DO
  ntotalobj4send=nobj4send(8)
  DO neighbor=1,8
     nobj4recv(neighbor)=nobj4recv(neighbor-1)+nobj4recv(neighbor)
  END DO
  ntotalobj4recv=nobj4recv(8)

  ! pack objects4send in stack
  ALLOCATE(obj4send(2,ntotalobj4send))
  obj4send=0
  ALLOCATE(locindex4sendobj(ntotalobj4send))
  locindex4sendobj=0

  countobj=1
  DO neighbor=1,8
     DO myobj=1,nobj4proc
        IF(proc4obj(neighbor,myobj)==1) THEN
           locindex4sendobj(countobj)=myobj ! local object index
           obj4send(1,countobj)=obj4proc(myobj) ! global object index
           obj4send(2,countobj)=nvertex4proc(myobj)-nvertex4proc(myobj-1) ! number of vertices
           countobj=countobj+1
        END IF
     END DO
  END DO

  ! exchange object info with neighbors
  ALLOCATE(obj4recv(2,ntotalobj4recv))
  obj4recv=0
  DO neighbor=2,8,2
     countsend=nobj4send(neighbor)-nobj4send(neighbor-1)
     IF(countsend/=0) THEN
        CALL MPI_SEND(obj4send(1,nobj4send(neighbor-1)+1),2*countsend,MPI_INTEGER,nearproc(neighbor),0,comm2d,ierr)
     END IF
     countrecv=nobj4recv(neighbor-1)-nobj4recv(neighbor-2)
     IF(countrecv/=0) THEN
        CALL MPI_RECV(obj4recv(1,nobj4recv(neighbor-2)+1),2*countrecv,MPI_INTEGER,nearproc(neighbor-1),0,comm2d,status,ierr)
     END IF
     countsend=nobj4send(neighbor-1)-nobj4send(neighbor-2)
     IF(countsend/=0) THEN
        CALL MPI_SEND(obj4send(1,nobj4send(neighbor-2)+1),2*countsend,MPI_INTEGER,nearproc(neighbor-1),0,comm2d,ierr)
     END IF
     countrecv=nobj4recv(neighbor)-nobj4recv(neighbor-1)
     IF(countrecv/=0) THEN
        CALL MPI_RECV(obj4recv(1,nobj4recv(neighbor-1)+1),2*countrecv,MPI_INTEGER,nearproc(neighbor),0,comm2d,status,ierr)
     END IF
  END DO
  DEALLOCATE(obj4send)

  ! determine number of panels and vertices to be received
  ALLOCATE(nvertex4recv(0:ntotalobj4recv))
  nvertex4recv=0
  DO neighbor=1,8
     DO countobj=nobj4recv(neighbor-1)+1,nobj4recv(neighbor)
        nvertex4recv(countobj)=nvertex4recv(countobj-1)+obj4recv(2,countobj)
     END DO
  END DO
  ntotalvertex4recv=nvertex4recv(ntotalobj4recv)

  ! exchange vertices, panels and curvatures
  ALLOCATE(newvertex(3,ntotalvertex4recv))

  ! vertices
  DO neighbor=2,8,2
     DO countobj=nobj4send(neighbor-1)+1,nobj4send(neighbor)
        myobj=locindex4sendobj(countobj)
        myvertex=nvertex4proc(myobj-1)+1
        countvertex=nvertex4proc(myobj)-nvertex4proc(myobj-1)
        IF(countvertex/=0) THEN
           CALL MPI_SEND(vertex(1,myvertex),3*countvertex,MPI_DOUBLE_PRECISION,nearproc(neighbor),0,comm2d,ierr)
        END IF
     END DO
     DO countobj=nobj4recv(neighbor-2)+1,nobj4recv(neighbor-1)
        myvertex=nvertex4recv(countobj-1)+1
        countvertex=nvertex4recv(countobj)-nvertex4recv(countobj-1)
        IF(countvertex/=0) THEN
           CALL MPI_RECV(newvertex(1,myvertex),3*countvertex,MPI_DOUBLE_PRECISION,nearproc(neighbor-1),0,comm2d,status,ierr)
        END IF
     END DO
     DO countobj=nobj4send(neighbor-2)+1,nobj4send(neighbor-1)
        myobj=locindex4sendobj(countobj)
        myvertex=nvertex4proc(myobj-1)+1
        countvertex=nvertex4proc(myobj)-nvertex4proc(myobj-1)
        IF(countvertex/=0) THEN
           CALL MPI_SEND(vertex(1,myvertex),3*countvertex,MPI_DOUBLE_PRECISION,nearproc(neighbor-1),0,comm2d,ierr)
        END IF
     END DO
     DO countobj=nobj4recv(neighbor-1)+1,nobj4recv(neighbor)
        myvertex=nvertex4recv(countobj-1)+1
        countvertex=nvertex4recv(countobj)-nvertex4recv(countobj-1)
        IF(countvertex/=0) THEN
           CALL MPI_RECV(newvertex(1,myvertex),3*countvertex,MPI_DOUBLE_PRECISION,nearproc(neighbor),0,comm2d,status,ierr)
        END IF
     END DO
  END DO
  DEALLOCATE(locindex4sendobj)

  ! save old object data
  ALLOCATE(oldobj4proc(nobj4proc))
  ALLOCATE(oldnvertex4proc(0:nobj4proc))
  oldobj4proc=obj4proc
  oldnvertex4proc=nvertex4proc
  DEALLOCATE(obj4proc)
  DEALLOCATE(nvertex4proc)
  nlocvertex=oldnvertex4proc(nobj4proc)
  ALLOCATE(oldvertex(3,nlocvertex))
  oldvertex=vertex

  DEALLOCATE(vertex)
  DEALLOCATE(xs)
  DEALLOCATE(ys)

  ! count number of received repeated objects
  nobj4repeat=0
  DO myobj=1,ntotalobj4recv
     DO countobj=1,myobj-1
        IF(obj4recv(1,countobj)/=0.AND.obj4recv(1,myobj)==obj4recv(1,countobj)) THEN
           obj4recv(1,myobj)=0
           nobj4repeat=nobj4repeat+1
           EXIT
        END IF
     END DO
  END DO

  ! determine number of objects out of processor myid (@0)
  nobj4delete=0
  DO myobj=1,nobj4proc
     nobj4delete=nobj4delete-proc4obj(0,myobj)
  END DO

  ! update local object info
  oldnobj4proc=nobj4proc
  nobj4proc=oldnobj4proc-nobj4delete+ntotalobj4recv-nobj4repeat
  ALLOCATE(obj4proc(nobj4proc))
  ALLOCATE(nvertex4proc(0:nobj4proc))
  obj4proc=0
  nvertex4proc=0
  countobj=1
  DO myobj=1,oldnobj4proc
     IF(proc4obj(0,myobj)==-1) THEN
        oldobj4proc(myobj)=0
     ELSE
        obj4proc(countobj)=oldobj4proc(myobj)
        nvertex4proc(countobj)=nvertex4proc(countobj-1)+(oldnvertex4proc(myobj)-oldnvertex4proc(myobj-1))
        countobj=countobj+1
     END IF
  END DO
  DO myobj=1,ntotalobj4recv
     IF(obj4recv(1,myobj)/=0) THEN
        obj4proc(countobj)=obj4recv(1,myobj)
        nvertex4proc(countobj)=nvertex4proc(countobj-1)+obj4recv(2,myobj)
        countobj=countobj+1
     END IF
  END DO
  DEALLOCATE(oldobj4proc)
  nlocvertex=nvertex4proc(nobj4proc)

  ! update local vertices, panels and curvatures
  ALLOCATE(vertex(3,nlocvertex))
  ALLOCATE(xs(nlocvertex))
  ALLOCATE(ys(nlocvertex))
  countobj=1
  DO myobj=1,oldnobj4proc
     IF(proc4obj(0,myobj)/=-1) THEN
        vertex(:,nvertex4proc(countobj-1)+1:nvertex4proc(countobj))=&
             oldvertex(:,oldnvertex4proc(myobj-1)+1:oldnvertex4proc(myobj))
        countobj=countobj+1
     END IF
  END DO

  DEALLOCATE(proc4obj)
  DEALLOCATE(oldnvertex4proc)
  DEALLOCATE(oldvertex)
  DO myobj=1,ntotalobj4recv
     IF(obj4recv(1,myobj)/=0) THEN
        vertex(:,nvertex4proc(countobj-1)+1:nvertex4proc(countobj))=&
             newvertex(:,nvertex4recv(myobj-1)+1:nvertex4recv(myobj))
        countobj=countobj+1
     END IF
  END DO

  DEALLOCATE(obj4recv)
  DEALLOCATE(nvertex4recv)
  DEALLOCATE(newvertex)

  ! determine neighboring processors for each object
  ALLOCATE(proc4obj(0:8,nobj4proc))
  proc4obj=0

  DO myobj=1,nobj4proc
     DO myobj_global = 1, nobj
        IF ( obj4proc(myobj) == object_list(myobj_global) ) THEN
           CALL objUpdate(myobj_global,0,objMove)
           xsc0     = objMove(1)
           ysc0     = objMove(2)
           theta0   = objMove(3)
        ENDIF
     ENDDO
     DO neighbor=0,8
        myproc=nearproc(neighbor)
        DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
           xs(myvertex)=xsc0+vertex(1,myvertex)*COS(theta0)-vertex(2,myvertex)*SIN(theta0)
           ys(myvertex)=ysc0+vertex(1,myvertex)*SIN(theta0)+vertex(2,myvertex)*COS(theta0)
           xs(myvertex+1)=xsc0+vertex(1,myvertex+1)*COS(theta0)-vertex(2,myvertex+1)*SIN(theta0)
           ys(myvertex+1)=ysc0+vertex(1,myvertex+1)*SIN(theta0)+vertex(2,myvertex+1)*COS(theta0)

           A(1)=xs(myvertex)
           A(2)=ys(myvertex)
           B(1)=xs(myvertex+1)
           B(2)=ys(myvertex+1)
           corner0=allcorner0(:,myproc)
           corner1=allcorner1(:,myproc)
           IF(panelinPE(A,B,corner0,corner1)) THEN
              proc4obj(neighbor,myobj)=1
              EXIT ! regard the object on myproc if one of its panels on myproc
           ENDIF
        END DO
     END DO
  END DO

  time1 = MPI_WTIME()
  movingTime = movingTime + time1-time0

END SUBROUTINE exchangeobject
