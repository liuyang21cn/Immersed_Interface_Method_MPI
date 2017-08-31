PROGRAM main

  USE MPI
  USE para
  USE Euler
  USE Lagrange
  DOUBLE PRECISION:: time0, time1,ttotal,t1,t2,t3,t4

  CALL MPI_INIT(ierr)
  time0 = MPI_WTIME()

  CALL readpara
  CALL startMPI
  CALL Cartesiangrid

  IF(nobj>0) THEN
     CALL locateobject
     CALL loadobject
     CALL objectInitial
  END IF
  CALL allocateEuler
  IF(nstep==0) CALL outputObj

  CALL setuphypre
  CALL initial

  ! march in time
  nstart=0
  nstart=nstart+1
  nend=nstart+nstep-1
  tstart=t
  DO nn=1, nstep
     CALL cfl
     IF(myid==0) PRINT *, 'nn = ',nn, ' t =',t,' dt =',dt
     CALL rk
     IF(isingular.EQV..TRUE.) THEN
        IF(MOD(nn,2)==1 .AND. myid==0) THEN
           WRITE(86,100)t,(cd(myobj),cl(myobj),myobj=1,nobj)
        ENDIF
     ENDIF
     IF(nn==nstep) CALL outputObj
  ENDDO

  CALL run_output
  CALL outputEuler
100 FORMAT(1x,20e16.6)

  time1 = MPI_WTIME()
  CALL MPI_Reduce(time1-time0,ttotal,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,comm2d,ierr)
  CALL MPI_Reduce(ptime,     t1,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,comm2d,ierr)
  CALL MPI_Reduce(hypreTime, t2,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,comm2d,ierr)
  CALL MPI_Reduce(ghostTime, t3,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,comm2d,ierr)
  CALL MPI_Reduce(movingTime,t4,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,comm2d,ierr)


  IF(myid.EQ.0) THEN
     PRINT*, 'Computational time =', ttotal, 'seconds'
     PRINT*, 'principle of p     =', t1, 'seconds'
     PRINT*, 'hypre time         =', t2, 'seconds'
     PRINT*, 'ghost layer time   =', t3, 'seconds'
     PRINT*, 'object moving time =', t4, 'seconds'
  ENDIF
  CALL MPI_FINALIZE(ierr)

END PROGRAM main
