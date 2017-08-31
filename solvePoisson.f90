SUBROUTINE solvePoisson(flag)

  USE MPI
  USE para
  USE Euler
  USE stretchxy

  ! flag = 1, solve pressure; else solve stream function
  INTEGER,INTENT(in):: flag
  INTEGER:: iters,status(MPI_STATUS_SIZE)
  INTEGER, DIMENSION(1:2):: indices
  DOUBLE PRECISION:: resid,pnorm8,pnorm,sum1,sum_p,aver_p,summ,summ1
  DOUBLE PRECISION:: time0,time1,time2,time3

  IF(flag==0) tol = 1.0d-13
  IF(myid==0)  PRINT*,'current tol = ',tol
  ! assign values to bvector and xvector
  DO j = jpbeg,jpend
     indices(2) = j
     DO i = ipbeg,ipend
        indices(1) = i
        IF(flag.EQ.1) THEN
           CALL HYPRE_StructVectorSetValues(bvector,indices,rhsp(i,j),ierr)
           CALL HYPRE_StructVectorSetValues(xvector,indices,p(i,j),ierr)
        ELSE
           CALL HYPRE_StructVectorSetValues(bvector,indices,rhsph(i,j),ierr)
           CALL HYPRE_StructVectorSetValues(xvector,indices,ph(i,j),ierr)
        ENDIF
     END DO
  END DO
  CALL HYPRE_StructVectorAssemble(bvector,ierr)
  CALL HYPRE_StructVectorAssemble(xvector,ierr)

100 CONTINUE

  time2 = MPI_WTIME()
  ! Step 5: setup the solver
  CALL HYPRE_StructSMGCreate(MPI_COMM_WORLD,solver,ierr)
  CALL HYPRE_StructSMGSetMemoryUse(solver, 0, ierr);
  CALL HYPRE_StructSMGSetMaxIter(solver, maxit,ierr)
  IF (flag==1) THEN
     CALL HYPRE_StructSMGSetTol(solver, tol,ierr)
  ELSE
     CALL HYPRE_StructSMGSetTol(solver,tol,ierr)
  ENDIF
  CALL HYPRE_StructSMGSetRelChange(solver, relch,ierr)
  CALL HYPRE_StructSMGSetZeroGuess(solver, ierr)
  CALL HYPRE_StructSMGSetNumPreRelax(solver, npre,ierr)
  CALL HYPRE_StructSMGSetNumPostRelax(solver, npost,ierr)
  CALL HYPRE_StructSMGSetPrintLevel(solver, 0,ierr)
  CALL HYPRE_StructSMGSetLogging(solver, 1,ierr)


  ! solve the system: Lmatrix*xvector=bvector
  CALL HYPRE_StructSMGSetup(solver,Lmatrix,bvector,xvector,ierr)
  CALL HYPRE_StructSMGSolve(solver,Lmatrix,bvector,xvector,ierr)

  ! return the norm of the final relative residual and the number of iterations
  CALL HYPRE_StructSMGGetFinalRelative(solver,resid,ierr)
  CALL HYPRE_StructSMGGetNumIterations(solver,iters,ierr)

  time3 = MPI_WTIME()
  hypreTime = hypreTime + time3-time2

  IF(myid==0) THEN
     PRINT*
     IF(flag==1) PRINT*, "Pressure Poisson Iterations = ", iters
     IF(flag==0) PRINT*, "Stream Function Iterations = ", iters
     PRINT*, "Final Relative Residual Norm = ", resid
     PRINT*
  ENDIF

  IF(iters.EQ.maxit) THEN
     IF(myid==0) WRITE(*,*)'Hyper reaches to maximum iterations.'
     tol=tol*5
     IF(myid==0) PRINT*,'new tol = ', tol
     IF(tol>1.0d-3) THEN
        IF(myid==0) PRINT*,'tol too big'
        STOP
     ENDIF
     CALL HYPRE_StructSMGDestroy(solver,ierr)
     go to 100
  ENDIF

  ! arrange the solution in the 2d array
  DO j = jpbeg,jpend
     indices(2) = j
     DO i = ipbeg,ipend
        indices(1) = i
        IF(flag.EQ.1) THEN
           CALL HYPRE_StructVectorGetValues(xvector,indices,p(i,j),ierr)
        ELSE
           CALL HYPRE_StructVectorGetValues(xvector,indices,ph(i,j),ierr)
        ENDIF
     END DO
  END DO
  CALL HYPRE_StructSMGDestroy(solver,ierr)

  ! exchange pressure
  time0 = MPI_WTIME()

  IF(flag.EQ.1) THEN
     ! i-direction
     CALL MPI_SENDRECV(p(ipend-nxghost+1, jpbeg-nyghost),1,pitype,right,0,&
          p(ipbeg-nxghost, jpbeg-nyghost),1,pitype,left, 0,&
          comm2d,status,ierr)
     CALL MPI_SENDRECV(p(ipbeg,jpbeg-nyghost),1,pitype,left, 0,&
          p(ipend+1,  jpbeg-nyghost),1,pitype,right,0,&
          comm2d,status,ierr)
     ! j-direction
     CALL MPI_SENDRECV(p(ipbeg-nxghost, jpend-nyghost+1),1,pjtype,back, 0,&
          p(ipbeg-nxghost, jpbeg-nyghost),1,pjtype,front,0,&
          comm2d,status,ierr)
     CALL MPI_SENDRECV(p(ipbeg-nxghost, jpbeg),1,pjtype,front,0,&
          p(ipbeg-nxghost, jpend+1),1,pjtype,back, 0,&
          comm2d,status,ierr)
  ENDIF

  time1 = MPI_WTIME()
  ghostTime = ghostTime + time1-time0


END SUBROUTINE solvePoisson
