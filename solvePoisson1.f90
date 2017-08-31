SUBROUTINE solvePoisson(flag)

  USE MPI
  USE para
  USE Euler
  USE stretchxy

  ! flag = 1, solve pressure else solve stream function
  INTEGER,INTENT(in):: flag
  INTEGER:: iters,status(MPI_STATUS_SIZE)
  INTEGER, DIMENSION(1:2):: indices
  DOUBLE PRECISION:: resid,pnorm8,pnorm,sum1,sum_p,aver_p,summ,summ1

  ! PCG solver with SMG preconditioner

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

  ! Step 5: setup the solver
  CALL HYPRE_StructPCGCreate(MPI_COMM_WORLD, solver, ierr)
  CALL HYPRE_StructPCGSetMaxIter(solver, maxit, ierr )
  CALL HYPRE_StructPCGSetTol(solver, tol, ierr)
  CALL HYPRE_StructPCGSetTwoNorm(solver, 1, ierr )
  CALL HYPRE_StructPCGSetRelChange(solver, 0, ierr )
  ! CALL HYPRE_StructPCGSetPrintLevel(solver, 2, ierr ) !/* PRINT each CG iteration */
  CALL HYPRE_StructPCGSetLogging(solver, 1, ierr)

  ! 349 /* USE symmetric SMG as preconditioner */
  CALL HYPRE_StructSMGCreate(MPI_COMM_WORLD, precond, ierr)
  CALL HYPRE_StructSMGSetMemoryUse(precond, 0, ierr)
  CALL HYPRE_StructSMGSetMaxIter(precond, 1, ierr)
  CALL HYPRE_StructSMGSetTol(precond, 0.0, ierr)
  CALL HYPRE_StructSMGSetNonZeroGuess(precond, ierr)
  CALL HYPRE_StructSMGSetNumPreRelax(precond, npre, ierr)
  CALL HYPRE_StructSMGSetNumPostRelax(precond, npost, ierr)

  ! /* Set the preconditioner and solve */
  call HYPRE_StructPCGSetPrecond(solver, 0, precond, ierr)
  CALL HYPRE_StructPCGSetup(solver,Lmatrix,bvector,xvector,ierr)
  CALL HYPRE_StructPCGSolve(solver,Lmatrix,bvector,xvector,ierr)

  ! /* Get some info on the run */
  CALL HYPRE_StructPCGGetNumIterations(solver,iters,ierr)
  CALL HYPRE_StructPCGGetFinalRelative(solver,resid,ierr)
  IF(myid==0) THEN
     PRINT*
     IF(flag==1) PRINT*, "Pressure Poisson Iterations = ", iters
     IF(flag==0) PRINT*, "Stream Function Iterations = ", iters
     PRINT*, "Final Relative Residual Norm = ", resid
     PRINT*
  ENDIF

  IF(iters.EQ.maxit) THEN
     WRITE(*,*)'Hyper reaches to maximum iterations.'
     STOP
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

  ! exchange pressure
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

  ! clean up pressure
  ! sum1=zero
  ! sum_p=zero
  ! DO i=ipbeg,ipend
  !    DO j=jpbeg,jpend
  !       IF(flag.EQ.1) THEN
  !          sum1=sum1+p(i,j)
  !       ELSE
  !          sum1=sum1+ph(i,j)
  !       ENDIF
  !    ENDDO
  ! ENDDO
  ! CALL MPI_Allreduce(sum1,sum_p,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm2d,ierr)
  ! aver_p=sum_p/nx/ny
  ! IF(myid.EQ.0) THEN
  !    IF(flag.EQ.1) THEN
  !       PRINT*, 'aver_p',aver_p
  !    ELSE
  !       PRINT*,'aver stream',aver_p
  !    ENDIF
  ! ENDIF
  ! IF(flag.EQ.1) THEN
  !    p=p-aver_p
  ! ELSE
  !    ph=ph-aver_p
  ! ENDIF

  CALL HYPRE_StructSMGDestroy(precond,ierr)
  CALL HYPRE_StructPCGDestroy(solver, ierr)


END SUBROUTINE solvePoisson
