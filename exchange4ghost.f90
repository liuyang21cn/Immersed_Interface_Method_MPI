SUBROUTINE exchange4ghost

  USE MPI
  USE para
  USE Euler

  INTEGER:: status(MPI_STATUS_SIZE)
  DOUBLE PRECISION:: time0,time1
  time0 = MPI_WTIME()

  ! u
  ! i-direction
  CALL MPI_SENDRECV(u(iuend-nxghost,   jubeg-nyghost),1,uitype,right,0,&
       u(iubeg-nxghost+1, jubeg-nyghost),1,uitype,left, 0,&
       comm2d,status,ierr)
  CALL MPI_SENDRECV(u(iubeg+1,jubeg-nyghost),1,uitype,left, 0,&
       u(iuend,  jubeg-nyghost),1,uitype,right,0,&
       comm2d,status,ierr)
  ! j-direction
  CALL MPI_SENDRECV(u(iubeg-nxghost, juend-nyghost  ),1,ujtype,back, 0,&
       u(iubeg-nxghost, jubeg-nyghost+1),1,ujtype,front,0,&
       comm2d,status,ierr)
  CALL MPI_SENDRECV(u(iubeg-nxghost, jubeg+1),1,ujtype,front,0,&
       u(iubeg-nxghost, juend  ),1,ujtype,back, 0,&
       comm2d,status,ierr)

  ! v
  ! i-direction
  CALL MPI_SENDRECV(v(ivend-nxghost,    jvbeg-nyghost),1,vitype,right,0,&
       v(ivbeg-nxghost+1,  jvbeg-nyghost),1,vitype,left, 0,&
       comm2d,status,ierr)
  CALL MPI_SENDRECV(v(ivbeg+1,  jvbeg-nyghost),1,vitype,left, 0,&
       v(ivend,    jvbeg-nyghost),1,vitype,right,0,&
       comm2d,status,ierr)
  ! j-direction
  CALL MPI_SENDRECV(v(ivbeg-nxghost,jvend-nyghost  ),1,vjtype,back, 0,&
       v(ivbeg-nxghost,jvbeg-nyghost+1),1,vjtype,front,0,&
       comm2d,status,ierr)
  CALL MPI_SENDRECV(v(ivbeg-nxghost,jvbeg+1),1,vjtype,front,0,&
       v(ivbeg-nxghost,jvend  ),1,vjtype,back, 0,&
       comm2d,status,ierr)

  time1 = MPI_WTIME()
  ghostTime = ghostTime + time1-time0

END SUBROUTINE exchange4ghost
