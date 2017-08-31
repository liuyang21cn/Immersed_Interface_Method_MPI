SUBROUTINE startMPI

  USE MPI
  USE para
  INTEGER:: nxloc,nyloc,rank
  INTEGER, DIMENSION(1:2):: nearcoords

  INTERFACE
     SUBROUTINE indexrange(n,np,id,ibeg,iend)
       INTEGER, INTENT(in):: n,np,id
       INTEGER, INTENT(out):: ibeg,iend
     END SUBROUTINE indexrange

     SUBROUTINE coords2rank(mycoords,myrank)
       INTEGER, DIMENSION(1:2), INTENT(in):: mycoords
       INTEGER, INTENT(out):: myrank
     END SUBROUTINE coords2rank
  END INTERFACE

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

  IF (myid==0) THEN
     WRITE(*,'(A,i5,A)') 'Starting MPI with ', nprocs,' processors.'
  ENDIF

  IF (nprocs/=px*py) THEN
     IF (myid==0) THEN
        PRINT *, ' Error: Incorrect processor layout'
        PRINT *, ' nprocs = ',nprocs
        PRINT *, ' px = ',px,', py = ',py
     ENDIF
     CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
  ENDIF

  dims = (/px, py/)
  periods = .FALSE.
  IF(ipbcx0==0) THEN
     periods(1)=.TRUE.
  END IF
  IF(jpbcy0==0) THEN
     periods(2)=.TRUE.
  END IF

  CALL MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,.TRUE.,comm2d,ierr)
  CALL MPI_CART_COORDS(comm2d,myid,2,coords,ierr)
  CALL MPI_CART_SHIFT(comm2d,0,1,left,right,ierr)
  CALL MPI_CART_SHIFT(comm2d,1,1,front,back,ierr)

  ! 8 near processors surrounding myid, 4 pairs of opposite processors for send and recv
  nearproc(0)=myid
  nearproc(1)=left
  nearproc(2)=right
  nearproc(3)=front
  nearproc(4)=back

  i=coords(1)
  j=coords(2)

  ! ij-plane with 4 corners
  ! pair 1
  nearcoords=(/i-1,j-1/)
  CALL coords2rank(nearcoords,nearproc(5))
  nearcoords=(/i+1,j+1/)
  CALL coords2rank(nearcoords,nearproc(6))
  ! pair 2
  nearcoords=(/i-1,j+1/)
  CALL coords2rank(nearcoords,nearproc(7))
  nearcoords=(/i+1,j-1/)
  CALL coords2rank(nearcoords,nearproc(8))

  ! index ranges of the local subdomain
  ! p
  CALL indexrange(nx,px,coords(1),ipbeg,ipend)
  CALL indexrange(ny,py,coords(2),jpbeg,jpend)

  ! u
  iubeg=ipbeg-1
  IF(ipend.EQ.nx) THEN
     iuend=ipend
  ELSE
     iuend=ipend+1
  ENDIF
  jubeg=jpbeg-1
  juend=jpend+1

  ! v
  ivbeg=ipbeg-1
  ivend=ipend+1
  jvbeg=jpbeg-1
  IF(jpend.EQ.ny) THEN
     jvend=jpend
  ELSE
     jvend=jpend+1
  ENDIF

  ! data types for data exchange
  !p
  nxloc = (ipend+nxghost)-(ipbeg-nxghost)+1
  nyloc = (jpend+nyghost)-(jpbeg-nyghost)+1

  CALL MPI_TYPE_VECTOR(nyloc,nxghost,nxloc,MPI_DOUBLE_PRECISION,pitype,ierr)
  CALL MPI_TYPE_COMMIT(pitype,ierr)
  CALL MPI_TYPE_VECTOR(1,nyghost*nxloc,nxloc*nyloc,MPI_DOUBLE_PRECISION,pjtype,ierr)
  CALL MPI_TYPE_COMMIT(pjtype,ierr)

  !u
  nxloc = (iuend+nxghost)-(iubeg-nxghost)+1
  nyloc = (juend+nyghost)-(jubeg-nyghost)+1

  CALL MPI_TYPE_VECTOR(nyloc,nxghost,nxloc,MPI_DOUBLE_PRECISION,uitype,ierr)
  CALL MPI_TYPE_COMMIT(uitype,ierr)
  CALL MPI_TYPE_VECTOR(1,nyghost*nxloc,nxloc*nyloc,MPI_DOUBLE_PRECISION,ujtype,ierr)
  CALL MPI_TYPE_COMMIT(ujtype,ierr)

  !v
  nxloc = (ivend+nxghost)-(ivbeg-nxghost)+1
  nyloc = (jvend+nyghost)-(jvbeg-nyghost)+1

  CALL MPI_TYPE_VECTOR(nyloc,nxghost,nxloc,MPI_DOUBLE_PRECISION,vitype,ierr)
  CALL MPI_TYPE_COMMIT(vitype,ierr)
  CALL MPI_TYPE_VECTOR(1,nyghost*nxloc,nxloc*nyloc,MPI_DOUBLE_PRECISION,vjtype,ierr)
  CALL MPI_TYPE_COMMIT(vjtype,ierr)

END SUBROUTINE startMPI
