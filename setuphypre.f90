SUBROUTINE setuphypre
  USE MPI
  USE para
  USE Euler

  INTEGER, DIMENSION(1:2):: ilower,iupper,offsets,indices
  INTEGER, DIMENSION(0:4):: entries
  DOUBLE PRECISION, DIMENSION(0:4):: values
  DOUBLE PRECISION:: bcvalue1,bcvalue2
  ! DiscreteLaplacian(p) = sp => Lmatrix*p = sp*dxi*dxi := rhs4p
  ! Lmatrix = DiscreteLaplacian*dxi*dxi

  ! set the local extents
  ilower = (/ipbeg,jpbeg/)
  iupper = (/ipend,jpend/)

  ! Step 1: create the 3d grid
  CALL HYPRE_StructGridCreate(MPI_COMM_WORLD,2,grid,ierr)
  CALL HYPRE_StructGridSetExtents(grid,ilower,iupper,ierr)
  CALL HYPRE_StructGridAssemble(grid,ierr)

  ! Step 2: set up the stencil
  CALL HYPRE_StructStencilCreate(2,5,stencil,ierr)
  offsets = (/0,0/)  ! center
  CALL HYPRE_StructStencilSetElement(stencil,0,offsets,ierr)
  offsets = (/-1,0/) ! left
  CALL HYPRE_StructStencilSetElement(stencil,1,offsets,ierr)
  offsets = (/1,0/)  ! right
  CALL HYPRE_StructStencilSetElement(stencil,2,offsets,ierr)
  offsets = (/0,-1/) ! front
  CALL HYPRE_StructStencilSetElement(stencil,3,offsets,ierr)
  offsets = (/0,1/)  ! back
  CALL HYPRE_StructStencilSetElement(stencil,4,offsets,ierr)

  ! Step 3: setup the matrix Lmatrix
  CALL HYPRE_StructMatrixCreate(MPI_COMM_WORLD,grid,stencil,Lmatrix,ierr)
  CALL HYPRE_StructMatrixInitialize(Lmatrix,ierr)
  entries = (/0,1,2,3,4/) ! stencil entries
  ! set the entry of the matrix Lmatrix point by point
  DO j = jpbeg,jpend
     indices(2) = j
     DO i = ipbeg,ipend
        indices(1) = i
        values(1) = xicx(i)*xifx(i-1)
        values(2) = xicx(i)*xifx(i)
        values(3) = r2dxideta*etacy(j)*etafy(j-1)
        values(4) = r2dxideta*etacy(j)*etafy(j)
        values(0) = -(values(1)+values(2)+values(3)+values(4))
        CALL HYPRE_StructMatrixSetValues(Lmatrix,indices,5,entries,values,ierr)
     ENDDO
  ENDDO

  ! set the boundary grid points based on BCs
  ! left
  IF(ipbeg==1) THEN
     i=1
     SELECT CASE(ipbcx0)
     CASE(0)
     CASE(1)
        indices(1) = i
        DO j = jpbeg,jpend
           indices(2) =  j
           bcvalue1 = ZERO
           CALL HYPRE_StructMatrixSetValues(Lmatrix,indices,1,1,bcvalue1,ierr)
        END DO
     CASE(2)
        indices(1) = i
        DO j = jpbeg,jpend
           indices(2) =  j
           bcvalue1 = ZERO
           bcvalue2 = xicx(i)*xifx(i)+xicx(i)*xifx(i-1)
           CALL HYPRE_StructMatrixSetValues(Lmatrix,indices,1,1,bcvalue1,ierr)
           CALL HYPRE_StructMatrixSetValues(Lmatrix,indices,1,2,bcvalue2,ierr)
        END DO
     END SELECT
  END IF
  ! right
  IF(ipend==nx) THEN
     i=nx
     SELECT CASE(ipbcx1)
     CASE(0)
     CASE(1)
        indices(1) = i
        DO j = jpbeg,jpend
           indices(2) =  j
           bcvalue1 = ZERO
           CALL HYPRE_StructMatrixSetValues(Lmatrix,indices,1,2,bcvalue1,ierr)
        END DO
     CASE(2)
        indices(1) = i
        DO j = jpbeg,jpend
           indices(2) =  j
           bcvalue1 = ZERO
           bcvalue2 = xicx(i)*xifx(i-1)+xicx(i)*xifx(i)
           CALL HYPRE_StructMatrixSetValues(Lmatrix,indices,1,2,bcvalue1,ierr)
           CALL HYPRE_StructMatrixSetValues(Lmatrix,indices,1,1,bcvalue2,ierr)
        END DO
     END SELECT
  END IF
  ! front
  IF(jpbeg==1) THEN
     j=1
     SELECT CASE(jpbcy0)
     CASE(0)
     CASE(1)
        indices(2) = j
        DO i = ipbeg,ipend
           indices(1) =  i
           bcvalue1 = ZERO
           CALL HYPRE_StructMatrixSetValues(Lmatrix,indices,1,3,bcvalue1,ierr)
        END DO
     CASE(2)
        indices(2) = j
        DO i = ipbeg,ipend
           indices(1) =  i
           bcvalue1 = ZERO
           bcvalue2 = r2dxideta*(etacy(j)*etafy(j)+etacy(j)*etafy(j-1))
           CALL HYPRE_StructMatrixSetValues(Lmatrix,indices,1,3,bcvalue1,ierr)
           CALL HYPRE_StructMatrixSetValues(Lmatrix,indices,1,4,bcvalue2,ierr)
        END DO
     END SELECT
  END IF
  ! back
  IF(jpend==ny) THEN
     j=ny
     SELECT CASE(jpbcy1)
     CASE(0)
     CASE(1)
        indices(2) = j
        DO i = ipbeg,ipend
           indices(1) =  i
           bcvalue1 = ZERO
           CALL HYPRE_StructMatrixSetValues(Lmatrix,indices,1,4,bcvalue1,ierr)
        END DO
     CASE(2)
        indices(2) = j
        DO i = ipbeg,ipend
           indices(1) =  i
           bcvalue1 = ZERO
           bcvalue2 = r2dxideta*(etacy(j)*etafy(j-1)+etacy(j)*etafy(j))
           CALL HYPRE_StructMatrixSetValues(Lmatrix,indices,1,4,bcvalue1,ierr)
           CALL HYPRE_StructMatrixSetValues(Lmatrix,indices,1,3,bcvalue2,ierr)
        END DO
     END SELECT
  END IF

  CALL HYPRE_StructMatrixAssemble(Lmatrix,ierr)

  ! Step 4: setup the solution and rhs vectors
  CALL HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,bvector,ierr)
  CALL HYPRE_StructVectorInitialize(bvector,ierr)
  CALL HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,xvector,ierr)
  CALL HYPRE_StructVectorInitialize(xvector,ierr)

END SUBROUTINE setuphypre
