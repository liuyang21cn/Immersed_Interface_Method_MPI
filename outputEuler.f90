SUBROUTINE outputEuler

  USE para
  USE Euler
  USE Lagrange
  CHARACTER(5) filesuffix

  CALL HYPRE_StructVectorDestroy(xvector,ierr)
  CALL HYPRE_StructVectorDestroy(bvector,ierr)
  CALL HYPRE_StructMatrixDestroy(Lmatrix,ierr)
  CALL HYPRE_StructGridDestroy(grid,ierr)

  DEALLOCATE(u)
  DEALLOCATE(ucc)
  DEALLOCATE(uce)
  DEALLOCATE(uee)
  DEALLOCATE(un)
  DEALLOCATE(rhsu1)
  DEALLOCATE(rhsu2)
  DEALLOCATE(rhsu3)
  DEALLOCATE(rhsu4)

  DEALLOCATE(v)
  DEALLOCATE(vcc)
  DEALLOCATE(vee)
  DEALLOCATE(vec)
  DEALLOCATE(vn)
  DEALLOCATE(rhsv1)
  DEALLOCATE(rhsv2)
  DEALLOCATE(rhsv3)
  DEALLOCATE(rhsv4)

  DEALLOCATE(p)
  DEALLOCATE(rhsp)
  DEALLOCATE(d)
  DEALLOCATE(dn)
  DEALLOCATE(ph)
  DEALLOCATE(rhsph)
  DEALLOCATE(iou)
  DEALLOCATE(iov)
  DEALLOCATE(iop)

  DEALLOCATE(xic)
  DEALLOCATE(xicx)
  DEALLOCATE(xicxx)
  DEALLOCATE(xc)
  DEALLOCATE(xif)
  DEALLOCATE(xifx)
  DEALLOCATE(xifxx)
  DEALLOCATE(xf)
  DEALLOCATE(etac)
  DEALLOCATE(etacy)
  DEALLOCATE(etacyy)
  DEALLOCATE(yc)
  DEALLOCATE(etaf)
  DEALLOCATE(etafy)
  DEALLOCATE(etafyy)
  DEALLOCATE(yf)
  DEALLOCATE(allcorner0)
  DEALLOCATE(allcorner1)
  DEALLOCATE(ProcsInfo)
  DEALLOCATE(ppp)
  DEALLOCATE(qqq)

  IF(myid==0 .AND. isingular.EQV..TRUE.) THEN
     DEALLOCATE(cd)
     DEALLOCATE(cl)
  ENDIF

  ! objects properties
  IF(objAllocated.EQV..TRUE.) THEN
     DEALLOCATE(xsc)
     DEALLOCATE(ysc)
     DEALLOCATE(theta)
     DEALLOCATE(xsct)
     DEALLOCATE(ysct)
     DEALLOCATE(thetat)
     DEALLOCATE(xsctt)
     DEALLOCATE(ysctt)
     DEALLOCATE(thetatt)
     DEALLOCATE(taox)
     DEALLOCATE(taoy)
     DEALLOCATE(us)
     DEALLOCATE(vs)
     DEALLOCATE(xs)
     DEALLOCATE(ys)
     objAllocated=.FALSE.
  ENDIF


END SUBROUTINE outputEuler
