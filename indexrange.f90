SUBROUTINE indexrange(n,np,id,ibeg,iend)
  INTEGER, INTENT(in):: n,np,id
  INTEGER, INTENT(out):: ibeg,iend
  INTEGER:: nloc,nr

  nloc = FLOOR(DBLE(n/np))
  nr = MOD(n,np)

  ibeg = id*nloc+1+MIN(id,nr)
  IF(id < nr) THEN
     nloc = nloc+1
  ENDIF
  iend = ibeg+nloc-1
  IF((id.EQ.np-1).OR.(iend > n)) THEN
     iend = n
  ENDIF

END SUBROUTINE indexrange
