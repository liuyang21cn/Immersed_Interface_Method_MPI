SUBROUTINE allocateEuler
  USE para
  USE Euler

  ! u,v
  ALLOCATE(u(iubeg-nxghost:iuend+nxghost,jubeg-nyghost:juend+nyghost))
  ALLOCATE(ucc(iubeg-1:iuend+1,jpbeg-2:jpend+2))
  ALLOCATE(uce(ipbeg:ipend,jpbeg-1:jpend))
  ALLOCATE(uee(ipbeg-1:ipend,jpbeg-1:jpend))
  ALLOCATE(un(iubeg:iuend,jubeg:juend))
  ALLOCATE(rhsu1(iubeg:iuend,jubeg:juend))
  ALLOCATE(rhsu2(iubeg:iuend,jubeg:juend))
  ALLOCATE(rhsu3(iubeg:iuend,jubeg:juend))
  ALLOCATE(rhsu4(iubeg:iuend,jubeg:juend))
  ALLOCATE(rhsu(iubeg:iuend,jubeg:juend,1:4))

  ALLOCATE(v(ivbeg-nxghost:ivend+nxghost,jvbeg-nyghost:jvend+nyghost))
  ALLOCATE(vcc(ipbeg-2:ipend+2,jvbeg-1:jvend+1))
  ALLOCATE(vee(ipbeg-1:ipend,jpbeg-1:jpend))
  ALLOCATE(vec(ipbeg-1:ipend,jpbeg:jpend))
  ALLOCATE(vn(ivbeg:ivend,jvbeg:jvend))
  ALLOCATE(rhsv1(ivbeg:ivend,jvbeg:jvend))
  ALLOCATE(rhsv2(ivbeg:ivend,jvbeg:jvend))
  ALLOCATE(rhsv3(ivbeg:ivend,jvbeg:jvend))
  ALLOCATE(rhsv4(ivbeg:ivend,jvbeg:jvend))
  ALLOCATE(rhsv(ivbeg:ivend,jvbeg:jvend,1:4))

  u=zero
  ucc=zero
  uce=zero
  uee=zero
  un=zero
  rhsu1=zero
  rhsu2=zero
  rhsu3=zero
  rhsu4=zero
  rhsu=zero
  rhsv=zero

  v=zero
  vcc=zero
  vee=zero
  vec=zero
  vn=zero
  rhsv1=zero
  rhsv2=zero
  rhsv3=zero
  rhsv4=zero

  ! p,d
  ALLOCATE(p(ipbeg-nxghost:ipend+nxghost,jpbeg-nxghost:jpend+nxghost))
  ALLOCATE(rhsp(ipbeg:ipend,jpbeg:jpend))
  ALLOCATE(d(ipbeg-1:ipend+1,jpbeg-1:jpend+1))
  ALLOCATE(rhsph(ipbeg-1:ipend+1,jpbeg-1:jpend+1))
  ALLOCATE(ph(ipbeg-1:ipend+1,jpbeg-1:jpend+1))
  ALLOCATE(dn(ipbeg:ipend,jpbeg:jpend))
  p=zero
  rhsp=zero
  d=zero
  dn=zero

  IF(ipbeg.EQ.1) THEN
     idbeg=ipbeg
  ELSE
     idbeg=ipbeg-1
  ENDIF

  IF(ipend.EQ.nx) THEN
     idend=ipend
  ELSE
     idend=ipend+1
  ENDIF

  IF(jpbeg.EQ.1) THEN
     jdbeg=jpbeg
  ELSE
     jdbeg=jpbeg-1
  ENDIF

  IF(jpend.EQ.ny) THEN
     jdend=jpend
  ELSE
     jdend=jpend+1
  ENDIF
  ALLOCATE(ux(idbeg:idend,jdbeg:jdend))
  ALLOCATE(uy(idbeg:idend,jdbeg:jdend))
  ALLOCATE(vx(idbeg:idend,jdbeg:jdend))
  ALLOCATE(vy(idbeg:idend,jdbeg:jdend))
  ux=zero
  uy=zero
  vx=zero
  vy=zero

  ! Raycrossing indicator
  ALLOCATE(iou(iubeg:iuend,jubeg:juend))
  ALLOCATE(iov(ivbeg:ivend,jvbeg:jvend))
  ALLOCATE(iop(ipbeg-1:ipend+1,jpbeg-1:jpend+1))


END SUBROUTINE allocateEuler
