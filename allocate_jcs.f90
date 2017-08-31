SUBROUTINE allocate_jcs(krk)

  USE para
  USE Lagrange
  USE Euler

  INTEGER:: krk

  ! principle jump conditions
  ALLOCATE(ujc(1:6, 1:nvertex4proc(nobj4proc)))
  ALLOCATE(vjc(1:6, 1:nvertex4proc(nobj4proc)))
  ALLOCATE(pjc(1:8, 1:nvertex4proc(nobj4proc)))
  ujc=zero
  vjc=zero
  pjc=zero
  
  ! allocate ujcn,vjcn to store dudn|+ dvdn|+ for draglift
  IF(krk==0) THEN
     ALLOCATE(ujcn(1:nvertex4proc(nobj4proc)))
     ALLOCATE(vjcn(1:nvertex4proc(nobj4proc)))
     ALLOCATE(taoxn(1:nvertex4proc(nobj4proc)))
     ALLOCATE(taoyn(1:nvertex4proc(nobj4proc)))
     JCsaved = .TRUE.
  ENDIF
  IF(krk==3) THEN
     IF(object_move.EQV..TRUE.) THEN
        IF(JCsaved.EQV..TRUE.) THEN
           DEALLOCATE(ujcn)
           DEALLOCATE(vjcn)
           DEALLOCATE(taoxn)
           DEALLOCATE(taoyn)
        ENDIF
        ALLOCATE(ujcn(1:nvertex4proc(nobj4proc)))
        ALLOCATE(vjcn(1:nvertex4proc(nobj4proc)))
        ALLOCATE(taoxn(1:nvertex4proc(nobj4proc)))
        ALLOCATE(taoyn(1:nvertex4proc(nobj4proc)))
        ujcn=zero
        vjcn=zero
        taoxn=zero
        taoyn=zero
        JCsaved = .TRUE.
     ENDIF
  ENDIF

  ! xf
  ALLOCATE(ujcxf(1:6,1:n_xf_int))
  ALLOCATE(vjcxf(1:6,1:n_xf_int))
  ! xc
  ALLOCATE(vjcxc(1:6,1:n_xc_int))
  ALLOCATE(ujcxc(1:6,1:n_xc_int))
  ALLOCATE(pjcxc(1:8,1:n_xc_int))
  ! yf
  ALLOCATE(ujcyf(1:6,1:n_yf_int))
  ALLOCATE(vjcyf(1:6,1:n_yf_int))
  ! yc
  ALLOCATE(vjcyc(1:6,1:n_yc_int))
  ALLOCATE(ujcyc(1:6,1:n_yc_int))
  ALLOCATE(pjcyc(1:8,1:n_yc_int))

  ujcxf=zero
  vjcxf=zero

  ujcxc=zero
  vjcxc=zero
  pjcxc=zero

  ujcyf=zero
  vjcyf=zero

  vjcyc=zero
  ujcyc=zero
  pjcyc=zero
END SUBROUTINE allocate_jcs
