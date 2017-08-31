SUBROUTINE Euler_reset(krk)

  USE Euler
  USE para

  INTEGER,INTENT(in):: krk

  d = zero
  rhsp = zero
  ux = zero
  uy = zero
  vx = zero
  vy = zero

  ucc=zero
  vcc=zero
  vee=zero
  uee=zero
  uce=zero
  vec=zero

  rhsph=zero
  ph=zero

  IF(krk.EQ.1) THEN
     rhsu1 = zero
     rhsv1 = zero
     dn = zero
     un = zero
     vn = zero
  ELSEIF ( krk.EQ.2 ) THEN
     rhsu2 = zero
     rhsv2 = zero
  ELSEIF ( krk.EQ.3 ) THEN
     rhsu3 = zero
     rhsv3 = zero
  ELSE
     rhsu4 = zero
     rhsv4 = zero
  ENDIF

END SUBROUTINE Euler_reset
