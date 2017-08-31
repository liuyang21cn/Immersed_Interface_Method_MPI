SUBROUTINE initial

  USE para
  USE Euler
  USE Lagrange

  nstart=0
  t=ZERO
  t0=t

  ! initial fluid field
  p=ZERO
  d=ZERO
  dn=ZERO
  u=one
  v=ZERO

  ! initial all right hand side
  rhsp=ZERO
  rhsu1=ZERO
  rhsv1=ZERO
  rhsu2=ZERO
  rhsv2=ZERO
  rhsu3=ZERO
  rhsv3=ZERO
  rhsu4=ZERO
  rhsv4=ZERO

  ! initial time counter
  hypreTime = ZERO
  ghostTime = ZERO
  movingTime = ZERO
  pTime = ZERO

  ! initial velocity with jump conditions
  CALL velocity_initial

END SUBROUTINE initial
