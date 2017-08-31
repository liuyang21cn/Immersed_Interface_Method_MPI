FUNCTION panelinPE(A,B,corner0,corner1)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(in):: A,B,corner0,corner1
  LOGICAL panelinPE
  INTEGER:: i
  INTEGER, DIMENSION(1:2):: m,n,l
  DOUBLE PRECISION, DIMENSION(1:2):: A2,B2,C2,P2
  DOUBLE PRECISION:: coord,y,x

  panelinPE=.FALSE.

  ! test if vertex A,B in box [corner0,corner1]
  IF((corner0(1)<=A(1).AND.corner1(1)>=A(1).AND.&
       corner0(2)<=A(2).AND.corner1(2)>=A(2)).OR.&
       (corner0(1)<=B(1).AND.corner1(1)>=B(1).AND.&
       corner0(2)<=B(2).AND.corner1(2)>=B(2))) THEN
     panelinPE=.TRUE.
     RETURN
  ENDIF

  !test if vertex point on edge of box [corner0,corner1]
  IF(A(1).EQ.B(1)) THEN
     IF((A(1).EQ.corner0(1)).OR.(A(1).EQ.corner1(1))) THEN
        panelinPE=.TRUE.
        RETURN
     ENDIF
  ENDIF
  IF(A(2).EQ.B(2)) THEN
     IF((A(2).EQ.corner0(2)).OR.(A(2).EQ.corner1(2))) THEN
        panelinPE=.TRUE.
        RETURN
     ENDIF
  ENDIF

  !test if panel cross processor but vertex not on it
  IF((A(1).NE.B(1)).AND.(A(2).NE.B(2))) THEN
     y=(B(2)-A(2))/(B(1)-A(1))*(corner0(1)-A(1))+A(2)
     IF((y.LT.corner1(2)).AND.(y.GT.corner0(2))) THEN
        panelinPE=.TRUE.
        RETURN
     ENDIF
     y=(B(2)-A(2))/(B(1)-A(1))*(corner1(1)-A(1))+A(2)
     IF((y.GT.corner0(2)).AND.(y.LT.corner1(2))) THEN
        panelinPE=.TRUE.
        RETURN
     ENDIF
     x=(corner0(2)-A(2))*(B(1)-A(1))/(B(2)-A(2))+A(1)
     IF((x.GT.corner0(1)).AND.(x.LT.corner1(1))) THEN
        panelinPE=.TRUE.
        RETURN
     ENDIF
     x=(corner1(2)-A(2))*(B(1)-A(1))/(B(2)-A(2))+A(1)
     IF((x.GT.corner0(1)).AND.(x.LT.corner1(1))) THEN
        panelinPE=.TRUE.
        RETURN
     ENDIF
  ENDIF

END FUNCTION panelinPE
