SUBROUTINE objUpdate(myobj_global,krk,objMove)

  USE para
  USE Lagrange

  INTEGER, INTENT(in):: myobj_global,krk
  DOUBLE PRECISION:: phi0,phi,phit,tc
  DOUBLE PRECISION, DIMENSION(1:9),INTENT(out)::objMove

  tc=pi/0.4d0
  IF(myobj_global==1) THEN

     objMove(3)=0.75*pi+0.25d0*pi*(SIN(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc)))
     objMove(1)=1.25d0*(COS(0.8d0*tmove)+1.0d0)*COS(pi/3.0d0)
     objMove(2)=1.25d0*(COS(0.8d0*tmove)+1.0d0)*SIN(pi/3.0d0)

     objMove(6)=0.25d0*pi*0.8d0*COS(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc))+ &
          0.25d0*pi*(SIN(0.8d0*tmove))*(EXP(-tmove/tc)/tc)
     objMove(4)=1.25d0*(-0.8d0*SIN(0.8d0*tmove))*COS(pi/3.0d0)
     objMove(5)=1.25d0*(-0.8d0*SIN(0.8d0*tmove))*SIN(pi/3.0d0)

     objMove(9)=-0.25d0*pi*0.64d0*SIN(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc))+ &
          0.5d0*pi*(0.8d0*COS(0.8d0*tmove))*(EXP(-tmove/tc)/tc)- &
          0.25d0*pi*(SIN(0.8d0*tmove))*(EXP(-tmove/tc)/tc/tc)
     objMove(7)=1.25d0*(-0.64d0*COS(0.8d0*tmove))*COS(pi/3.0d0)
     objMove(8)=1.25d0*(-0.64d0*COS(0.8d0*tmove))*SIN(pi/3.0d0)

  ELSEIF ( myobj_global==2 ) THEN

     objMove(3)=0.75d0*pi+0.25d0*pi*(SIN(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc)))
     objMove(1)=1.25d0*(COS(0.8d0*tmove)+1.0d0)*COS(pi/3.0d0)
     objMove(2)=-3.0d0+1.25d0*(COS(0.8d0*tmove)+1.0d0)*SIN(pi/3.0d0)

     objMove(6)=0.25d0*pi*0.8d0*COS(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc))+ &
          0.25d0*pi*(SIN(0.8d0*tmove))*(EXP(-tmove/tc)/tc)
     objMove(4)=1.25d0*(-0.8d0*SIN(0.8d0*tmove))*COS(pi/3.0d0)
     objMove(5)=1.25d0*(-0.8d0*SIN(0.8d0*tmove))*SIN(pi/3.0d0)

     objMove(9)=-0.25d0*pi*0.64d0*SIN(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc))+ &
          0.5d0*pi*(0.8d0*COS(0.8d0*tmove))*(EXP(-tmove/tc)/tc)- &
          0.25d0*pi*(SIN(0.8d0*tmove))*(EXP(-tmove/tc)/tc/tc)
     objMove(7)=1.25d0*(-0.64d0*COS(0.8d0*tmove))*COS(pi/3.0d0)
     objMove(8)=1.25d0*(-0.64d0*COS(0.8d0*tmove))*SIN(pi/3.0d0)

  ELSEIF ( myobj_global==3 ) THEN

     objMove(3)=0.75d0*pi+0.25d0*pi*(SIN(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc)))
     objMove(1)=-3.0d0+1.25d0*(COS(0.8d0*tmove)+1.0d0)*COS(pi/3.0d0)
     objMove(2)=1.25d0*(COS(0.8d0*tmove)+1.0d0)*SIN(pi/3.0d0)

     objMove(6)=0.25d0*pi*0.8d0*COS(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc))+ &
          0.25d0*pi*(SIN(0.8d0*tmove))*(EXP(-tmove/tc)/tc)
     objMove(4)=1.25d0*(-0.8d0*SIN(0.8d0*tmove))*COS(pi/3.0d0)
     objMove(5)=1.25d0*(-0.8d0*SIN(0.8d0*tmove))*SIN(pi/3.0d0)

     objMove(9)=-0.25d0*pi*0.64d0*SIN(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc))+ &
          0.5d0*pi*(0.8d0*COS(0.8d0*tmove))*(EXP(-tmove/tc)/tc)- &
          0.25d0*pi*(SIN(0.8d0*tmove))*(EXP(-tmove/tc)/tc/tc)
     objMove(7)=1.25d0*(-0.64d0*COS(0.8d0*tmove))*COS(pi/3.0d0)
     objMove(8)=1.25d0*(-0.64d0*COS(0.8d0*tmove))*SIN(pi/3.0d0)

  ELSEIF ( myobj_global==4 ) THEN

     objMove(3)=0.75d0*pi+0.25d0*pi*(SIN(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc)))
     objMove(1)=-3.0d0+1.25d0*(COS(0.8d0*tmove)+1.0d0)*COS(pi/3.0d0)
     objMove(2)=-3.0d0+1.25d0*(COS(0.8d0*tmove)+1.0d0)*SIN(pi/3.0d0)

     objMove(6)=0.25d0*pi*0.8d0*COS(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc))+ &
          0.25d0*pi*(SIN(0.8d0*tmove))*(EXP(-tmove/tc)/tc)
     objMove(4)=1.25d0*(-0.8d0*SIN(0.8d0*tmove))*COS(pi/3.0d0)
     objMove(5)=1.25d0*(-0.8d0*SIN(0.8d0*tmove))*SIN(pi/3.0d0)

     objMove(9)=-0.25d0*pi*0.64d0*SIN(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc))+ &
          0.5d0*pi*(0.8d0*COS(0.8d0*tmove))*(EXP(-tmove/tc)/tc)- &
          0.25d0*pi*(SIN(0.8d0*tmove))*(EXP(-tmove/tc)/tc/tc)
     objMove(7)=1.25d0*(-0.64d0*COS(0.8d0*tmove))*COS(pi/3.0d0)
     objMove(8)=1.25d0*(-0.64d0*COS(0.8d0*tmove))*SIN(pi/3.0d0)

  ENDIF

  ! 1024 flappers
  ! idx = MOD((myobj_global-1),32)
  ! jdx = INT((myobj_global-1)/32)
  !
  ! xc = 1.0d0 + 2.0d0*idx
  ! yc = 1.0d0 + 2.0d0*jdx
  ! tc=pi/0.4d0
  !
  ! objMove(3) = 0.75d0*pi+0.25d0*pi*(SIN(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc)))
  ! objMove(1) = xc + 1.25d0*(COS(0.8d0*tmove)+1.0d0)*COS(pi/3.0d0)
  ! objMove(2) = yc + 1.25d0*(COS(0.8d0*tmove)+1.0d0)*SIN(pi/3.0d0)
  !
  ! objMove(6) = 0.25d0*pi*0.8d0*COS(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc))+ &
  !      0.25d0*pi*(SIN(0.8d0*tmove))*(EXP(-tmove/tc)/tc)
  ! objMove(4) = 1.25d0*(-0.8d0*SIN(0.8d0*tmove))*COS(pi/3.0d0)
  ! objMove(5) = 1.25d0*(-0.8d0*SIN(0.8d0*tmove))*SIN(pi/3.0d0)
  !
  ! objMove(9) = -0.25d0*pi*0.64d0*SIN(0.8d0*tmove)*(1.0d0-EXP(-tmove/tc))+ &
  !      0.5d0*pi*(0.8d0*COS(0.8d0*tmove))*(EXP(-tmove/tc)/tc)- &
  !      0.25d0*pi*(SIN(0.8d0*tmove))*(EXP(-tmove/tc)/tc/tc)
  ! objMove(7) = 1.25d0*(-0.64d0*COS(0.8d0*tmove))*COS(pi/3.0d0)
  ! objMove(8) = 1.25d0*(-0.64d0*COS(0.8d0*tmove))*SIN(pi/3.0d0)
  !
  !
  ! ! object translating along a circle, up to 8 objects
  ! objMove(3)=(myobj_global-1)*pi/4.0d0+tmove*pi/4.0d0
  ! objMove(1)=2.0d0*COS(objMove(3))
  ! objMove(2)=2.0d0*SIN(objMove(3))
  !
  ! objMove(4)=-pi/2.0d0*SIN(objMove(3))
  ! objMove(5)= pi/2.0d0*COS(objMove(3))
  ! objMove(6)= pi/4.0d0
  ! objMove(7)= -pi*pi/8.0d0*COS(objMove(3))
  ! objMove(8)= -pi*pi/8.0d0*SIN(objMove(3))
  ! objMove(9)= 0.0d0


END SUBROUTINE objUpdate
