MODULE stretchxy
  PRIVATE
  PUBLIC:: xi2x,x2xi,dxidx,dxdxi,d2xidx2,d2xdxi2,eta2y,y2eta,detady,dydeta,d2etady2,d2ydeta2
  DOUBLE PRECISION, PARAMETER:: pi=3.14159265358979323d0
  DOUBLE PRECISION, PARAMETER:: e =2.71828182845904524d0
  DOUBLE PRECISION, PARAMETER:: e0=1.0d0
  DOUBLE PRECISION, PARAMETER:: e1=3.0d0
  DOUBLE PRECISION, PARAMETER:: cc=1.0d0/3.0d0
  DOUBLE PRECISION:: a=2.0d0*pi/(e1-e0),b=pi*(e1+e0)/(e1-e0)

CONTAINS

  ! tan

  ! x(xi)
  FUNCTION xi2x(ksi) RESULT(x)
    DOUBLE PRECISION, INTENT(in):: ksi
    DOUBLE PRECISION:: x

    x=2.0d0*tan(pi/4.0d0*ksi)
  END FUNCTION xi2x

  ! xi(x)
  FUNCTION x2xi(x) RESULT(ksi)
    DOUBLE PRECISION, INTENT(in):: x
    DOUBLE PRECISION:: ksi

    ksi =4.0d0/pi*atan(x/2.0d0)
  END FUNCTION x2xi

  ! dx/dxi (x)
  FUNCTION dxdxi(x)
    DOUBLE PRECISION, INTENT(in):: x
    DOUBLE PRECISION:: dxdxi

    dxdxi=2.0d0*pi/4.0d0*(1.0d0+x*x/4.0d0)
  END FUNCTION dxdxi

  ! d2x/dxi2 (x)
  FUNCTION d2xdxi2(x)
    DOUBLE PRECISION, INTENT(in):: x
    DOUBLE PRECISION:: d2xdxi2

    d2xdxi2=pi*pi/8.0d0*(1.0d0+x*x/4.0d0)*x
  END FUNCTION d2xdxi2

  ! dxi/dx (x)
  FUNCTION dxidx(x)
    DOUBLE PRECISION, INTENT(in):: x
    DOUBLE PRECISION:: dxidx

    dxidx=1.0d0/dxdxi(x)
    dxidx=8.0d0/pi/(4.0d0+x*x)
  END FUNCTION dxidx

  ! d2xi/dx2 (x)
  FUNCTION d2xidx2(x)
    DOUBLE PRECISION, INTENT(in):: x
    DOUBLE PRECISION:: d2xidx2

    d2xidx2=-d2xdxi2(x)/dxdxi(x)/dxdxi(x)/dxdxi(x)
    d2xidx2=-16.0d0*x/pi/(4.0d0+x*x)/(4.0d0+x*x)
  END FUNCTION d2xidx2

  ! y(eta)
  FUNCTION eta2y(ita) RESULT(y)
    DOUBLE PRECISION, INTENT(in):: ita
    DOUBLE PRECISION:: y

    y=2.0d0*tan(pi/4.0d0*ita)
  END FUNCTION eta2y

  ! eta(y)
  FUNCTION y2eta(y) RESULT(ita)
    DOUBLE PRECISION, INTENT(in):: y
    DOUBLE PRECISION:: ita

    ita =4.0d0/pi*atan(y/2.0d0)
  END FUNCTION y2eta

  ! dy/deta (y)
  FUNCTION dydeta(y)
    DOUBLE PRECISION, INTENT(in):: y
    DOUBLE PRECISION:: dydeta

    dydeta=2.0d0*pi/4.0d0*(1.0d0+y*y/4.0d0)
  END FUNCTION dydeta

  ! d2y/deta2 (y)
  FUNCTION d2ydeta2(y)
    DOUBLE PRECISION, INTENT(in):: y
    DOUBLE PRECISION:: d2ydeta2

    d2ydeta2=pi*pi/8.0d0*(1.0d0+y*y/4.0d0)*y
  END FUNCTION d2ydeta2

  ! deta/dy (y)
  FUNCTION detady(y)
    DOUBLE PRECISION, INTENT(in):: y
    DOUBLE PRECISION:: detady

    detady=1.0d0/dydeta(y)
    detady=8.0d0/pi/(4.0d0+y*y)
  END FUNCTION detady

  ! d2eta/dy2 (y)
  FUNCTION d2etady2(y)
    DOUBLE PRECISION, INTENT(in):: y
    DOUBLE PRECISION:: d2etady2

    d2etady2=-d2ydeta2(y)/dydeta(y)/dydeta(y)/dydeta(y)
    d2etady2=-16.0d0*y/pi/(4.0d0+y*y)/(4.0d0+y*y)
  END FUNCTION d2etady2

END MODULE stretchxy
