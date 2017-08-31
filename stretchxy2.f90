MODULE stretchxy
  PRIVATE
  PUBLIC:: xi2x,x2xi,dxidx,dxdxi,d2xidx2,d2xdxi2,eta2y,y2eta,detady,dydeta,d2etady2,d2ydeta2
  DOUBLE PRECISION, PARAMETER:: pi=3.14159265358979323d0
  DOUBLE PRECISION, PARAMETER:: e =2.71828182845904524d0
  DOUBLE PRECISION, PARAMETER:: e0=1.0d0
  DOUBLE PRECISION, PARAMETER:: e1=3.0d0
  DOUBLE PRECISION:: a=2.0d0*pi/(e1-e0),b=pi*(e1+e0)/(e1-e0)
  DOUBLE PRECISION, PARAMETER:: CC=2.0d0 ! ratio for cc=2.8 is 87.8821
  DOUBLE PRECISION, PARAMETER:: AA=0.5d0/TANH(CC),BB=0.5d0

CONTAINS
  ! tanh
  
  ! x(xi)
  FUNCTION xi2x(ksi) RESULT(x)
    DOUBLE PRECISION, INTENT(in):: ksi
    DOUBLE PRECISION:: x

    x=AA*TANH(CC*ksi) + BB
  END FUNCTION xi2x

  ! xi(x)
  FUNCTION x2xi(x) RESULT(ksi)
    DOUBLE PRECISION, INTENT(in):: x
    DOUBLE PRECISION:: ksi

    ksi=atanh((x-BB)/AA)/CC
  END FUNCTION x2xi

  ! dx/dxi (x)
  FUNCTION dxdxi(x)
    DOUBLE PRECISION, INTENT(in):: x
    DOUBLE PRECISION:: dxdxi

    dxdxi=AA*CC*(1.0d0-(x-BB)*(x-BB)/AA/AA)
  END FUNCTION dxdxi

  ! d2x/dxi2 (x)
  FUNCTION d2xdxi2(x)
    DOUBLE PRECISION, INTENT(in):: x
    DOUBLE PRECISION:: d2xdxi2

    d2xdxi2=-2.0d0*CC*CC*(x-BB)*(1.0d0-(x-BB)*(x-BB)/AA/AA)
  END FUNCTION d2xdxi2

  ! dxi/dx (x)
  FUNCTION dxidx(x)
    DOUBLE PRECISION, INTENT(in):: x
    DOUBLE PRECISION:: dxidx

    dxidx=1.0d0/dxdxi(x)
    dxidx=AA/CC/(AA*AA-(x-BB)*(x-BB))
  END FUNCTION dxidx

  ! d2xi/dx2 (x)
  FUNCTION d2xidx2(x)
    DOUBLE PRECISION, INTENT(in):: x
    DOUBLE PRECISION:: d2xidx2

    d2xidx2=-d2xdxi2(x)/dxdxi(x)/dxdxi(x)/dxdxi(x)
    d2xidx2=2.0d0*AA*(x-BB)/CC/(AA*AA-(x-BB)*(x-BB))/(AA*AA-(x-BB)*(x-BB))
  END FUNCTION d2xidx2

  ! y(eta)
  FUNCTION eta2y(ita) RESULT(y)
    DOUBLE PRECISION, INTENT(in):: ita
    DOUBLE PRECISION:: y

    y=AA*TANH(CC*ita) + BB
  END FUNCTION eta2y

  ! eta(y)
  FUNCTION y2eta(y) RESULT(ita)
    DOUBLE PRECISION, INTENT(in):: y
    DOUBLE PRECISION:: ita

    ita=atanh((y-BB)/AA)/CC
  END FUNCTION y2eta

  ! dy/deta (y)
  FUNCTION dydeta(y)
    DOUBLE PRECISION, INTENT(in):: y
    DOUBLE PRECISION:: dydeta

    dydeta=AA*CC*(1.0d0-(y-BB)*(y-BB)/AA/AA)
    ! dydeta=AA/CC/(AA*AA-(y-BB)*(y-BB))
  END FUNCTION dydeta

  ! d2y/deta2 (y)
  FUNCTION d2ydeta2(y)
    DOUBLE PRECISION, INTENT(in):: y
    DOUBLE PRECISION:: d2ydeta2

    d2ydeta2=-2.0d0*CC*CC*(y-BB)*(1.0d0-(y-BB)*(y-BB)/AA/AA)
  END FUNCTION d2ydeta2

  ! deta/dy (y)
  FUNCTION detady(y)
    DOUBLE PRECISION, INTENT(in):: y
    DOUBLE PRECISION:: detady

    detady=1.0d0/dydeta(y)
    detady=AA/CC/(AA*AA-(y-BB)*(y-BB))
  END FUNCTION detady

  ! d2eta/dy2 (y)
  FUNCTION d2etady2(y)
    DOUBLE PRECISION, INTENT(in):: y
    DOUBLE PRECISION:: d2etady2

    d2etady2=-d2ydeta2(y)/dydeta(y)/dydeta(y)/dydeta(y)
    d2etady2=2.0d0*AA*(y-BB)/CC/(AA*AA-(y-BB)*(y-BB))/(AA*AA-(y-BB)*(y-BB))
  END FUNCTION d2etady2

END MODULE stretchxy
