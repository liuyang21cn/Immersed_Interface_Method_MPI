MODULE Euler
  IMPLICIT NONE

  ! grid and coord transformation
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: xic,etac,xif,etaf
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: xc,xf,yc,yf
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: xicx,xicxx,etacy,etacyy
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: xifx,xifxx,etafy,etafyy
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: ppp,qqq
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: allcorner0,allcorner1

  ! MPI information
  INTEGER, DIMENSION(:,:), ALLOCATABLE:: ProcsInfo

  ! u,v,w,p,d
  DOUBLE PRECISION, DIMENSION(:,:),  ALLOCATABLE:: ucc,uee,uce,u,un,rhsu1,rhsu2,rhsu3,rhsu4
  DOUBLE PRECISION, DIMENSION(:,:),  ALLOCATABLE:: vcc,vee,vec,v,vn,rhsv1,rhsv2,rhsv3,rhsv4
  DOUBLE PRECISION, DIMENSION(:,:),  ALLOCATABLE:: p,rhsp,d,dn,rhsph,ph
  DOUBLE PRECISION, DIMENSION(:,:),  ALLOCATABLE:: ux,uy,vx,vy
  DOUBLE PRECISION, DIMENSION(:,:,:),ALLOCATABLE:: rhsu,rhsv

  ! in or out indicator
  INTEGER*2, DIMENSION(:,:), ALLOCATABLE:: iou,iov,iop

END MODULE Euler
