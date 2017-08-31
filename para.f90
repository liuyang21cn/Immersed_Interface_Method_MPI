MODULE para
  IMPLICIT NONE

  ! constants
  DOUBLE PRECISION, PARAMETER:: pi=3.14159265358979323d0
  DOUBLE PRECISION, PARAMETER:: e =2.71828182845904524d0
  DOUBLE PRECISION, PARAMETER:: ZERO=0.0d0
  DOUBLE PRECISION, PARAMETER:: HALF=0.5d0
  DOUBLE PRECISION, PARAMETER:: ONE=1.0d0
  DOUBLE PRECISION, PARAMETER:: TWO=2.0d0
  DOUBLE PRECISION, PARAMETER:: eps=1.0d-12

  ! (sub)domain and grid
  INTEGER:: nx,ny,nxghost,nyghost
  INTEGER:: ipbeg,ipend,jpbeg,jpend
  INTEGER:: iubeg,iuend,jubeg,juend
  INTEGER:: ivbeg,ivend,jvbeg,jvend
  INTEGER:: idbeg,idend,jdbeg,jdend
  INTEGER:: i,j,k,l,m,n
  INTEGER:: icfl,iread,iwrite,iplot,itout,ianimation
  LOGICAL:: isingular,object_move,JCsaved,objAllocated

  DOUBLE PRECISION:: Re,Re1
  DOUBLE PRECISION:: x0,x1,y0,y1,r2dxideta,xi0,xi1,eta0,eta1
  DOUBLE PRECISION:: dxi,halfdxi,twodxi,deta,halfdeta,twodeta
  DOUBLE PRECISION:: dx,halfdx,twodx,dy,halfdy,twody
  DOUBLE PRECISION:: MinDx,MinDy,MaxDx,MaxDy
  DOUBLE PRECISION:: dxmin,dymin,dxmax,dymax

  ! BCs
  INTEGER:: ipbcx0,ipbcx1,jpbcy0,jpbcy1
  INTEGER:: iubcx0,iubcx1,jubcy0,jubcy1
  INTEGER:: ivbcx0,ivbcx1,jvbcy0,jvbcy1

  ! time-marching
  INTEGER:: nstart,nend,nstep,nn
  DOUBLE PRECISION:: t,dt,cflc,cflv,dtcfl
  DOUBLE PRECISION:: dt0,dt1,dt2,tout,t0,tmove

  ! time counter
  DOUBLE PRECISION:: hypreTime, ghostTime, movingTime, pTime

  ! MPI
  INTEGER:: nprocs,px,py,ierr,comm2d
  INTEGER, DIMENSION(1:2):: dims, coords
  LOGICAL, DIMENSION(1:2):: periods
  INTEGER, DIMENSION(0:8):: nearproc ! 2d, there are 8 procs around selected one, if 3d, then 3*9-1 = 26
  INTEGER:: myid,left,right,front,back
  INTEGER:: uitype,ujtype,vitype,vjtype,pitype,pjtype

  ! hypre
  INTEGER:: maxit,relch,rlxtype,npre,npost
  INTEGER*8:: grid,stencil,Lmatrix,bvector,xvector,solver,precond
  DOUBLE PRECISION:: tol

END MODULE para
