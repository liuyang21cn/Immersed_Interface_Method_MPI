MODULE Lagrange
  IMPLICIT NONE

  ! local object
  INTEGER:: nobj,nobj4proc
  INTEGER, DIMENSION(:), ALLOCATABLE:: object_list
  INTEGER, DIMENSION(:), ALLOCATABLE:: obj4proc,nvertex4proc,nvertexlocal,vertexindexlocal
  INTEGER, DIMENSION(:,:), ALLOCATABLE:: proc4obj
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: vertex
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: taox,taoy,taoxn,taoyn

  ! object movement
  DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE:: ObjPosition
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE:: us,vs,xs,ys,cd,cl
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE:: xsc,ysc,theta
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE:: xsct,ysct,thetat
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE:: xsctt,ysctt,thetatt

  ! panel-gridline intersections
  INTEGER,DIMENSION(:,:),ALLOCATABLE:: xc_int_vertex,yc_int_vertex,xf_int_vertex,yf_int_vertex
  INTEGER :: n_xc_int,n_yc_int,n_xf_int,n_yf_int ! number of irregular pts
  DOUBLE PRECISION, dimension(:,:),ALLOCATABLE:: ujcxc,ujcxf,vjcxc,vjcxf,pjcxc,pjcyc,ujcyc,ujcyf,vjcyc,vjcyf
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: xf_int_info, yf_int_info, xc_int_info, yc_int_info

  ! Cartesian jump conditions
  DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE:: ujc,vjc,pjc
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: dudxJC,dudyJC,ujcn,vjcn
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: dvdxJC,dvdyJC

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: d2udx2JC,d2udxyJC,d2udy2JC
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: d2vdx2JC,d2vdxyJC,d2vdy2JC

  INTEGER, DIMENSION(:), ALLOCATABLE:: n_middle_pt,middle_pt_index
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: dpxjcM,dpyjcM,jcp_rhs

END MODULE Lagrange
