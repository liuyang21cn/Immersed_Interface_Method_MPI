SUBROUTINE run_output

  USE mpi
  USE para
  USE Euler
  USE Lagrange

  DOUBLE PRECISION,DIMENSION(ipbeg:ipend,jpbeg:jpend)::rhspsubarray,psubarray,usubarray,vsubarray,dsubarray
  INTEGER:: pfiletype,ufiletype,vfiletype,dfiletype,dfh,pfh,ufh,vfh, info, np,nu,nv
  INTEGER:: rhspfiletype,rhspfh
  INTEGER,DIMENSION(1:2)::psizes,psubsizes,pstarts
  INTEGER,DIMENSION(1:2)::dsizes,dsubsizes,dstarts
  INTEGER,DIMENSION(1:2)::rhspsizes,rhspsubsizes,rhspstarts
  INTEGER,DIMENSION(1:2)::usizes,usubsizes,ustarts
  INTEGER,DIMENSION(1:2)::vsizes,vsubsizes,vstarts
  INTEGER(kind=MPI_OFFSET_KIND),PARAMETER::zero_off=0

  INTEGER:: ih,jh,iref,igap,iend,jgap,jend,ii,jj,mark1,mark2,iexact
  DOUBLE PRECISION:: foo,unorm,vnorm,pnorm1,pnorm2,cp1,cp2,xcr,ycr,d0,pmax,pmin
  DOUBLE PRECISION,DIMENSION(ipbeg-1:ipend+1,jpbeg-1:jpend+1)::au,av,ap,au1,av1,ap1
  DOUBLE PRECISION:: o1,o2,r1,r2,a,b,rr,norm_u,norm_v,norm1_p,norm2_p
  INTEGER:: iumax,jumax,jvmax,ivmax,ipmax,jpmax

  np=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  psizes=(/nx,ny/)
  psubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  pstarts=(/ipbeg-1,jpbeg-1/)
  psubarray=p(ipbeg:ipend,jpbeg:jpend)

  np=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  rhspsizes=(/nx,ny/)
  rhspsubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  rhspstarts=(/ipbeg-1,jpbeg-1/)
  rhspsubarray=rhsp(ipbeg:ipend,jpbeg:jpend)/dxi/dxi

  nu=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  usizes=(/nx,ny/)
  usubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  ustarts=(/ipbeg-1,jpbeg-1/)
  usubarray=ucc(ipbeg:ipend,jpbeg:jpend)

  nv=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  vsizes=(/nx,ny/)
  vsubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  vstarts=(/ipbeg-1,jpbeg-1/)
  vsubarray=vcc(ipbeg:ipend,jpbeg:jpend)


  ! output pressure
  CALL MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       pfiletype,ierr)
  CALL MPI_Type_commit(pfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/p.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
  CALL MPI_FILE_SET_VIEW(pfh,zero_off,MPI_DOUBLE_PRECISION,pfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(pfh,psubarray,np,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(pfh,ierr)

  CALL MPI_TYPE_CREATE_SUBARRAY(2, rhspsizes, rhspsubsizes, rhspstarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       rhspfiletype,ierr)
  CALL MPI_Type_commit(rhspfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/rhsp.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,rhspfh,ierr)
  CALL MPI_FILE_SET_VIEW(rhspfh,zero_off,MPI_DOUBLE_PRECISION,rhspfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(rhspfh,rhspsubarray,np,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(rhspfh,ierr)

  ! streamFunction
  CALL MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       pfiletype,ierr)
  CALL MPI_Type_commit(pfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/ph.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
  CALL MPI_FILE_SET_VIEW(pfh,zero_off,MPI_DOUBLE_PRECISION,pfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(pfh,-ph(ipbeg:ipend,jpbeg:jpend),np,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(pfh,ierr)

  CALL MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       pfiletype,ierr)
  CALL MPI_Type_commit(pfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/wo.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
  CALL MPI_FILE_SET_VIEW(pfh,zero_off,MPI_DOUBLE_PRECISION,pfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(pfh,rhsph(ipbeg:ipend,jpbeg:jpend)/dxi/dxi,np,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(pfh,ierr)

  ! output ucc
  CALL MPI_TYPE_CREATE_SUBARRAY(2, usizes, usubsizes, ustarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       ufiletype,ierr)
  CALL MPI_Type_commit(ufiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/uc.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,ufh,ierr)
  CALL MPI_FILE_SET_VIEW(ufh,zero_off,MPI_DOUBLE_PRECISION,ufiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(ufh,usubarray,nu,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(ufh,ierr)

  ! output vcc
  CALL MPI_TYPE_CREATE_SUBARRAY(2, vsizes, vsubsizes, vstarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       vfiletype,ierr)
  CALL MPI_Type_commit(vfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/vc.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,vfh,ierr)
  CALL MPI_FILE_SET_VIEW(vfh,zero_off,MPI_DOUBLE_PRECISION,vfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(vfh,vsubarray,nv,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(vfh,ierr)

  ! u
  nu=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  usizes=(/nx,ny/)
  usubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  ustarts=(/ipbeg-1,jpbeg-1/)
  usubarray=u(ipbeg:ipend,jpbeg:jpend)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, usizes, usubsizes, ustarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       ufiletype,ierr)
  CALL MPI_Type_commit(ufiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/u.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,ufh,ierr)
  CALL MPI_FILE_SET_VIEW(ufh,zero_off,MPI_DOUBLE_PRECISION,ufiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(ufh,usubarray,nu,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(ufh,ierr)

  ! uee
  nu=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  usizes=(/nx,ny/)
  usubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  ustarts=(/ipbeg-1,jpbeg-1/)
  usubarray=uee(ipbeg:ipend,jpbeg:jpend)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, usizes, usubsizes, ustarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       ufiletype,ierr)
  CALL MPI_Type_commit(ufiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/uee.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,ufh,ierr)
  CALL MPI_FILE_SET_VIEW(ufh,zero_off,MPI_DOUBLE_PRECISION,ufiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(ufh,usubarray,nu,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(ufh,ierr)

  ! uce
  nu=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  usizes=(/nx,ny/)
  usubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  ustarts=(/ipbeg-1,jpbeg-1/)
  usubarray=uce(ipbeg:ipend,jpbeg:jpend)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, usizes, usubsizes, ustarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       ufiletype,ierr)
  CALL MPI_Type_commit(ufiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/uce.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,ufh,ierr)
  CALL MPI_FILE_SET_VIEW(ufh,zero_off,MPI_DOUBLE_PRECISION,ufiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(ufh,usubarray,nu,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(ufh,ierr)

  ! vee
  nv=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  vsizes=(/nx,ny/)
  vsubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  vstarts=(/ipbeg-1,jpbeg-1/)
  vsubarray=vee(ipbeg:ipend,jpbeg:jpend)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, vsizes, vsubsizes, vstarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       vfiletype,ierr)
  CALL MPI_Type_commit(vfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/vee.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,vfh,ierr)
  CALL MPI_FILE_SET_VIEW(vfh,zero_off,MPI_DOUBLE_PRECISION,vfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(vfh,vsubarray,nv,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(vfh,ierr)

  ! v
  nv=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  vsizes=(/nx,ny/)
  vsubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  vstarts=(/ipbeg-1,jpbeg-1/)
  vsubarray=v(ipbeg:ipend,jpbeg:jpend)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, vsizes, vsubsizes, vstarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       vfiletype,ierr)
  CALL MPI_Type_commit(vfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/v.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,vfh,ierr)
  CALL MPI_FILE_SET_VIEW(vfh,zero_off,MPI_DOUBLE_PRECISION,vfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(vfh,vsubarray,nv,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(vfh,ierr)

  ! vec
  nv=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  vsizes=(/nx,ny/)
  vsubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  vstarts=(/ipbeg-1,jpbeg-1/)
  vsubarray=vec(ipbeg:ipend,jpbeg:jpend)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, vsizes, vsubsizes, vstarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       vfiletype,ierr)
  CALL MPI_Type_commit(vfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/vec.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,vfh,ierr)
  CALL MPI_FILE_SET_VIEW(vfh,zero_off,MPI_DOUBLE_PRECISION,vfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(vfh,vsubarray,nv,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(vfh,ierr)

  ! iop
  np=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  psizes=(/nx,ny/)
  psubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  pstarts=(/ipbeg-1,jpbeg-1/)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
       MPI_ORDER_FORTRAN,MPI_INTEGER2,pfiletype,ierr)
  CALL MPI_Type_commit(pfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/iop.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
  CALL MPI_FILE_SET_VIEW(pfh,zero_off,MPI_INTEGER2,pfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(pfh,iop(ipbeg:ipend,jpbeg:jpend),np,MPI_INTEGER2,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(pfh,ierr)

  ! iou
  np=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  psizes=(/nx,ny/)
  psubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  pstarts=(/ipbeg-1,jpbeg-1/)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
       MPI_ORDER_FORTRAN,MPI_INTEGER2,pfiletype,ierr)
  CALL MPI_Type_commit(pfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/iou.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
  CALL MPI_FILE_SET_VIEW(pfh,zero_off,MPI_INTEGER2,pfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(pfh,iou(ipbeg:ipend,jpbeg:jpend),np,MPI_INTEGER2,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(pfh,ierr)

  ! iov
  np=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  psizes=(/nx,ny/)
  psubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  pstarts=(/ipbeg-1,jpbeg-1/)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
       MPI_ORDER_FORTRAN,MPI_INTEGER2,pfiletype,ierr)
  CALL MPI_Type_commit(pfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/iov.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
  CALL MPI_FILE_SET_VIEW(pfh,zero_off,MPI_INTEGER2,pfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(pfh,iov(ipbeg:ipend,jpbeg:jpend),np,MPI_INTEGER2,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(pfh,ierr)

  ! ux
  psubarray=ux(ipbeg:ipend,jpbeg:jpend)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       pfiletype,ierr)
  CALL MPI_Type_commit(pfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/ux.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
  CALL MPI_FILE_SET_VIEW(pfh,zero_off,MPI_DOUBLE_PRECISION,pfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(pfh,psubarray,np,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(pfh,ierr)

  ! uy
  psubarray=uy(ipbeg:ipend,jpbeg:jpend)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       pfiletype,ierr)
  CALL MPI_Type_commit(pfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/uy.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
  CALL MPI_FILE_SET_VIEW(pfh,zero_off,MPI_DOUBLE_PRECISION,pfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(pfh,psubarray,np,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(pfh,ierr)

  ! d
  np=(ipend-ipbeg+1)*(jpend-jpbeg+1)
  dsizes=(/nx,ny/)
  dsubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
  dstarts=(/ipbeg-1,jpbeg-1/)
  dsubarray=d(ipbeg:ipend,jpbeg:jpend)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, dsizes, dsubsizes, dstarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       dfiletype,ierr)
  CALL MPI_Type_commit(dfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/d.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,dfh,ierr)
  CALL MPI_FILE_SET_VIEW(dfh,zero_off,MPI_DOUBLE_PRECISION,dfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(dfh,dsubarray,np,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(dfh,ierr)

  ! vx
  psubarray=vx(ipbeg:ipend,jpbeg:jpend)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       pfiletype,ierr)
  CALL MPI_Type_commit(pfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/vx.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
  CALL MPI_FILE_SET_VIEW(pfh,zero_off,MPI_DOUBLE_PRECISION,pfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(pfh,psubarray,np,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(pfh,ierr)

  ! vy
  psubarray=vy(ipbeg:ipend,jpbeg:jpend)
  CALL MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
       pfiletype,ierr)
  CALL MPI_Type_commit(pfiletype,ierr)
  CALL MPI_Info_create(info,ierr)
  CALL MPI_Info_set(info,'collective_buffering','true',ierr)
  CALL MPI_File_open(comm2d,'OUTPUT/vy.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
  CALL MPI_FILE_SET_VIEW(pfh,zero_off,MPI_DOUBLE_PRECISION,pfiletype,'native',MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_ALL(pfh,psubarray,np,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  CALL MPI_FILE_CLOSE(pfh,ierr)


  o1=1.0d0
  o2=-1.0d0
  r1=0.5d0
  r2=2.0d0
  a=(o2*r2*r2-o1*r1*r1)/(r2*r2-r1*r1)
  b=(o1-o2)*r1*r1*r2*r2/(r2*r2-r1*r1)

  mark1=0
  mark2=0
  d0=0.0d0
  DO j=jpbeg,jpend
     DO i=ipbeg,ipend
        xcr=xc(i)-d0*SIN(t)
        ycr=yc(j)-d0*SIN(t)
        IF((xcr*xcr+ycr*ycr).GT.r1*r1*(1.0d0+0.01d0)) THEN
           IF(mark1.EQ.0) THEN
              cp1=-(0.5d0*a*a*(xcr*xcr+ycr*ycr)-&
                   0.5d0*b*b/(xcr*xcr+ycr*ycr)+&
                   a*b*LOG(xcr*xcr+ycr*ycr))-&
                   d0*SIN(t)*(xc(i)+yc(j))+p(i,j)
              mark1=1
           ENDIF
           au(i,j)=(a+b/(xcr*xcr+ycr*ycr))*(-ycr)+d0*COS(t)
           av(i,j)=(a+b/(xcr*xcr+ycr*ycr))*( xcr)+d0*COS(t)
           ap(i,j)=(0.5d0*a*a*(xcr*xcr+ycr*ycr)-0.5d0*b*b/(xcr*xcr+ycr*ycr)+&
                a*b*LOG(xcr*xcr+ycr*ycr))+d0*SIN(t)*(xc(i)+yc(j))+cp1
        ELSE
           au(i,j)=ucc(i,j)
           av(i,j)=vcc(i,j)
           ap(i,j)=p(i,j)
        ENDIF
     ENDDO
  ENDDO

  iumax=0
  jumax=0
  jvmax=0
  ivmax=0
  ipmax=0
  jpmax=0
  unorm=0.0d0
  vnorm=0.0d0
  pnorm1=0.0d0
  pnorm2=0.0d0
  pmax=0.0d0
  pmin=0.0d0
  DO j=jpbeg,jpend
     DO i=ipbeg,ipend
        au1(i,j)=ucc(i,j)-au(i,j)
        av1(i,j)=vcc(i,j)-av(i,j)
        ap1(i,j)=p(i,j)-ap(i,j)
        IF(unorm.GT.ABS(au1(i,j))) THEN
           iumax=iumax
           jumax=jumax
        ELSE
           iumax=i
           jumax=j
        ENDIF
        IF(vnorm.GT.ABS(av1(i,j))) THEN
           ivmax=ivmax
           jvmax=jvmax
        ELSE
           ivmax=i
           jvmax=j
        ENDIF
        IF(pnorm1.GT.ABS(ap1(i,j))) THEN
           ipmax=ipmax
           jpmax=jpmax
        ELSE
           ipmax=i
           jpmax=j
        ENDIF
        unorm=MAX(unorm,ABS(au1(i,j)))
        vnorm=MAX(vnorm,ABS(av1(i,j)))
        pnorm1=MAX(pnorm1,ABS(ap1(i,j)))
        pmax=MAX(pmax,ap1(i,j))
        pmin=MIN(pmin,ap1(i,j))
     ENDDO
  ENDDO
  pnorm2=0.5d0*(pmax-pmin)
  CALL MPI_Allreduce(unorm, norm_u, 1,MPI_DOUBLE_PRECISION,MPI_MAX,comm2d,ierr)
  CALL MPI_Allreduce(vnorm, norm_V, 1,MPI_DOUBLE_PRECISION,MPI_MAX,comm2d,ierr)
  CALL MPI_Allreduce(pnorm1,norm1_p,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm2d,ierr)
  CALL MPI_Allreduce(pnorm2,norm2_p,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm2d,ierr)

  IF(myid.EQ.0) THEN
     WRITE(*,*)
     WRITE(*,*)'unorm  = ',norm_u
     WRITE(*,*)'vnorm  = ',norm_V
     WRITE(*,*)'pnorm1 = ',norm1_p
     WRITE(*,*)'pnorm2 = ',norm2_p
  ENDIF

END SUBROUTINE run_output
