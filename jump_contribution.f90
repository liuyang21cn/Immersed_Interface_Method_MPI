SUBROUTINE jump_contribution(krk)

  USE para
  USE Lagrange
  USE Euler

  INTEGER:: index,flag
  INTEGER,INTENT(in):: krk
  DOUBLE PRECISION:: ucpdx,vcpdy,vcixfyf

  ucpdx = zero
  vcpdy = zero
  ! yf
  CALL interpolateuv(4)
  DO index = 1, n_yf_int
     CALL correction_interpolate(index,4)
     CALL dudv_surface(index,4)
     CALL correction_velocity(index,4,krk,zero)
  ENDDO

  ! yc
  CALL interpolateuv(3)
  DO index = 1, n_yc_int
     CALL correction_interpolate(index,3)
     CALL correction_strain(index,3)
     CALL dudv_surface(index,3)
     CALL jc_pressure(index,3)
     CALL correction_pressure(index,3,ucpdx)
     CALL correction_velocity(index,3,krk,ucpdx)
  ENDDO

  ! xf
  CALL interpolateuv(2)
  DO index = 1, n_xf_int
     CALL correction_interpolate(index,2)
     CALL dudv_surface(index,2)
     CALL correction_velocity(index,2,krk,zero)
  ENDDO

  ! xc
  CALL interpolateuv(1)
  DO index = 1, n_xc_int
     CALL correction_interpolate(index,1)
     CALL correction_strain(index,1)
     CALL dudv_surface(index,1)
     CALL jc_pressure(index,1)
     CALL correction_pressure(index,1,vcpdy)
     CALL correction_velocity(index,1,krk,vcpdy)
  ENDDO

END SUBROUTINE jump_contribution
