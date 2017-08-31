OBJ = allocateEuler.o Cartesiangrid.o cfl.o coords2rank.o\
      Euler.o exchange4ghost.o exchangeobject.o draglift.o objUpdate.o\
      indexrange.o initial.o interpolateuv.o Lagrange.o loadobject.o locateobject.o\
      outputEuler.o panelinPE.o para.o objectInitial.o\
      readpara.o rhs4p.o rk.o setuphypre.o startMPI.o stretchxy.o old_save.o\
      ubc.o vbc.o pbc.o rhs4uv.o updateuv.o run_output.o\
      velocity_initial.o Raycrossing.o surface_property.o jc_firstsecond.o interpolate.o\
      principle_jc_p.o allocate_jcs.o mac_distribute.o correction_interpolate.o jump_contribution.o\
 			uv_strain.o correction_strain.o dudv_surface.o correction_velocity.o pressure_distribute.o \
			correction_pressure.o velocity_reset.o free_jcs.o jc_pressure.o solvePoisson.o\
			rhs4p_clean.o Euler_reset.o Divergence_reset.o surface_move.o streamFunction.o outputObj.o\
			main.o
MOD1 = para.mod Euler.mod
MOD2 = para.mod Lagrange.mod
MOD3 = para.mod Euler.mod Lagrange.mod

FC = mpif90
FLAG = -g -O2
INC =

# Miguel
HYPRE = /home/yang/Downloads/hypre-2.11.2/src/hypre
#MBP
# HYPRE = /Users/YangLiu/Desktop/hypre-2.11.0/src/hypre
#iMac
#HYPRE = /Users/yang/Desktop/hypre-2.11.2/src/hypre
LIB = -L$(HYPRE)/lib -lHYPRE

# maneframe
# FC = mpif90
# FLAG = -g -O2 -fopenmp
# INC =
# LIB = -L$(HYPRE)/lib -lHYPRE
#
# HYPRE_DIR  = /grid/software/hypre/2.11.0/gcc-4.9.1
# BLASLIBS   = -L/grid/software/AtlasLib/3.10.2/gcc-4.9.1/lib -lf77blas -latlas
# LAPACKLIBS = -L/grid/software/AtlasLib/3.10.2/gcc-4.9.1/lib  -llapack -lcblas -latlas
# LIB        = -L$(HYPRE_DIR)/lib -lHYPRE -lm -lgfortran $(LAPACKLIBS) $(BLASLIBS)

main: $(OBJ)
	$(FC) $(FLAG) -o iim2dmpi.exe $(OBJ) $(LIB)
#	rm *.mod
#	rm *.o

# modules
para.mod: para.o para.f90
	$(FC) -c $(FLAG) para.f90
para.o: para.f90
	$(FC) -c $(FLAG) para.f90
Euler.mod: para.mod Euler.o Euler.f90
	$(FC) -c $(FLAG) Euler.f90
Euler.o: para.mod Euler.f90
	$(FC) -c $(FLAG) Euler.f90
Lagrange.mod: para.mod Lagrange.o Lagrange.f90
	$(FC) -c $(FLAG) Lagrange.f90
Lagrange.o: para.mod Lagrange.f90
	$(FC) -c $(FLAG) Lagrange.f90
stretchxy.mod: stretchxy.o stretchxy.f90
	$(FC) -c $(FLAG) stretchxy.f90
stretchxy.o: stretchxy.f90
	$(FC) -c $(FLAG) stretchxy.f90

# dependencies
allocateEuler.o: $(MOD1) indexrange.f90 allocateEuler.f90
	$(FC) -c $(FLAG) allocateEuler.f90
interpolate.o: interpolate.f90
	$(FC) -c $(FLAG) interpolate.f90
Cartesiangrid.o: $(MOD1) stretchxy.mod Cartesiangrid.f90
	$(FC) -c $(FLAG) Cartesiangrid.f90
cfl.o: $(MOD1) cfl.f90
	$(FC) -c $(FLAG) cfl.f90
coords2rank.o: para.mod coords2rank.f90
	$(FC) -c $(FLAG) coords2rank.f90
exchange4ghost.o: $(MOD1) exchange4ghost.f90
	$(FC) -c $(FLAG) exchange4ghost.f90
exchangeobject.o: $(MOD3) panelinPE.f90 exchangeobject.f90 objUpdate.f90
	$(FC) -c $(FLAG) exchangeobject.f90
outputObj.o: $(MOD2) outputObj.f90
	$(FC) -c $(FLAG) outputObj.f90
objectInitial.o: $(MOD2) objectInitial.f90
	$(FC) -c $(FLAG) objectInitial.f90
streamFunction.o: $(MOD3) interpolateuv.f90 objUpdate.f90 streamFunction.f90
	$(FC) -c $(FLAG) streamFunction.f90
vbc.o: $(MOD1) vbc.f90
	$(FC) -c $(FLAG) vbc.f90
ubc.o: $(MOD1) ubc.f90
	$(FC) -c $(FLAG) ubc.f90
pbc.o: $(MOD1) pbc.f90
	$(FC) -c $(FLAG) pbc.f90
allocate_jcs.o: $(MOD3) allocate_jcs.f90
	$(FC) -c $(FLAG) allocate_jcs.f90
jc_firstsecond.o: $(MOD3) jc_firstsecond.f90 interpolate.f90
	$(FC) -c $(FLAG) jc_firstsecond.f90
principle_jc_p.o: $(MOD3) principle_jc_p.f90 interpolate.f90
	$(FC) -c $(FLAG) principle_jc_p.f90
surface_property.o: $(MOD2) surface_property.f90 objUpdate.f90
	$(FC) -c $(FLAG) surface_property.f90
objUpdate.o: $(MOD3) objUpdate.f90
	$(FC) -c $(FLAG) objUpdate.f90
surface_move.o: $(MOD3) exchangeobject.f90 surface_move.f90 objUpdate.f90
	$(FC) -c $(FLAG) surface_move.f90
uv_strain.o: $(MOD3) uv_strain.f90
	$(FC) -c $(FLAG) uv_strain.f90
correction_strain.o: $(MOD3) correction_strain.f90
	$(FC) -c $(FLAG) correction_strain.f90
mac_distribute.o: $(MOD3) mac_distribute.f90
	$(FC) -c $(FLAG) mac_distribute.f90
pressure_distribute.o: $(MOD3) pressure_distribute.f90
	$(FC) -c $(FLAG) pressure_distribute.f90
correction_pressure.o: $(MOD3) correction_pressure.f90
	$(FC) -c $(FLAG) correction_pressure.f90
dudv_surface.o: $(MOD3) dudv_surface.f90 interpolate.f90
	$(FC) -c $(FLAG) dudv_surface.f90
jc_pressure.o: $(MOD3) jc_pressure.f90
	$(FC) -c $(FLAG) jc_pressure.f90
free_jcs.o: $(MOD3) free_jcs.f90
	$(FC) -c $(FLAG) free_jcs.f90
correction_velocity.o: $(MOD3) correction_velocity.f90
	$(FC) -c $(FLAG) correction_velocity.f90
correction_interpolate.o: $(MOD3) stretchxy.mod correction_interpolate.f90
	$(FC) -c $(FLAG) correction_interpolate.f90
interpolateuv.o: $(MOD3) interpolateuv.f90
	$(FC) -c $(FLAG) interpolateuv.f90
jump_contribution.o: $(MOD3) jump_contribution.f90 correction_interpolate.f90 interpolateuv.f90\
	 													 correction_strain.f90 uv_strain.f90 jc_pressure.f90\
														 correction_pressure.f90 dudv_surface.f90 correction_velocity.f90
	$(FC) -c $(FLAG) jump_contribution.f90
velocity_reset.o: $(MOD3) velocity_reset.f90
	$(FC) -c $(FLAG) velocity_reset.f90
velocity_initial.o: $(MOD3) velocity_initial.f90 solvePoisson.f90 exchange4ghost.f90 principle_jc_p.f90 \
														allocate_jcs.f90 ubc.f90 vbc.f90 pbc.f90 interpolateuv.f90 surface_property.f90\
														jc_firstsecond.f90 mac_distribute.f90 Raycrossing.f90 correction_pressure.f90\
														pressure_distribute.f90 velocity_reset.f90 free_jcs.f90 rhs4p_clean.f90\
														surface_move.f90
	$(FC) -c $(FLAG) velocity_initial.f90
initial.o: $(MOD3) ubc.f90 vbc.f90 initial.f90
	$(FC) -c $(FLAG) initial.f90
old_save.o: $(MOD3) old_save.f90
	$(FC) -c $(FLAG) old_save.f90
loadobject.o: $(MOD3) panelinPE.f90 loadobject.f90
	$(FC) -c $(FLAG) loadobject.f90
locateobject.o: $(MOD3) panelinPE.f90 locateobject.f90
	$(FC) -c $(FLAG) locateobject.f90
run_output.o:$(MOD3) run_output.f90 streamFunction.f90
	$(FC) -c $(FLAG) run_output.f90
draglift.o:$(MOD3) draglift.f90 interpolate.f90
	$(FC) -c $(FLAG) draglift.f90
main.o: $(MOD3) readpara.f90 startMPI.f90 Cartesiangrid.f90 exchangeobject.f90\
                locateobject.f90 loadobject.f90 allocateEuler.f90 objectInitial.f90\
                setuphypre.f90 initial.f90 run_output.f90 draglift.f90\
                cfl.f90 rk.f90 outputEuler.f90 mac_distribute.f90\
								surface_move.f90 surface_property.f90 outputObj.f90 Raycrossing.f90 main.f90
	$(FC) -c $(FLAG) main.f90
outputEuler.o: $(MOD3) outputEuler.f90
	$(FC) -c $(FLAG) outputEuler.f90
readpara.o: $(MOD2) para.mod readpara.f90
	$(FC) -c $(FLAG) readpara.f90
Divergence_reset.o: $(MOD3) Divergence_reset.f90
	$(FC) -c $(FLAG) Divergence_reset.f90
rhs4p.o: $(MOD1) rhs4p.f90 pbc.f90 Divergence_reset.f90 rhs4p_clean.f90
	$(FC) -c $(FLAG) rhs4p.f90
Raycrossing.o: $(MOD3) stretchxy.mod Raycrossing.f90
	$(FC) -c $(FLAG) Raycrossing.f90
rhs4uv.o: $(MOD1) rhs4uv.f90
	$(FC) -c $(FLAG) rhs4uv.f90
updateuv.o: $(MOD1) updateuv.f90 ubc.f90 vbc.f90
	$(FC) -c $(FLAG) updateuv.f90
Euler_reset.o: $(MOD1) Euler_reset.f90
	$(FC) -c $(FLAG) Euler_reset.f90
rk.o: $(MOD3) rk.f90 old_save.f90 uv_strain.f90 rhs4p.f90 surface_property.f90 Raycrossing.f90\
	 						allocate_jcs.f90 jc_firstsecond.f90 principle_jc_p.f90 mac_distribute.f90\
							jump_contribution.f90 solvePoisson.f90 rhs4uv.f90 updateuv.f90 free_jcs.f90\
							exchange4ghost.f90 interpolateuv.f90 Euler_reset.f90 surface_move.f90 draglift.f90
	$(FC) -c $(FLAG) rk.f90
setuphypre.o: $(MOD1) setuphypre.f90
	$(FC) -c $(FLAG) setuphypre.f90
rhs4p_clean.o: $(MOD3) rhs4p_clean.f90
	$(FC) -c $(FLAG) rhs4p_clean.f90
solvePoisson.o: $(MOD1) solvePoisson.f90
	$(FC) -c $(FLAG) solvePoisson.f90
startMPI.o: $(MOD1) indexrange.f90 coords2rank.f90 startMPI.f90
	$(FC) -c $(FLAG) startMPI.f90

# rules for others
%.o: %.f90
	$(FC) -c $(FLAG) $<

# clean all
clean:
	rm -f *.mod
	rm -f $(OBJ)
	rm *~
	rm *.exe
	clear

realclean:
	rm -f *.exe
	clear
