# List of Solvers
RIEMANN_SOLVERS = \
	rp1_acoustics \
	rp1_advection \
	rp1_burgers \
	rp1_euler_with_efix \
	rp1_nonlinear_elasticity_fwave \
	rp1_reactive_euler_with_efix \
	rp1_shallow_roe_with_efix \
	rp1_layered_shallow_water \
	rp1_traffic \
	rp2_acoustics \
	rp2_advection \
	rp2_euler_5wave \
	rp2_euler_4wave \
	rp2_euler_mapgrid \
	rp2_kpp \
	rp2_psystem \
	rp2_shallow_roe_with_efix \
	rp2_shallow_sphere \
	rp2_vc_acoustics \
	rp2_vc_advection \
	rp3_vc_acoustics \
	rp3_euler

SHARED_OBJECTS = $(addsuffix .so,$(RIEMANN_SOLVERS))

COMPILER=gnu95

# Generic targets
rp1_%.so: $(CLAW)/riemann/src/rp1_%.f90
	f2py --fcompiler=$(COMPILER) -m $(basename $(notdir $@)) -c $^

rp2_%.so: $(CLAW)/riemann/src/rpn2_%.f90 $(CLAW)/riemann/src/rpt2_%.f90
	f2py --fcompiler=$(COMPILER) -m $(basename $(notdir $@)) -c $^

rp3_%.so: $(CLAW)/riemann/src/rpn3_%.f90 $(CLAW)/riemann/src/rpt3_%.f90 $(CLAW)/riemann/src/rptt3_%.f90
	f2py --fcompiler=$(COMPILER) -m $(basename $(notdir $@)) -c $^
	
# Phony targets
.PHONY: all clean new
all: $(SHARED_OBJECTS)

clean:
	-rm -f $(SHARED_OBJECTS)
	
new:
	$(MAKE) clean
	$(MAKE) all

# Special rules
rp1_layered_shallow_water.so: $(CLAW)/riemann/src/rp1_layered_shallow_water.f90
	f2py --fcompiler=$(COMPILER) --link-lapack_opt -m $(basename $(notdir $@)) -c $^

rp2_kpp.so: $(CLAW)/riemann/src/rpn2_kpp.f90 $(CLAW)/riemann/src/rpt2_dummy.f90
	f2py --fcompiler=$(COMPILER) -m $(basename $(notdir $@)) -c $^

rp2_euler_mapgrid.so: $(CLAW)/riemann/src/rpn2_euler_mapgrid.f90 $(CLAW)/riemann/src/rpt2_euler_mapgrid.f90 $(CLAW)/riemann/src/euler_roe_solver_mapgrid.f90 $(CLAW)/riemann/src/getquadinfo_mapgrid.f90
	f2py --fcompiler=$(COMPILER) -m $(basename $(notdir $@)) -c $^
	
rp2_euler_4wave.so: $(CLAW)/riemann/src/rpn2_euler_4wave.f90 $(CLAW)/riemann/src/rpt2_euler.f90
	f2py --fcompiler=$(COMPILER) -m $(basename $(notdir $@)) -c $^
	
# This Riemann solver is not compatible with PyClaw - Fix some day but will
# leave the rule here for the future
rp2_layered_shallow_water.so: $(CLAW)/riemann/src/rpn2_layered_shallow_water.f90 $(CLAW)/riemann/src/rpt2_layered_shallow_water.f90 $(CLAW)/riemann/src/geoclaw_riemann_utils.f
	f2py --fcompiler=$(COMPILER) -m $(basename $(notdir $@)) -c $^
