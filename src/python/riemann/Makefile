# List of Solvers
RIEMANN_SOLVERS = \
	rp1_acoustics \
	rp1_advection \
	rp1_burgers \
	rp1_euler_with_efix \
	rp1_nonlinear_elasticity_fwave \
	rp1_reactive_euler_with_efix \
	rp1_shallow_roe_with_efix \
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
	rp3acv
	
SHARED_OBJECTS = $(addsuffix .so,$(RIEMANN_SOLVERS))


# Generic targets
rp1_%.so: $(RIEMANN)/src/rp1_%.f90
	f2py -m $(basename $(notdir $@)) -c $^

rp2_%.so: $(RIEMANN)/src/rpn2_%.f90 $(RIEMANN)/src/rpt2_%.f90
	f2py -m $(basename $(notdir $@)) -c $^

rp3_%.so: $(RIEMANN)/src/rpn3_%.f90 $(RIEMANN)/src/rpt3_%.f90 $(RIEMANN)/src/rptt3_%.f90
	f2py -m $(basename $(notdir $@)) -c $^
	
# Phony targets
.PHONY: all clean new
all: $(SHARED_OBJECTS)

clean:
	-rm -f $(SHARED_OBJECTS)
	
new:
	$(MAKE) clean
	$(MAKE) all

# Special rules	
rp2_kpp.so: $(RIEMANN)/src/rpn2_kpp.f90 $(RIEMANN)/src/rpt2_dummy.f90
	f2py -m $(basename $(notdir $@)) -c $^

rp2_euler_mapgrid.so: $(RIEMANN)/src/rpn2_euler_mapgrid.f90 $(RIEMANN)/src/rpt2_euler_mapgrid.f90 $(RIEMANN)/src/euler_roe_solver_mapgrid.f90 $(RIEMANN)/src/getquadinfo_mapgrid.f90
	f2py -m $(basename $(notdir $@)) -c $^
	
rp2_euler_4wave.so: $(RIEMANN)/src/rpn2_euler_4wave.f90 $(RIEMANN)/src/rpt2_euler.f90
	f2py -m $(basename $(notdir $@)) -c $^
	
rp3acv.so: $(RIEMANN)/src/rpn3acv.f90 $(RIEMANN)/src/rpt3acv.f90 $(RIEMANN)/src/rptt3acv.f90
	f2py -m $(basename $(notdir $@)) -c $^
