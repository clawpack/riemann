all:
	make rp1_acoustics.so
	make rp1_advection.so
	make rp1_burgers.so 
	make rp1_euler_with_efix.so
	make rp1_nonlinear_elasticity_fwave.so 
	make rp1_reactive_euler_with_efix.so 
	make rp1_shallow_roe_with_efix.so 
	make rp2_acoustics.so 
	make rp2_advection.so 
	make rp2_euler_5wave.so 
	make rp2_euler_mapgrid.so 
	make rp2_kpp.so 
	make rp2_psystem.so 
	make rp2_shallow_roe_with_efix.so 
	make rp2_shallow_sphere.so 
	make rp2_vc_acoustics.so 
	make rp2_vc_advection.so 
	make rp3_vc_acoustics.so 
	make rp3acv.so

rp1_acoustics.so: $(RIEMANN)/src/rp1_acoustics.f90
	f2py -m rp1_acoustics -c  $(RIEMANN)/src/rp1_acoustics.f90

rp1_advection.so: $(RIEMANN)/src/rp1_advection.f90
	f2py -m rp1_advection -c  $(RIEMANN)/src/rp1_advection.f90

rp1_burgers.so: $(RIEMANN)/src/rp1_burgers.f90
	f2py -m rp1_burgers -c  $(RIEMANN)/src/rp1_burgers.f90

rp1_euler_with_efix.so: $(RIEMANN)/src/rp1_euler_with_efix.f90
	f2py -m rp1_euler_with_efix -c  $(RIEMANN)/src/rp1_euler_with_efix.f90

rp1_nonlinear_elasticity_fwave.so: $(RIEMANN)/src/rp1_nonlinear_elasticity_fwave.f90
	f2py -m rp1_nonlinear_elasticity_fwave -c  $(RIEMANN)/src/rp1_nonlinear_elasticity_fwave.f90

rp1_reactive_euler_with_efix.so: $(RIEMANN)/src/rp1_reactive_euler_with_efix.f90
	f2py -m rp1_reactive_euler_with_efix -c  $(RIEMANN)/src/rp1_reactive_euler_with_efix.f90

rp1_shallow_roe_with_efix.so: $(RIEMANN)/src/rp1_shallow_roe_with_efix.f90
	f2py -m rp1_shallow_roe_with_efix -c  $(RIEMANN)/src/rp1_shallow_roe_with_efix.f90

rp2_acoustics.so: $(RIEMANN)/src/rpn2_acoustics.f90 $(RIEMANN)/src/rpt2_acoustics.f90
	f2py -m rp2_acoustics -c  $(RIEMANN)/src/rpn2_acoustics.f90 $(RIEMANN)/src/rpt2_acoustics.f90

rp2_advection.so: $(RIEMANN)/src/rpn2_advection.f90 $(RIEMANN)/src/rpt2_advection.f90
	f2py -m rp2_advection -c  $(RIEMANN)/src/rpn2_advection.f90 $(RIEMANN)/src/rpt2_advection.f90

rp2_euler_5wave.so: $(RIEMANN)/src/rpn2_euler_5wave.f90 $(RIEMANN)/src/rpt2_euler_5wave.f90
	f2py -m rp2_euler_5wave -c  $(RIEMANN)/src/rpn2_euler_5wave.f90 $(RIEMANN)/src/rpt2_euler_5wave.f90

rp2_euler_mapgrid.so: $(RIEMANN)/src/rpn2_euler_mapgrid.f90 $(RIEMANN)/src/rpt2_euler_mapgrid.f90 $(RIEMANN)/src/euler_roe_solver_mapgrid.f90 $(RIEMANN)/src/getquadinfo_mapgrid.f90
	f2py -m rp2_euler_mapgrid -c  $(RIEMANN)/src/rpn2_euler_mapgrid.f90 $(RIEMANN)/src/rpt2_euler_mapgrid.f90 $(RIEMANN)/src/euler_roe_solver_mapgrid.f90 $(RIEMANN)/src/getquadinfo_mapgrid.f90

rp2_kpp.so: $(RIEMANN)/src/rpn2_kpp.f90 $(RIEMANN)/src/rpt2_dummy.f90
	f2py -m rp2_kpp -c  $(RIEMANN)/src/rpn2_kpp.f90 $(RIEMANN)/src/rpt2_dummy.f90

rp2_psystem.so: $(RIEMANN)/src/rpn2_psystem.f90 $(RIEMANN)/src/rpt2_psystem.f90
	f2py -m rp2_psystem -c  $(RIEMANN)/src/rpn2_psystem.f90 $(RIEMANN)/src/rpt2_psystem.f90

rp2_shallow_roe_with_efix.so: $(RIEMANN)/src/rpn2_shallow_roe_with_efix.f90 $(RIEMANN)/src/rpt2_shallow_roe_with_efix.f90
	f2py -m rp2_shallow_roe_with_efix -c  $(RIEMANN)/src/rpn2_shallow_roe_with_efix.f90 $(RIEMANN)/src/rpt2_shallow_roe_with_efix.f90

rp2_shallow_sphere.so: $(RIEMANN)/src/rpn2_shallow_sphere.f90 $(RIEMANN)/src/rpt2_shallow_sphere.f90
	f2py -m rp2_shallow_sphere -c  $(RIEMANN)/src/rpn2_shallow_sphere.f90 $(RIEMANN)/src/rpt2_shallow_sphere.f90

rp2_vc_acoustics.so: $(RIEMANN)/src/rpn2_vc_acoustics.f90 $(RIEMANN)/src/rpt2_vc_acoustics.f90
	f2py -m rp2_vc_acoustics -c  $(RIEMANN)/src/rpn2_vc_acoustics.f90 $(RIEMANN)/src/rpt2_vc_acoustics.f90

rp2_vc_advection.so: $(RIEMANN)/src/rpn2_vc_advection.f90 $(RIEMANN)/src/rpt2_vc_advection.f90
	f2py -m rp2_vc_advection -c  $(RIEMANN)/src/rpn2_vc_advection.f90 $(RIEMANN)/src/rpt2_vc_advection.f90

rp3_vc_acoustics.so: $(RIEMANN)/src/rpn3_vc_acoustics.f90 $(RIEMANN)/src/rpt3_vc_acoustics.f90  $(RIEMANN)/src/rptt3_vc_acoustics.f90 
	f2py -m rp3_vc_acoustics -c  $(RIEMANN)/src/rpn3_vc_acoustics.f90 $(RIEMANN)/src/rpt3_vc_acoustics.f90 $(RIEMANN)/src/rptt3_vc_acoustics.f90

rp3acv.so: $(RIEMANN)/src/rpn3acv.f90 $(RIEMANN)/src/rpt3acv.f90 $(RIEMANN)/src/rptt3acv.f90
	f2py -m rp3_vc_acoustics -c $(RIEMANN)/src/rpn3acv.f90 $(RIEMANN)/src/rpt3acv.f90 $(RIEMANN)/src/rptt3acv.f90


