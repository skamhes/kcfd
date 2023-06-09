##########################################################
# Makefile for EDU2D-CCFV-Euler-EXPLCT
##########################################################
 PROGRAM = kcfd
##########################################################
# Suffix Rule for f90
# Note: use "gfortran -O2" for best performance, but
#       don't use it until you're sure bugs are removed.
##########################################################
.SUFFIXES : .o .f90
.f90.o:
#	gfortran -O0 -g -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace -fall-intrinsics -c $<
	gfortran -O2 -pg -c $<
#	gfortran -O3 -c $<
##########################################################
SDIR = .

OBCTS = $(SDIR)/kcfd_module_input_parameter.o\
	$(SDIR)/kcfd_module_common_data.o\
	$(SDIR)/kcfd_module_ccfv_data_grid.o\
	$(SDIR)/kcfd_module_ccfv_data_soln.o\
	$(SDIR)/kcfd_module_ccfv_gradient.o\
	$(SDIR)/kcfd_module_write_files.o\
	$(SDIR)/kcfd_module_flux.o\
	$(SDIR)/kcfd_derivative_data_df5.o\
	$(SDIR)/kcfd_flux_functions_ddt.o\
	$(SDIR)/kcfd_module_flux_jac_interface.o\
	$(SDIR)/kcfd_module_bc_states.o\
	$(SDIR)/kcfd_module_ccfv_limiter.o\
	$(SDIR)/kcfd_module_ccfv_residual.o\
	$(SDIR)/kcfd_module_jacobian.o\
	$(SDIR)/kcfd_module_linear_solver.o\
	$(SDIR)/kcfd_module_steady_solver.o\
	$(SDIR)/kcfd_main.o
##########################################################
# Make executable "mg" 
# Note: use "gfortran -O2" for best performance, but
#       don't use it until you're sure bugs are removed.
##########################################################
$(PROGRAM): $(OBCTS)
#	gfortran -O0 -g -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace -fall-intrinsics -o $@ $(OBCTS)
	gfortran -O3 -o $@ $(OBCTS)
#	gfortran -O2 -pg $@ $(OBCTS)
##########################################################
# Clean up
##########################################################
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.mod0
	rm -f kcfd
