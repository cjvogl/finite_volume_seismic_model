
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = amrclaw                  # Clawpack package to use
EXE = xamr                          # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot_surface.py   # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first
RESTART ?= False                    # Should = clawdata.restart in setrun

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?= -fopenmp

# ---------------------------------
# package sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/2d
include $(AMRLIB)/Makefile.amr_2d

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

MODULES = \
  fault_module.f90

SOURCES = \
  qinit.f \
  setprob.f90 \
  bc2amr.f \
  setaux.f90 \
  mapc2p.f90 \
  b4step2.f90 \
  rpn2_vc_elasticity_mapped_fault.f90 \
  rpt2_vc_elasticity_mapped_fault.f90

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


#-------------------------------------------------------------------
# Add make options for paper

FIGURE1_LIST = setrun_fault.py setplot_fault.py setfault_verydeep.py
FIGURE2_LIST = setrun_both.py setplot_both.py setfault_baseline.py
FIGURE3_LIST = setrun_surface.py setplot_surface.py setfault_baseline.py \
	setfault_deep.py setfault_steep.py setfault_thin.py
FIGURE4_LIST = setrun_surface.py setplot_surface.py setfault_smooth.py \
	setfault_bimodal.py

figures: ;
	rm -rf .output* .plots* .figure1 .figure3 .figure4
	$(MAKE) new
	$(MAKE) .figures

.figures: $(MAKEFILE_LIST) $(FIGURE1_LIST) $(FIGURE2_LIST) $(FIGURE3_LIST);
	$(MAKE) .exe
	$(MAKE) .figure1
	$(MAKE) .figure2
	$(MAKE) .figure3
	$(MAKE) .figure4

.figure1: $(FIGURE1_LIST);
	$(MAKE) .output_fault
	$(MAKE) .plots_fault
	touch .figure1

.figure2: $(FIGURE2_LIST);
	$(MAKE) .output_both
	$(MAKE) .plots_both
	touch .figure2

.figure3: $(FIGURE3_LIST);
	$(MAKE) .output_baseline
	$(MAKE) .plots_baseline
	$(MAKE) .output_deep
	$(MAKE) .plots_deep
	$(MAKE) .output_steep
	$(MAKE) .plots_steep
	$(MAKE) .output_thin
	$(MAKE) .plots_thin
	touch .figure3

.figure4: $(FIGURE4_LIST);
	$(MAKE) .output_smooth
	$(MAKE) .plots_smooth
	$(MAKE) .output_bimodal
	$(MAKE) .plots_bimodal
	touch .figure4

.output_fault: setrun_fault.py setfault_verydeep.py;
	$(CLAW_PYTHON) setfault_verydeep.py
	$(CLAW_PYTHON) setrun_fault.py $(CLAW_PKG)
	$(CLAW_PYTHON) $(CLAW)/clawutil/src/python/clawutil/runclaw.py $(EXE) \
	output_fault $(OVERWRITE) $(RESTART) . $(GIT_STATUS) $(NOHUP) \
	$(NICE)
	touch .output_fault

.output_both: setrun_both.py setfault_baseline.py;
	$(CLAW_PYTHON) setfault_baseline.py
	$(CLAW_PYTHON) setrun_both.py $(CLAW_PKG)
	$(CLAW_PYTHON) $(CLAW)/clawutil/src/python/clawutil/runclaw.py $(EXE) \
	output_both $(OVERWRITE) $(RESTART) . $(GIT_STATUS) $(NOHUP) \
	$(NICE)
	touch .output_both

.output_baseline: setrun_surface.py setfault_baseline.py;
	$(CLAW_PYTHON) setfault_baseline.py
	$(CLAW_PYTHON) setrun_surface.py $(CLAW_PKG)
	$(CLAW_PYTHON) $(CLAW)/clawutil/src/python/clawutil/runclaw.py $(EXE) \
	output_baseline $(OVERWRITE) $(RESTART) . $(GIT_STATUS) $(NOHUP) \
	$(NICE)
	touch .output_baseline

.output_deep: setrun_surface.py setfault_deep.py;
	$(CLAW_PYTHON) setfault_deep.py
	$(CLAW_PYTHON) setrun_surface.py $(CLAW_PKG)
	$(CLAW_PYTHON) $(CLAW)/clawutil/src/python/clawutil/runclaw.py $(EXE) \
	output_deep $(OVERWRITE) $(RESTART) . $(GIT_STATUS) $(NOHUP) \
	$(NICE)
	touch .output_deep

.output_steep: setrun_surface.py setfault_steep.py;
	$(CLAW_PYTHON) setfault_steep.py
	$(CLAW_PYTHON) setrun_surface.py $(CLAW_PKG)
	$(CLAW_PYTHON) $(CLAW)/clawutil/src/python/clawutil/runclaw.py $(EXE) \
	output_steep $(OVERWRITE) $(RESTART) . $(GIT_STATUS) $(NOHUP) \
	$(NICE)
	touch .output_steep

.output_thin: setrun_surface.py setfault_thin.py;
	$(CLAW_PYTHON) setfault_thin.py
	$(CLAW_PYTHON) setrun_surface.py $(CLAW_PKG)
	$(CLAW_PYTHON) $(CLAW)/clawutil/src/python/clawutil/runclaw.py $(EXE) \
	output_thin $(OVERWRITE) $(RESTART) . $(GIT_STATUS) $(NOHUP) \
	$(NICE)
	touch .output_thin

.output_smooth: setrun_surface.py setfault_smooth.py
	$(CLAW_PYTHON) setfault_smooth.py
	$(CLAW_PYTHON) setrun_surface.py $(CLAW_PKG)
	$(CLAW_PYTHON) $(CLAW)/clawutil/src/python/clawutil/runclaw.py $(EXE) \
	output_smooth $(OVERWRITE) $(RESTART) . $(GIT_STATUS) $(NOHUP) \
	$(NICE)
	touch .output_smooth

.output_bimodal: setrun_surface.py setfault_bimodal.py
	$(CLAW_PYTHON) setfault_bimodal.py
	$(CLAW_PYTHON) setrun_surface.py $(CLAW_PKG)
	$(CLAW_PYTHON) $(CLAW)/clawutil/src/python/clawutil/runclaw.py $(EXE) \
	output_bimodal $(OVERWRITE) $(RESTART) . $(GIT_STATUS) $(NOHUP) \
	$(NICE)
	touch .output_bimodal

.plots_fault: .output_fault setplot_fault.py;
	$(PLOTCMD) output_fault plots_fault setplot_fault.py
	touch .plots_fault

.plots_both: .output_both setplot_both.py;
	$(PLOTCMD) output_both plots_both setplot_both.py
	touch .plots_both

.plots_baseline: .output_baseline setplot_surface.py;
	$(PLOTCMD) output_baseline plots_baseline setplot_surface.py
	touch .plots_baseline

.plots_deep: .output_deep setplot_surface.py;
	$(PLOTCMD) output_deep plots_deep setplot_surface.py
	touch .plots_deep

.plots_steep: .output_steep setplot_surface.py;
	$(PLOTCMD) output_steep plots_steep setplot_surface.py
	touch .plots_steep

.plots_thin: .output_thin setplot_surface.py;
	$(PLOTCMD) output_thin plots_thin setplot_surface.py
	touch .plots_thin

.plots_smooth: .output_smooth setplot_surface.py;
	$(PLOTCMD) output_smooth plots_smooth setplot_surface.py
	touch .plots_smooth

.plots_bimodal: .output_bimodal setplot_surface.py;
	$(PLOTCMD) output_bimodal plots_bimodal setplot_surface.py
	touch .plots_bimodal
