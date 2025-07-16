# Specify a particular compiler with "make COMPILER=pgi", etc.
# Specify debugging flags with "make MODE=debug"
# If floating-point consistency is required then set the flag "CONS=1"
#   - This prevents the code from using extended 80-bit precision which is
#     incapable of generating exactly reproducible results in parallel
# Alternatively, these options can be specified as environment variables
# eg. "export COMPILER=gfortran" can be added to $HOME/.bashrc


# Compiler specific flags

# Note: you MUST specify a COMPILER option. None are specified by default.

ifeq ($(strip $(COMPILER)),)
  MAKECMDGOALS = error
  CMP_GUESS = $(shell ../scripts/infer_compiler.sh)
error:
	@echo '*** ERROR ***'
	@echo 'You MUST set a value for the COMPILER variable'
	@printf " *?* It looks like you're using "
	@printf $(CMP_GUESS)
	@echo " *?*"
	@printf ' so try "make COMPILER='
	@printf $(CMP_GUESS)
	@echo '"'
	@printf 'Alternatively, you can add "export COMPILER='
	@printf $(CMP_GUESS)
	@echo '" to $$HOME/.bashrc'
	@exit 1
endif

D = -D

ifeq ($(strip $(COMPILER)), auto)
  override COMPILER := $(shell ../scripts/infer_compiler.sh)
  NEWMAKEFLAGS := "$(filter-out COMPILER=auto, $(MAKEFLAGS))"
  $(warning Inferred compiler as $(COMPILER))
endif

# PGI
# ===
ifeq ($(strip $(COMPILER)),pgi)
  FFLAGS = -r8 -fast -fastsse -O3 -Mipa=fast,inline -Minfo # Optimised
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -Mbounds -g                                     # Debug
  endif
  MODULEFLAG = -module $(OBJDIR)
endif

# Intel
# =====
ifeq ($(strip $(COMPILER)),intel)
  FFLAGS = -O3 -g -stand f03
  ifeq ($(strip $(CONS)),1)
    FLTCONS = -pc64 -fltconsistency
  endif
  #FFLAGS = -O3 -heap-arrays 64 -ipo -xHost # Optimised (B)
  #FFLAGS = -O3 -heap-arrays 64 -ipo -xAVX  # Optimised (W)
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -fpe0 -nothreads -traceback -fltconsistency \
             -C -g -heap-arrays 64 -warn all -warn errors -fpic
    ifeq ($(strip $(SYSTEM)),Darwin)
      FFLAGS += -Wl,-no_pie
    endif
  endif
  MODULEFLAG = -module $(OBJDIR)
endif

# gfortran
# ========
ifeq ($(strip $(COMPILER)),gfortran)
  FFLAGS = -O3 -g -std=f2003

  GNUVER := $(shell gfortran -dumpversion | head -1 \
      | sed 's/[^0-9\.]*\([0-9\.]\+\).*/\1/')
  GNUMAJOR := $(shell echo $(GNUVER) | cut -f1 -d\.)
  GNUMINOR := $(shell echo $(GNUVER) | cut -f2 -d\.)

  # Allow for 99 minor revisions
  GNUVER := $(shell expr 100 \* $(GNUMAJOR) + $(GNUMINOR))

  # gfortran-4.9
  GNUGE49 := $(shell expr $(GNUVER) \>= 409)
  ifeq "$(GNUGE49)" "1"
    FFLAGS += -ffpe-summary=invalid,zero,overflow
  endif

  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -g -std=f2003 -Wall -Wextra -pedantic -Wno-unused-parameter \
             -Wno-unused-dummy-argument -Wno-integer-division \
             -ffree-line-length-80 -fcheck=all -fbacktrace -Werror \
             -ffpe-trap=invalid,zero,overflow -fbounds-check
    #FFLAGS += -ffpe-trap=underflow,denormal

    # gfortran-4.3
    GNUGE43 := $(shell expr $(GNUVER) \>= 403)
    ifeq "$(GNUGE43)" "1"
      FFLAGS += -fbacktrace -fdump-core

      # gfortran-4.6
      GNUGE46 := $(shell expr $(GNUVER) \>= 406)
      ifeq "$(GNUGE46)" "1"
        FFLAGS += -Wno-unused-dummy-argument

        # gfortran-4.8
        GNUGE48 := $(shell expr $(GNUVER) \>= 408)
        ifeq "$(GNUGE48)" "1"
          FFLAGS += -Wno-target-lifetime -fno-realloc-lhs
          FFLAGS += -Wno-maybe-uninitialized
        endif
      endif
    endif
  endif
  MODULEFLAG = -I$(OBJDIR) -J$(OBJDIR)
  ifneq ($(wildcard /usr/include/.),)
    MODULEFLAG += -I/usr/include
  endif
  INFO_FLAGS = -Wno-conversion -fno-range-check
endif

# g95
# ========
ifeq ($(strip $(COMPILER)),g95)
  FFLAGS = -O3
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -g                                        # Debug
  endif
  MODULEFLAG = -fmod=$(OBJDIR)
endif

# IBM Bluegene
# ============
ifeq ($(strip $(COMPILER)),ibm)
  FFLAGS = -O5 -qhot -qipa # Optimised
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -C -g -qfullpath -qinfo #-qkeepparm -qflttrap \
          -qnosmp -qxflag=dvz -Q! -qnounwind -qnounroll # Debug
    #FFLAGS = -O0 -qarch=qp -qtune=qp
    #FFLAGS = -qthreaded -qsmp=noauto -qsmp=omp # Hybrid stuff
  endif
  MODULEFLAG = -I$(OBJDIR) -qmoddir=$(OBJDIR)
  MPIF90 ?= mpixlf90_r

  # IBM compiler needs a -WF to recognise preprocessor directives
  D = -WF,-D
endif

# ARCHER (also works for HECToR)
# ========
ifeq ($(strip $(COMPILER)),archer)
  FFLAGS = -O3 -hipa2
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -g -ea -ec -eC -eD -eI -en -hfp_trap -Ktrap=fp -m0 -M1438,7413
  endif
  MODULEFLAG = -em -I/usr/include -I$(OBJDIR) -J$(OBJDIR)
  MPIF90 ?= ftn
endif

# Don't compile encoded source if MODE=debug or ENC=no
# Do compile encoded source if MODE=debug and ENC=yes
ifeq ($(strip $(ENC)),no)
  ENCODED_SOURCE = epoch_source_info_dummy.o
else
  ifneq ($(strip $(ENC)),yes)
    ifeq ($(strip $(MODE)),debug)
      ENCODED_SOURCE = epoch_source_info_dummy.o
    endif
  endif
endif

MPIF90 ?= mpif90
FFLAGS += -I$(SDF)/include
FFLAGS += $(MODULEFLAG)
LDFLAGS = $(FFLAGS) -L$(SDF)/lib -lsdf

# Set some of the build parameters
TARGET = epoch2d

# Set pre-processor defines
DEFINES := $(DEFINE)

# The following are a list of pre-processor defines which can be added to
# the above line modifying the code behaviour at compile time.

# Set every pseudoparticle in a species to represent the same number of
# real particles.
#DEFINES += $(D)PER_SPECIES_WEIGHT

# Enable support for zero-current particles (this will make the particle push
# slightly slower).
# Zero-current particles which behave like normal particles but do not affect
# the field, current or each other.
# WARNING: These particles may not behave in the way you might expect. See the
# documentation for further details.
#DEFINES += $(D)ZERO_CURRENT_PARTICLES

# Don't enable support for particle probes (in order to make the particle
# push slightly faster).
# Particle probes are diagnostics which gather information on particles which
# cross a plane in the simulation.
#DEFINES += $(D)NO_PARTICLE_PROBES

# Probes will also output the time particles cross the probe
#DEFINES += $(D)PROBE_TIME

# Use second order particle weighting (default is third order).
#DEFINES += $(D)PARTICLE_SHAPE_TOPHAT

# Use fifth order particle weighting (default is third order).
#DEFINES += $(D)PARTICLE_SHAPE_BSPLINE3

# Include a unique global particle ID. The first flag defines the ID using
# an 8-byte integer, the second uses 4-bytes.
#DEFINES += $(D)PARTICLE_ID
#DEFINES += $(D)PARTICLE_ID4

# Include QED routines
#DEFINES += $(D)PHOTONS

# Use the Trident process for pair production
#DEFINES += $(D)TRIDENT_PHOTONS

# Include bremsstrahlung routines
#DEFINES += $(D)BREMSSTRAHLUNG

# Use Intel-specific 'mm_prefetch' calls to load next particle in the list
# into cache ahead of time. Gives ~10% speed-up with Intel 12 compilers.
# WARNING: Under some circumstances this causes SEGFAULTS or incorrect results
#DEFINES += $(D)PREFETCH

# Turn on debugging.
#DEFINES += $(D)PARSER_DEBUG $(D)PARTICLE_DEBUG $(D)MPI_DEBUG $(D)SIMPLIFY_DEBUG

# Don't generate any output at all. Useful for benchmarking.
#DEFINES += $(D)NO_IO

# Specify charge and mass per particle rather than per species.
# This is a legacy option which will be removed soon.
#DEFINES += $(D)PER_PARTICLE_CHARGE_MASS

# Perform checks on evaluated deck expressions.
# This slows down the code but may be required if floating point exceptions
# are enabled.
#DEFINES += $(D)PARSER_CHECKING

# Disable the use of isatty C call to enable coloured welcome banner
#DEFINES += $(D)NO_USE_ISATTY

# Compile the delta-f version of EPOCH
#DEFINES += $(D)DELTAF_METHOD

# Debugging option for delta-f
#DEFINES += $(D)DELTAF_DEBUG

# The time integrated work exerted on each particle requires another variable
# that has to be calculated every time step. This will slow down the simulation!
#DEFINES += $(D)WORK_DONE_INTEGRATED

# Higuera-Cary push
#DEFINES += $(D)HC_PUSH

# Disable MPI3 features
#DEFINES += $(D)NO_MPI3

# Set options for debugging deck
# Sets warning messages to be errors
#DEFINES += $(D)DECK_DEBUG

# If you are having trouble compiling the encoded_source module or
# don't want embedded source code in restart dumps, uncomment the
# following line.
#ENCODED_SOURCE = epoch_source_info_dummy.o


# Automatically add -DPARSER_CHECKING for debug builds
ifeq ($(strip $(MODE)),debug)
  DEFINES += $(D)PARSER_CHECKING
  DEFINES += $(D)DECK_DEBUG
endif


# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

ENCODED_SOURCE ?= epoch_source_info.o

all: main

SDF := ../SDF/FORTRAN
SDFMOD = $(SDF)/include/sdf.mod
SRCDIR = src
OBJDIR = obj
BINDIR = bin
INCDIR = $(SRCDIR)/include
FC = $(MPIF90)
DATE := $(shell date +%s)
MACHINE := $(shell uname -n)
PREPROFLAGS = $(DEFINES) $(D)_COMMIT='"$(COMMIT)"' $(D)_DATE=$(DATE) \
  $(D)_MACHINE='"$(MACHINE)"' -I$(INCDIR)

FC_INFO := $(shell ${FC} --version 2>/dev/null \
    || ${FC} -V 2>/dev/null | grep '[a-zA-Z]' | head -n 1)

SRCFILES = antennae.f90 background_collisions.F90 balance.F90 \
  bethe_heitler.F90 boundary.F90 bremsstrahlung.F90 calc_df.F90 \
  collision_ionise.F90 collisions.F90 constants.F90 current_smooth.F90 \
  custom_deck.f90 custom_laser.f90 custom_parser.f90 deck.F90 \
  deck_antenna_block.f90 deck_boundaries_block.f90 \
  deck_bremsstrahlung_block.F90 deck_collision_block.F90 \
  deck_constant_block.f90 deck_control_block.F90 deck_dist_fn_block.f90 \
  deck_fields_block.f90 deck_injector_block.F90 deck_io_block.F90 \
  deck_io_global_block.F90 deck_laser_block.f90 deck_part_from_file_block.F90 \
  deck_particle_probe_block.F90 deck_qed_block.F90 deck_species_block.F90 \
  deck_stencil_block.f90 deck_subset_block.F90 deck_window_block.f90 \
  deltaf_loader.F90 diagnostics.F90 dist_fn.F90 epoch2d.F90 evaluate.F90 \
  evaluator_blocks.F90 fields.f90 file_injectors.F90 finish.f90 helper.F90 \
  ic_module.f90 injectors.F90 ionise.F90 iterators.F90 laser.f90 \
  mpi_routines.F90 mpi_subtype_control.f90 numerics.f90 particle_id_hash.F90 \
  particle_migration.F90 particle_pointer_advance.f90 particle_temperature.F90 \
  particles.F90 partlist.F90 photons.F90 prefetch.F90 probes.F90 \
  random_generator.f90 redblack_module.f90 setup.F90 shape_functions.F90 \
  shared_data.F90 shunt.F90 simple_io.F90 split_particle.F90 stack.f90 \
  strings.f90 strings_advanced.f90 terminal_controls.F90 timer.f90 \
  tokenizer_blocks.f90 utilities.f90 version_data.F90 welcome.F90 window.F90

OBJFILES := $(SRCFILES:.f90=.o)
OBJFILES := $(OBJFILES:.F90=.o)

INCLUDES = $(INCDIR)/particle_head.inc $(INCDIR)/particle_to_grid.inc \
  $(INCDIR)/fields_at_particle.inc $(INCDIR)/fields_at_particle_head.inc \
  $(INCDIR)/bspline3/b_part.inc $(INCDIR)/bspline3/e_part.inc \
  $(INCDIR)/bspline3/gx.inc $(INCDIR)/bspline3/gxfac.inc \
  $(INCDIR)/bspline3/hx_dcell.inc \
  $(INCDIR)/tophat/b_part.inc $(INCDIR)/tophat/e_part.inc \
  $(INCDIR)/tophat/gx.inc $(INCDIR)/tophat/gxfac.inc \
  $(INCDIR)/tophat/hx_dcell.inc \
  $(INCDIR)/triangle/b_part.inc $(INCDIR)/triangle/e_part.inc \
  $(INCDIR)/triangle/gx.inc $(INCDIR)/triangle/gxfac.inc \
  $(INCDIR)/triangle/hx_dcell.inc

OBJFILES := $(OBJFILES) $(ENCODED_SOURCE)

FULLTARGET = $(BINDIR)/$(TARGET)

SOURCE_ALL := Makefile $(SRCFILES) $(INCLUDES) pack.py pack.sh \
  $(SRCDIR)/COMMIT $(SRCDIR)/gen_commit_string.sh epoch_source_info_dummy.f90

VPATH = $(SRCDIR):$(SRCDIR)/deck:$(SRCDIR)/housekeeping:$(SRCDIR)/io:\
  $(SRCDIR)/parser:$(SRCDIR)/physics_packages:$(SRCDIR)/user_interaction:\
  $(SDF)/src:$(OBJDIR)

PACK_SDF = $(SRCDIR)/pack.sh

GIT_WORK_TREE = ".."
PACK_PREFIX = epoch
PACK_SOURCE_CODE = 1
PACK_GIT_DIFF = 1
PACK_GIT_DIFF_FROM_ORIGIN = 1
GENERATE_CHECKSUM = 1
F77_OUTPUT = 0
PACK_OPTS = $(GIT_WORK_TREE) $(PACK_PREFIX) $(PACK_SOURCE_CODE) \
    $(PACK_GIT_DIFF) $(PACK_GIT_DIFF_FROM_ORIGIN) $(GENERATE_CHECKSUM) \
    $(F77_OUTPUT)

-include $(SRCDIR)/COMMIT

ifeq ($(DONE_COMMIT),)
main: commit
else
main: $(FULLTARGET)
endif

# Rule to build the fortran files

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $(OBJDIR)/$@ $<

%.o: %.F90
	$(FC) -c $(FFLAGS) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	@mkdir -p $(BINDIR)
	$(FC) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(LDFLAGS)

$(OBJDIR)/epoch_source_info.f90: $(SOURCE_ALL)
	sh $(PACK_SDF) $(PACK_OPTS) $@ "$(FC_INFO)" "$(FFLAGS)" $^
epoch_source_info.o: $(OBJDIR)/epoch_source_info.f90 $(SOURCE_ALL)
	$(FC) -c $(FFLAGS) $(INFO_FLAGS) -o $(OBJDIR)/$@ $<

$(SDFMOD):
	$(MAKE) MAKEFLAGS=$(NEWMAKEFLAGS) COMPILER=$(COMPILER) -C $(SDF)

rm:
	@rm -rf $(BINDIR) $(OBJDIR)

clean: tidy
cleanall: tidy
	@rm -rf $(BINDIR)

tidy:
	@rm -rf $(OBJDIR) *~ *.pbs.* *.sh.* $(SRCDIR)/*~ *.log
	$(MAKE) -C $(SDF) cleanall

datatidy:
	@rm -rf Data/*

tarball:
	@cd ..; sh ./scripts/make_tarball.sh

visit:
	@cd $(SDF)/../VisIt; ./build

visitclean:
	@cd $(SDF)/../VisIt; make clean; ./build -c; \
	  rm -rf .depend *.d *Info.C *Info.h CMake* cmake* Makefile

sdfutils:
	@cd $(SDF)/../C; make
	@cd $(SDF)/../utilities; ./build.sh

sdfutilsclean:
	@cd $(SDF)/../C; make clean
	@cd $(SDF)/../utilities; ./build.sh -c

sdf_reader_example: $(BINDIR)/sdf_reader_example

$(BINDIR)/sdf_reader_example: sdf_reader_example.f90 $(SDFMOD)
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -o $@ $< $(LDFLAGS)

$(OBJFILES): | $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

commit: FORCE
	@sh $(SRCDIR)/gen_commit_string.sh && $(MAKE) $(MAKECMDGOALS) DONE_COMMIT=1

FORCE:

.PHONY: commit clean cleanall tidy datatidy visit visitclean main FORCE \
        sdf_reader_example

# All the dependencies

antennae.o: antennae.f90 evaluate.o shared_data.o shunt.o
background_collisions.o: background_collisions.F90 calc_df.o partlist.o setup.o
balance.o: balance.F90 boundary.o constants.o mpi_subtype_control.o \
  redblack_module.o timer.o utilities.o
bethe_heitler.o: bethe_heitler.F90 calc_df.o particles.o partlist.o setup.o
boundary.o: boundary.F90 injectors.o laser.o mpi_subtype_control.o \
  particle_id_hash.o particle_temperature.o partlist.o utilities.o
bremsstrahlung.o: bremsstrahlung.F90 bethe_heitler.o calc_df.o partlist.o \
  setup.o
calc_df.o: calc_df.F90 boundary.o
collision_ionise.o: collision_ionise.F90 calc_df.o collisions.o prefetch.o
collisions.o: collisions.F90 calc_df.o prefetch.o
constants.o: constants.F90 $(SDFMOD)
current_smooth.o: current_smooth.F90 boundary.o constants.o shape_functions.o
custom_deck.o: custom_deck.f90 constants.o
custom_laser.o: custom_laser.f90 shared_data.o
custom_parser.o: custom_parser.f90 shared_data.o stack.o tokenizer_blocks.o
deck.o: deck.F90 bremsstrahlung.o custom_deck.o deck_antenna_block.o \
  deck_boundaries_block.o deck_bremsstrahlung_block.o deck_collision_block.o \
  deck_constant_block.o deck_control_block.o deck_dist_fn_block.o \
  deck_fields_block.o deck_injector_block.o deck_io_block.o \
  deck_io_global_block.o deck_laser_block.o deck_part_from_file_block.o \
  deck_particle_probe_block.o deck_qed_block.o deck_species_block.o \
  deck_stencil_block.o deck_subset_block.o deck_window_block.o photons.o \
  strings.o utilities.o version_data.o $(SDFMOD)
deck_antenna_block.o: deck_antenna_block.f90 antennae.o evaluate.o shunt.o \
  strings_advanced.o
deck_boundaries_block.o: deck_boundaries_block.f90 strings_advanced.o \
  utilities.o
deck_bremsstrahlung_block.o: deck_bremsstrahlung_block.F90 strings_advanced.o \
  utilities.o
deck_collision_block.o: deck_collision_block.F90 collisions.o strings_advanced.o
deck_constant_block.o: deck_constant_block.f90 evaluate.o shunt.o
deck_control_block.o: deck_control_block.F90 fields.o strings_advanced.o \
  timer.o $(SDFMOD)
deck_dist_fn_block.o: deck_dist_fn_block.f90 dist_fn.o strings_advanced.o \
  utilities.o
deck_fields_block.o: deck_fields_block.f90 simple_io.o strings_advanced.o
deck_injector_block.o: deck_injector_block.F90 file_injectors.o injectors.o \
  shunt.o strings_advanced.o utilities.o
deck_io_block.o: deck_io_block.F90 strings_advanced.o utilities.o
deck_io_global_block.o: deck_io_global_block.F90 strings_advanced.o utilities.o
deck_laser_block.o: deck_laser_block.f90 laser.o strings_advanced.o utilities.o
deck_part_from_file_block.o: deck_part_from_file_block.F90 \
  deck_species_block.o strings_advanced.o utilities.o
deck_particle_probe_block.o: deck_particle_probe_block.F90 probes.o \
  strings_advanced.o utilities.o
deck_qed_block.o: deck_qed_block.F90 strings_advanced.o utilities.o
deck_species_block.o: deck_species_block.F90 partlist.o setup.o simple_io.o \
  strings_advanced.o utilities.o
deck_stencil_block.o: deck_stencil_block.f90 fields.o strings_advanced.o
deck_subset_block.o: deck_subset_block.F90 particle_id_hash.o shunt.o \
  strings_advanced.o utilities.o
deck_window_block.o: deck_window_block.f90 strings_advanced.o
deltaf_loader.o: deltaf_loader.F90 shared_data.o
diagnostics.o: diagnostics.F90 antennae.o calc_df.o constants.o deck.o \
  deck_io_block.o dist_fn.o evaluate.o iterators.o particle_id_hash.o probes.o \
  setup.o strings.o timer.o version_data.o window.o $(SDFMOD) $(ENCODED_SOURCE)
dist_fn.o: dist_fn.F90 mpi_subtype_control.o particles.o $(SDFMOD)
epoch2d.o: epoch2d.F90 antennae.o background_collisions.o balance.o \
  bremsstrahlung.o calc_df.o collision_ionise.o collisions.o current_smooth.o \
  deck.o diagnostics.o fields.o finish.o helper.o ic_module.o injectors.o \
  ionise.o mpi_routines.o particle_migration.o particles.o photons.o setup.o \
  split_particle.o welcome.o window.o
	$(FC) -c $(FFLAGS) $(FLTCONS) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<
evaluate.o: evaluate.F90 shunt.o stack.o
evaluator_blocks.o: evaluator_blocks.F90 custom_parser.o stack.o strings.o
fields.o: fields.f90 boundary.o
	$(FC) -c $(FFLAGS) $(FLTCONS) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<
file_injectors.o: file_injectors.F90 partlist.o utilities.o
finish.o: finish.f90 antennae.o background_collisions.o collisions.o \
  constants.o deck.o diagnostics.o dist_fn.o injectors.o ionise.o laser.o \
  probes.o setup.o window.o
helper.o: helper.F90 balance.o boundary.o deltaf_loader.o partlist.o simple_io.o
ic_module.o: ic_module.f90 helper.o shared_data.o
injectors.o: injectors.F90 evaluate.o file_injectors.o particle_temperature.o \
  partlist.o random_generator.o utilities.o
ionise.o: ionise.F90 boundary.o numerics.o partlist.o random_generator.o \
  utilities.o
iterators.o: iterators.F90 particle_id_hash.o particle_pointer_advance.o \
  partlist.o
laser.o: laser.f90 custom_laser.o evaluate.o
mpi_routines.o: mpi_routines.F90 helper.o
mpi_subtype_control.o: mpi_subtype_control.f90 shared_data.o
numerics.o: numerics.f90 constants.o
particle_id_hash.o: particle_id_hash.F90 constants.o random_generator.o \
  shared_data.o utilities.o
particle_migration.o: particle_migration.F90 calc_df.o partlist.o prefetch.o
particle_pointer_advance.o: particle_pointer_advance.f90 shared_data.o
particle_temperature.o: particle_temperature.F90 constants.o evaluate.o \
  random_generator.o
particles.o: particles.F90 boundary.o partlist.o prefetch.o
partlist.o: partlist.F90 particle_id_hash.o random_generator.o shared_data.o
photons.o: photons.F90 partlist.o
prefetch.o: prefetch.F90 shared_data.o
probes.o: probes.F90 partlist.o $(SDFMOD)
random_generator.o: random_generator.f90
redblack_module.o: redblack_module.f90 partlist.o
setup.o: setup.F90 antennae.o balance.o boundary.o constants.o fields.o \
  helper.o injectors.o laser.o mpi_routines.o mpi_subtype_control.o \
  particle_id_hash.o shunt.o split_particle.o timer.o version_data.o welcome.o \
  window.o $(SDFMOD)
shape_functions.o: shape_functions.F90 constants.o
shared_data.o: shared_data.F90 constants.o $(SDFMOD)
shunt.o: shunt.F90 evaluator_blocks.o tokenizer_blocks.o utilities.o
simple_io.o: simple_io.F90 boundary.o mpi_subtype_control.o
split_particle.o: split_particle.F90 boundary.o
stack.o: stack.f90 constants.o
strings.o: strings.f90 shared_data.o
strings_advanced.o: strings_advanced.f90 evaluate.o shunt.o
terminal_controls.o: terminal_controls.F90
timer.o: timer.f90 constants.o
tokenizer_blocks.o: tokenizer_blocks.f90 strings.o
utilities.o: utilities.f90 constants.o shared_data.o
version_data.o: version_data.F90 $(SRCDIR)/COMMIT
welcome.o: welcome.F90 shared_data.o terminal_controls.o version_data.o
window.o: window.F90 boundary.o evaluate.o partlist.o
