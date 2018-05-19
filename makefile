# For NERSC Edison and Cori, you must first load the cray-netcdf module:
#   module load cray-netcdf


ifdef NERSC_HOST
        HOSTNAME = $(NERSC_HOST)
else
        HOSTNAME="laptop"
endif

ifeq ($(HOSTNAME),edison)
	FC = ftn
	# In the next line, we prmote real to double, as gs2 does, using -r8
	EXTRA_COMPILE_FLAGS = -r8
	EXTRA_LINK_FLAGS =
else ifeq ($(HOSTNAME),cori)
	FC = ftn
	# In the next line, we prmote real to double, as gs2 does, using -r8
	EXTRA_COMPILE_FLAGS = -r8
	EXTRA_LINK_FLAGS =
else
	# Options for my macbook laptop with macports
	FC = mpif90
	# In the next line, we prmote real to double, as gs2 does, using -fdefault-real-8 -fdefault-double-8
	EXTRA_COMPILE_FLAGS = -I/opt/local/include -ffree-line-length-none -fdefault-real-8 -fdefault-double-8
	EXTRA_LINK_FLAGS =  -L/opt/local/lib -lnetcdff  -lnetcdf
endif


# The variable LIBSTELL_DIR should either be "mini_libstell", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a directory containing libstell .mod files
# elsewhere on your system.
LIBSTELL_DIR = mini_libstell
#LIBSTELL_DIR=/Users/mattland/stellopt/LIBSTELL/Release

# The variable LIBSTELL_FOR_GS2 should either be "mini_libstell/mini_libstell.a", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a libstell.a library elsewhere on your system.
LIBSTELL_FOR_GS2 = mini_libstell/mini_libstell.a
#LIBSTELL_FOR_GS2=/Users/mattland/stellopt/LIBSTELL/Release/libstell.a

# End of system-dependent variable assignments

export

.PHONY: all clean

all: test_vmec_to_gs2_geometry_interface

include makefile.depend

%.o: %.f90
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

%.o: %.f
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

test_vmec_to_gs2_geometry_interface: $(OBJ_FILES)
	$(FC) -o test_vmec_to_gs2_geometry_interface $(OBJ_FILES) $(LIBSTELL_FOR_GS2) $(EXTRA_LINK_FLAGS)

mini_libstell/mini_libstell.a:
	$(MAKE) -C mini_libstell

clean:
	rm -f *.o *.mod *.MOD *~ test_vmec_to_gs2_geometry_interface
	cd mini_libstell; rm -f *.o *.mod *.MOD *.a

test_make:
	@echo HOSTNAME is $(HOSTNAME)
	@echo FC is $(FC)
	@echo EXTRA_COMPILE_FLAGS is $(EXTRA_COMPILE_FLAGS)
	@echo EXTRA_LINK_FLAGS is $(EXTRA_LINK_FLAGS)
