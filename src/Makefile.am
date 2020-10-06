#AM_FCFLAGS=-I$(top_srcdir)/src

# Some Fortran compilers can't build shared libraries, so sometimes we
# want to build a shared C library and a static Fortran library.  If so,
# pass the -static flag to the library linker.
#if FORTRAN_SHARED_CONDITIONAL
#else
#   AM_LDFLAGS+=-static
#endif

# http://stackoverflow.com/questions/20230827/how-to-set-include-paths-with-autotools


bin_PROGRAMS = irt3d irtpl ierf3d test_grib_hybrid 

irt3d_SOURCES = irt3d.f90 syms.f90 string_utility.f90 cubicweight.f90 io_types.f90 time.f90 sort.f90 \
		math.f90 pchip.f90 tools_grib.f90 tools_atmo.f90 tools_netcdf.f90 tools_globe.f90 \
		tools_topo.f90 eom.f90

irtpl_SOURCES = irtpl.f90 syms.f90 string_utility.f90 cubicweight.f90 time.f90 sort.f90 math.f90 pchip.f90 \
		tools_grib.f90 tools_atmo.f90 tools_netcdf.f90 tools_globe.f90 tools_topo.f90 \
		irtpl_read.f90 export.f90 irtpl_plot.f90

ierf3d_SOURCES = ierf3d.f90 syms.f90 string_utility.f90 cubicweight.f90 time.f90 sort.f90 math.f90 pchip.f90 \
		tools_grib.f90 tools_atmo.f90 tools_netcdf.f90 tools_globe.f90 tools_topo.f90 eom.f90 \
		ierf3d_netcdf.f90

test_grib_hybrid_SOURCES = test_grib_hybrid.f90 io_types.f90 syms.f90 grib_hybrid.f90 time.f90

# Hardcode the dependencies of these files. There isn't a known way of
# determining this automagically (like we do with the C files). So, when
# doing a parallel make, some modules could be made way before the
# modules they depend upon are actually made. *sigh*

io_types.o : io_types.f90
string_utility.o : string_utility.f90
#tools_complex.o : tools_complex.f90
syms.o : syms.f90 io_types.o
#byteswap.o : byteswap.f90 io_types.o
time.o : time.f90 io_types.o
sort.o : sort.f90 io_types.o
file_utility.o : file_utility.f90 io_types.o
math.o : math.f90 io_types.o syms.o

irt3d_LDFLAGS = -lnetcdf -lnetcdff

pchip.o : pchip.f90
cubicweight.o : cubicweight.f90
tools_grib.o : tools_grib.f90 time.o
tools_atmo.o : tools_atmo.f90 tools_grib.o cubicweight.o math.o syms.o
tools_netcdf.o : tools_netcdf.f90 time.o
tools_globe.o : tools_globe.f90 tools_netcdf.o
tools_topo.o : tools_topo.f90 tools_globe.o tools_netcdf.o cubicweight.o math.o syms.o
tools_netcdf.o : tools_netcdf.f90
absorption.o : absorption.f90 syms.o
eom.o : eom.f90 syms.o math.o tools_atmo.o
irt3d.o : irt3d.f90 io_types.o syms.o string_utility.o tools_atmo.o tools_topo.o tools_netcdf.o \
	cubicweight.o math.o eom.o

irtpl_LDFLAGS = -lnetcdf -lnetcdff

export.o : export.f90 syms.o math.o tools_atmo.o tools_topo.o irtpl_read.o sort.o
irtpl_read.o : irtpl_read.f90 math.o
irtpl_plot.o : irtpl_plot.f90 syms.o tools_atmo.o irtpl_read.o
irtpl.o : irtpl.f90 syms.o string_utility.o tools_atmo.o tools_topo.o tools_netcdf.o \
	irtpl_read.o irtpl_plot.o export.o math.o cubicweight.o

ierf3d_LDFLAGS = -lnetcdf -lnetcdff

ierf3d_netcdf.o : ierf3d_netcdf.f90
ierf3d.o : ierf3d.f90 io_types.o syms.o string_utility.o time.o tools_atmo.o tools_topo.o tools_netcdf.o \
	cubicweight.o math.o eom.o ierf3d_netcdf.o

test_grib_hybrid_LDFLAGS = 

grib_hybrid.o : grib_hybrid.f90 syms.o io_types.o time.o
test_grib_hybrid.o : test_grib_hybrid.f90 io_types.o grib_hybrid.o


CLEANFILES =	*.mod
EXTRA_DIST =	irt_argc.f90 irt_check.f90 irt_grib.f90 irt_marker.f90 irt_nof.f90 irt_params.f90 \
		irt_filename.f90 irt_ranges.f90 ierf3d_defaults.f90 ierf3d_parse_command_args.f90 \
		ierf3d_parse_file_args.f90 ierf3d_check_args.f90 ierf3d_grib.f90 ierf_ranges.f90 \
		ierf3d_help.f90 irt3d_raytrace.f90 ierf3d_raytrace.f90 irt_toolkit_version.f90
