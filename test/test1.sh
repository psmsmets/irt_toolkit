#!/bin/bash
# ***************************************************************************
#
#                          Infrasound Ray Tracer 3D
#                                     IRT3D
#
# ***************************************************************************
#
# Setup and execution file for the IRT3D ray tracer.
#
# Update all following variables to meet your ray trace specifications. If
# output data is stored, the plots can be (re)generated afterwards. Just
# comment the run $irt3d line.
# 
# Note :
#  
# Be aware that output type 2 and 3 (storing the rays) generate large amounts
# of data. For a first run, use type 1 (which is faster and only stores
# reflections).
#
#
# version 1.4 (February 2017)
#
# ---
# Pieter Smets (pieter.smets@knmi.nl)


# --------------------------------------------
# 
#               FILE LOCATIONS
# 
# --------------------------------------------

prefix="test1b"
tmpfile="$prefix.irt3d.tmp"
export GMT_CEFF='../data/tempc.cpt'
export GMT_WAXLER='../data/waxler.cpt'

# --------------------------------------------
# 
#              GENERAL OPTIONS
# 
# --------------------------------------------

# ECMWF atmospheric profile parameters
echo "test.grib" > $tmpfile		# atmospheric specifications file
echo "0 70 .5 0" >> $tmpfile		# altitude: min, max, and stepsize (km)
echo "no" >> $tmpfile			# reverse horizontal wind (no/yes)

# Topography
echo "local_etopo1.nc" >> $tmpfile		# topography file (no/yes/'location')

# Source location
echo "50.644 5.576 0.0" >> $tmpfile	# source: latitude (deg), longitude (deg), altitude (km)

# Ray trace options
echo "19.5 20.5 .5" >> $tmpfile		# elevation range (min, max and stepsize)
                                        # 0º, 90º (zenith)
echo "343.42796" >> $tmpfile		# azimuth range (min, max and stepsize)
                                        # 0º N, 90ºE, 180ºS, 270ºW
echo "1.0" >> $tmpfile			# transmission loss frequency (Hz)
echo "m" >> $tmpfile			# transmission loss RE unit (m/km)
echo "1.0" >> $tmpfile			# rk4 stepsize (km^2/s)

# Output
echo "$prefix" >> $tmpfile		# data prefix output name
echo "3" >> $tmpfile			# ray data output type (1=reflections,2=rays,3=both)
echo "350.42796 5.0 0.5 2 F" >> $tmpfile	# vertical cross-section: bearing (deg), horizontal and vertical resolution (km),
					# number of rays to skip in plot (-)
					# type 'no' for no output
echo "yes" >> $tmpfile			# verbose (no/yes)

# Stop conditions
echo "d 500" >> $tmpfile		
					# Stop conditions for ray tracing. First argument determines the type stop
					# condition, followed by values. If left blank, default is ll for segmented profiles
					# and d 1000.0 km for continuous profiles. Possible options:
					#	a) ll latmin latmax lonmin lonmax
					# 	b) ld latmin latmax dmax
					# 	c) dl dmax latmin latmax
					# 	d) d dmax


# --------------------------------------------
# 
#               PLOT OPTIONS
# 
# --------------------------------------------

# Additional markers to plot		# latitude (deg), longitude (deg), gmt symbol (a,c,d,g,h,i,n,s,t)
					# --> list can be as long as you want
echo "052.099184 05.176384 t" >> $tmpfile	# DBN
echo "053.117180 04.872100 t" >> $tmpfile	# TEX
echo "052.060724 05.874939 t" >> $tmpfile	# DIA
echo "051.970275 04.926201 t" >> $tmpfile	# CIA
echo "052.908682 06.868268 t" >> $tmpfile	# LOFAR
echo "-54.58057 -067.30923 d" >> $tmpfile	# I02AR
echo "-34.59761 0116.35669 d" >> $tmpfile	# I04AU
echo "-42.49080 0147.68063 d" >> $tmpfile	# I05AU
echo "-12.14645 0096.82032 d" >> $tmpfile	# I06AU
echo "-19.93482 0134.32953 d" >> $tmpfile	# I07AU
echo "-16.21523 -068.45345 d" >> $tmpfile	# I08BO
echo "-15.63797 -048.01642 d" >> $tmpfile	# I09BR
echo "050.20147 -096.02685 d" >> $tmpfile	# I10CA
echo "015.25729 -023.18388 d" >> $tmpfile	# I11CV
echo "-27.12726 -109.36265 d" >> $tmpfile	# I13CL
echo "-33.65379 -078.79598 d" >> $tmpfile	# I14CL
echo "006.67036 -004.85691 d" >> $tmpfile	# I17CI
echo "077.47556 -069.28776 d" >> $tmpfile	# I18DK
echo "011.47401 0043.17308 d" >> $tmpfile	# I19DJ
echo "-08.86783 -140.15907 d" >> $tmpfile	# I21FR
echo "-22.18445 0166.84592 d" >> $tmpfile	# I22FR
echo "-49.34578 0070.24159 d" >> $tmpfile	# I23FR
echo "-17.74929 -149.29582 d" >> $tmpfile	# I24FR
echo "048.85162 0013.71313 d" >> $tmpfile	# I26DE
echo "-70.66197 -008.32127 d" >> $tmpfile	# I27DE
echo "035.30776 0140.31376 d" >> $tmpfile	# I30JP
echo "050.40697 0058.03482 d" >> $tmpfile	# I31KZ
echo "-01.24216 0036.82721 d" >> $tmpfile	# I32KE
echo "-19.01086 0047.30502 d" >> $tmpfile	# I33MG
echo "047.80172 0106.41012 d" >> $tmpfile	# I34MN
echo "-19.19135 0017.57678 d" >> $tmpfile	# I35NA
echo "-43.91662 -176.48337 d" >> $tmpfile	# I36NZ
echo "069.07408 0018.60763 d" >> $tmpfile	# I37NO
echo "007.53547 0134.54704 d" >> $tmpfile	# I39PW
echo "-26.34230 -057.31188 d" >> $tmpfile	# I41PY
echo "039.04233 -028.00554 d" >> $tmpfile	# I42PT
echo "056.72136 0037.21759 d" >> $tmpfile	# I43RU
echo "053.10580 0157.71390 d" >> $tmpfile	# I44RU
echo "044.19990 0131.97730 d" >> $tmpfile	# I45RU
echo "053.94872 0084.81891 d" >> $tmpfile	# I46RU
echo "-28.62112 0025.23523 d" >> $tmpfile	# I47ZA
echo "035.80523 0009.32202 d" >> $tmpfile	# I48TN
echo "-37.08995 -012.33192 d" >> $tmpfile	# I49GB
echo "-07.93774 -014.37517 d" >> $tmpfile	# I50GB
echo "032.36154 -064.69874 d" >> $tmpfile	# I51GB
echo "-07.37781 0072.48416 d" >> $tmpfile	# I52GB
echo "064.87500 -147.86114 d" >> $tmpfile	# I53US
echo "-77.73149 0167.58742 d" >> $tmpfile	# I55US
echo "048.26408 -117.12567 d" >> $tmpfile	# I56US
echo "033.60585 -116.45328 d" >> $tmpfile	# I57US
echo "019.59153 -155.89360 d" >> $tmpfile	# I59US
#echo "064.77145 -146.88665 d" >> $tmpfile	# ILAR
#echo "068.30652 -133.52543 d" >> $tmpfile	# INK
#echo "050.20100 -096.02700 d" >> $tmpfile	# ISM
#echo "054.13160 -165.97920 t" >> $tmpfile	# Akutan, Alaska
#echo "058.77050 -153.67090 t" >> $tmpfile	# Fourpeaked, Alaska
#echo "053.46840 -168.17390 t" >> $tmpfile	# Okmok, Alaska

# --------------------------------------------
# 
#                     RUN
# 
# --------------------------------------------

../src/./irt3d $tmpfile
../src/./irtpl $tmpfile

# --------------------------------------------
# 
#                     CLEAN
# 
# --------------------------------------------

rm $tmpfile $prefix*.xy $prefix*.xyz $prefix*.gmt $prefix*.out

# --------------------------------------------
