#!/usr/bin/env python

import os.path
import sys
from ctypes import *
from read_write import read_input, write_output


BLEN = 200
kboltz = 0.0019872067    # boltzman constant in kcal/mol/K
mvsq2e = 2390.05736153349  # m*v^2 in kcal/mol

class mdsys(Structure):
	_fields_ = [ ("natoms", c_int), ("nfi",c_int), ("nsteps", c_int),
				 ("dt",c_double), ("mass",c_double), ("epsilon",c_double), ("sigma",c_double), ("box",c_double), ("rcut",c_double),
				 ("ekin", c_double),("epot",c_double), ("temp", c_double),
				 ("r_x",POINTER(c_double)),("r_y", POINTER(c_double)),("r_z", POINTER(c_double)),
				 ("v_x", POINTER(c_double)),("v_y", POINTER(c_double)),("v_z", POINTER(c_double)),
			     ("f_x", POINTER(c_double)),("f_y", POINTER(c_double)),("f_z", POINTER(c_double))]

# import DSO

if len(sys.argv) < 2:
	print "Expected input file!"
	exit(1)
input_file = sys.argv[1]
if not os.path.exists(input_file):
    print "Input file '" + input_file + "' does not seem to exist!"
    exit(1)

dynamic_library = "../lib/ljmd-serial.so"
if not os.path.exists(dynamic_library):
    print "Dynamic library not found at '"+dynamic_library+"' "
    exit(1)

dso = CDLL(dynamic_library)

params = read_input(input_file)

natom = int(params[0])
mass = c_double(float(params[1]))
epsilon = c_double(float(params[2]))
sigma = c_double(float(params[3]))
rcut = c_double(float(params[4]))
box = c_double(float(params[5]))
restart_f = params[6]
trajectory_f = params[7]
energies_f = params[8]
n_steps = int(params[9])
MD_time_step = float(params[10])
output_print_freq = int(params[11])

arrayType = c_double*natom

r_x = arrayType()
r_y = arrayType()
r_z = arrayType()

v_x = arrayType()
v_y = arrayType()
v_z = arrayType()

f_x = arrayType()
f_y = arrayType()
f_z = arrayType()

sys = mdsys(natoms=natom, nfi=0, nsteps=n_steps,
		dt=MD_time_step, mass=mass, epsilon=epsilon, sigma=sigma, box=box, rcut=rcut,
		ekin =0, epot=0, temp=0,
		r_x=r_x, r_y=r_y, r_z=r_z,
		v_x=v_x, v_y=v_y, v_z=v_z,
		f_x=f_x, f_y=f_y, f_z=f_z
		)

if not os.path.exists(restart_f):
    print "Input file '" + restart_f  + "' does not seem to exist!"
    exit(1)
res_file = open(restart_f, "rw+")
line = " "
line = " "
i = 0
while i < natom:
	 line = res_file.readline()
	 r_x[i] = float(line.split()[0])
	 r_y[i] = float(line.split()[1])
	 r_z[i] = float(line.split()[2])
	 i= i+1
i = 0

while i < natom:
	 line = res_file.readline()
	 v_x[i] = float(line.split()[0])
	 v_y[i] = float(line.split()[1])
	 v_z[i] = float(line.split()[2])
	 i= i+1

dso.azzero(sys.f_x, sys.natoms)
dso.azzero(sys.f_y, sys.natoms)
dso.azzero(sys.f_z, sys.natoms)

#initialize forces and energies.
sys.nfi = 0

dso.force(byref(sys))

dso.force(byref(sys))
dso.ekin(byref(sys))

erg = open(energies_f, "w")
traj = open(trajectory_f, "w")

print "Starting simulation with "+str(sys.natoms)+" atoms for "+str(sys.nsteps)+" steps."
print "     NFI            TEMP            EKIN                 EPOT              ETOT"
write_output(sys, erg, traj)

#**************************************************
#  main MD loop
#print range(0, sys.nsteps+1)
for sys.nfi in range(0, sys.nsteps+1):
	
	if not sys.nfi%output_print_freq:
		write_output(sys, erg, traj)
	dso.velverlet_1(byref(sys))
	dso.force(byref(sys))
	dso.velverlet_2(byref(sys))
	dso.ekin(byref(sys))

print "Simulation Done."
#**************************************************
