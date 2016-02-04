#!/usr/bin/env python

from ctypes import *


class mdsys(Structure):
	_fields_ = [ ("natoms", c_int), ("nfi",c_int), ("nsteps", c_int),
				 ("dt",c_double), ("mass",c_double), ("epsilon",c_double), ("sigma",c_double), ("box",c_double), ("rcut",c_double),
				 ("ekin", c_double),("epot",c_double), ("temp", c_double),
				 ("r_x",POINTER(c_double)),("r_y", POINTER(c_double)),("r_z", POINTER(c_double)),
				 ("v_x", POINTER(c_double)),("v_y", POINTER(c_double)),("v_z", POINTER(c_double)),
			     ("f_x", POINTER(c_double)),("f_y", POINTER(c_double)),("f_z", POINTER(c_double))]

# import DSO

dso = CDLL("../lib/ljmd-serial.so")

arrayType = c_double*2

r_x = arrayType()
r_y = arrayType()
r_z = arrayType()	


r_x[0] = c_double(6.67103294321331)
r_y[0] = c_double(-10.6146871435653)
r_z[0] = c_double(12.6336939877734)
r_x[1] = c_double(1.06574058650169)
r_y[1] = c_double(-3.33432278188177)
r_z[1] = c_double(-2.59038677851747)

v_x = arrayType()
v_y = arrayType()
v_z = arrayType()

   
v_x[0] = c_double(-1.5643224621482283e-03)
v_y[0] = c_double(4.8497508563925346e-04)
v_z[0] = c_double(-4.3352481732883966e-04)
v_x[1] = c_double(4.1676710257651452e-04)
v_y[1] = c_double(2.2858522230176587e-05)
v_z[1] = c_double(-6.1985040462745732e-04)

f_x = arrayType()
f_y = arrayType()
f_z = arrayType()



sys = mdsys(natoms=2, nfi=0, nsteps=1,
		dt=0, mass=39.948, epsilon=0.2379, sigma=3.405, box=17.1580, rcut=30,
		ekin =0, epot=0, temp=0,
		r_x=r_x, r_y=r_y, r_z=r_z,
		v_x=v_x, v_y=v_y, v_z=v_z,
		f_x=f_x, f_y=f_y, f_z=f_z
		)
dso.azzero(sys.f_x, sys.natoms)
dso.azzero(sys.f_y, sys.natoms)
dso.azzero(sys.f_z, sys.natoms)
print "Start Testing"
print "********************************************"
print "Calling Forces"
dso.force(byref(sys))

assert(abs(sys.f_x[0]+0.000822) <= 1e-6)
assert(abs(sys.f_y[0]-0.001067) <= 1e-6)
assert(abs(sys.f_z[0]-0.000284) <= 1e-6)
assert(abs(sys.f_x[1]-0.000822) <= 1e-6)
assert(abs(sys.f_y[1]+0.001067) <= 1e-6)
assert(abs(sys.f_z[1]+0.000284) <= 1e-6)

print("Force test passed")
print "********************************************"
print "Calling ekin"

dso.ekin(byref(sys))

assert(abs(sys.ekin-0.163682) <= 1e-6)
assert(abs(sys.temp-54.911863) <= 1e-6)
print "Ekin test passed"
print "********************************************"
print "Calling velverlet_1"

dso.velverlet_1(byref(sys))

assert(abs(sys.r_x[0]-6.671033)<=1e-6);
assert(abs(sys.r_y[0]+10.614687)<=1e-6);
assert(abs(sys.r_z[0]-12.633694)<=1e-6);
assert(abs(sys.v_x[0]+0.001564)<=1e-6);
assert(abs(sys.v_y[0]-0.000485)<=1e-6);
assert(abs(sys.v_z[0]+0.000434)<=1e-6);
assert(abs(sys.r_x[1]-1.065741)<=1e-6);
assert(abs(sys.r_y[1]+3.334323)<=1e-6);
assert(abs(sys.r_z[1]+2.590387)<=1e-6);
assert(abs(sys.v_x[1]-0.000417)<=1e-6);
assert(abs(sys.v_y[1]-0.000023)<=1e-6);
assert(abs(sys.v_z[1]+0.000620)<=1e-6);

print "Velverlet_1 test passed"
print "********************************************"
print "Done"

