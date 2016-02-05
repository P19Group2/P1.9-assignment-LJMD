#!/usr/bin/python

def read_input(input_file):
    input_param = open(input_file, "rw+")
    n_atom = input_param.readline().split()[0]
    mass = input_param.readline().split()[0]
    epsilon = input_param.readline().split()[0]
    sigma = input_param.readline().split()[0]
    rcut = input_param.readline().split()[0]
    box = input_param.readline().split()[0]
    restart_f = input_param.readline().split()[0]
    trajectory_f = input_param.readline().split()[0]
    energies_f = input_param.readline().split()[0]
    n_steps = input_param.readline().split()[0]
    MD_time_step = input_param.readline().split()[0]
    output_print_freq = input_param.readline().split()[0]

    return n_atom, mass, epsilon,sigma, rcut, box,restart_f, trajectory_f,energies_f, n_steps, MD_time_step, output_print_freq


def write_output(sys, erg, traj):

    print str(sys.nfi)+" "+str(sys.temp)+" "+str(sys.ekin)+" "+str(sys.epot)+" "+str(sys.ekin+sys.epot)
    erg.write(str([sys.nfi, sys.temp, sys.ekin, sys.epot, sys.ekin+sys.epot])[1:-1]+"\n")
    traj.write(str(sys.natoms)+"\n nfi="+str(sys.nfi)+" etot="+str(sys.ekin+sys.epot)+"\n")
    for i in range(0, sys.natoms):
        traj.write("Ar "+str(sys.r_x[i])+str(sys.r_y[i])+str(sys.r_z[i])+"\n")
