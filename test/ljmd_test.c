#include "../src/mdsys.h"
#include "../src/force.h"
#include "../src/utilities.h"
#include "../src/velverlet.h"
#include<math.h>
#include<assert.h>
#include<stdio.h>

const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

int main(){

    mdsys_t sys;
    sys.natoms=2;
    sys.mass=39.948;
    sys.epsilon=0.2379;
    sys.sigma=3.405;
    sys.rcut=30;
    sys.box=17.1580;
    sys.nsteps=1;
    
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.natoms*sizeof(double));
    sys.fy=(double *)malloc(sys.natoms*sizeof(double));
    sys.fz=(double *)malloc(sys.natoms*sizeof(double));
    
    sys.rx[0]=6.67103294321331;
    sys.ry[0]=-10.6146871435653;
    sys.rz[0]=12.6336939877734;
    sys.rx[1]=1.06574058650169;
    sys.ry[1]=-3.33432278188177;
    sys.rz[1]=-2.59038677851747;
    
    sys.vx[0]=-1.5643224621482283e-03;
    sys.vy[0]=4.8497508563925346e-04;
    sys.vz[0]=-4.3352481732883966e-04;
    sys.vx[1]=4.1676710257651452e-04;
    sys.vy[1]=2.2858522230176587e-05;
    sys.vz[1]=-6.1985040462745732e-04;
    
    azzero(sys.fx, sys.natoms);
    azzero(sys.fy, sys.natoms);
    azzero(sys.fz, sys.natoms);

    sys.nfi=0;
    
    force(&sys);
    assert(abs(sys.fx[0]+0.000822)<=1e-7);
    assert(abs(sys.fy[0]-0.001067)<=1e-7);
    assert(abs(sys.fz[0]-0.000284)<=1e-7);
    assert(abs(sys.fx[1]-0.000822)<=1e-7);
    assert(abs(sys.fy[1]+0.001067)<=1e-7);
    assert(abs(sys.fz[1]+0.000284)<=1e-7);
    printf("Force test passed\n");
    
    ekin(&sys);
    assert(abs(sys.ekin-0.163682)<=1e-7);
    assert(abs(sys.temp-54.911863)<=1e-7);
    printf("Ekin test passed\n");
    
    velverlet_1(&sys);
    assert(abs(sys.rx[0]-6.671033)<=1e-7);
    assert(abs(sys.ry[0]+10.614687)<=1e-7);
    assert(abs(sys.rz[0]-12.633694)<=1e-7);
    assert(abs(sys.vx[0]+0.001564)<=1e-7);
    assert(abs(sys.vy[0]-0.000485)<=1e-7);
    assert(abs(sys.vz[0]+0.000434)<=1e-7);
    assert(abs(sys.rx[1]-1.065741)<=1e-7);
    assert(abs(sys.ry[1]+3.334323)<=1e-7);
    assert(abs(sys.rz[1]+2.590387)<=1e-7);
    assert(abs(sys.vx[1]-0.000417)<=1e-7);
    assert(abs(sys.vy[1]-0.000023)<=1e-7);
    assert(abs(sys.vz[1]+0.000620)<=1e-7);
    printf("Velverlet_1 test passed\n");    
    
    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);

    return 0;
}

