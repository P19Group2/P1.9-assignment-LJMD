#include "force.h"

/* helper function: apply minimum image convention */
double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

// helper function: zero out an array
void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

/* compute forces */
void force(mdsys_t *sys) 
{
    double r,ffac;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    for(i=0; i < (sys->natoms); ++i) {
        for(j=0; j < (sys->natoms); ++j) {

            /* particles have no interactions with themselves */
            if (i==j) continue;
            
            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            r = sqrt(rx*rx + ry*ry + rz*rz);
      
            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
		double rinv = 1/r;
                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma*rinv,12.0)*rinv
                                         +6*pow(sys->sigma*rinv,6.0)*rinv);
                
                sys->epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma*rinv,12.0)
                                               -pow(sys->sigma*rinv,6.0));

                sys->fx[i] += rx*rinv*ffac;
                sys->fy[i] += ry*rinv*ffac;
                sys->fz[i] += rz*rinv*ffac;
            }
        }
    }
}
