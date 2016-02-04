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
    double rsq,ffac;
    double rx,ry,rz;
    double b = 0.5*sys->box;
    double csq = sys->rcut*sys->rcut;            
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
            rx=pbc(sys->rx[i] - sys->rx[j], b);
            ry=pbc(sys->ry[i] - sys->ry[j], b);
            rz=pbc(sys->rz[i] - sys->rz[j], b);
            rsq = rx*rx + ry*ry + rz*rz;
      
            /* compute force and energy if within cutoff */
            if (rsq < csq) {
		double rsqinv = 1/rsq;
		double r6 = rsqinv*rsqinv*rsqinv;
		double s = sys->sigma*sys->sigma;
		double p = s*s*s*r6;
		double p6 = 4.0*sys->epsilon*p;
                double p12 = p6*p;

		ffac = rsqinv*(12*p12-6*p6);
                
                sys->epot += 0.5*(p12-p6);

                sys->fx[i] += rx*ffac;
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;
            }
        }
    }
}
