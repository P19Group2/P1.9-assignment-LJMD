#include "force.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
/* helper function: apply minimum image convention */
inline double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

// helper function: zero out an array
inline void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

#ifdef _OPENMP
/* compute MORSE forces */
void force_morse(mdsys_t *sys) 
{
    double pot;
    double b = 0.5*sys->box;
    double s = sys->sigma*sys->sigma;
    double p = s*s*s;
    double csq = sys->rcut*sys->rcut;            
    int i;
    int natoms=sys->natoms;

    /* zero energy and forces */
    pot=0.0;
    azzero(sys->fx,natoms);
    azzero(sys->fy,natoms);
    azzero(sys->fz,natoms);

    #pragma omp parallel for default(shared) private(i) reduction(+:pot)
    for(i = 0; i < (natoms - 1); ++i) {
    /* private versions of sys->rx[i] */
    int j;
    double rxi,ryi,rzi;
    rxi = sys->rx[i];
    ryi = sys->ry[i];
    rzi = sys->rz[i];
    for(j = (1 + i); j < natoms; ++j) {
    double rx,ry,rz,rsq;
    /* particles have no interactions with themselves */
    if (i==j) continue;

    /* get distance between particle i and j */
    rx=pbc(rxi - sys->rx[j], b);
    ry=pbc(ryi - sys->ry[j], b);
    rz=pbc(rzi - sys->rz[j], b);
    rs = sqrt(rsqi);
    rsq = rx*rx + ry*ry + rz*rz;
            
	/* compute force and energy if within cutoff */
    if (rsq < csq) {
				
				esp = exp((-1)*sys->a*(rs-sys->rzero));
				espq = exp((-1)*sys->a*(rs-sys->rzero)*(rs-sys->rzero));
				
				ffac = 2*sys->De * (1-esp)*sys->a*esp;
				pot += sys->De*(1-espq)
                
				sys->fx[i] += rx*ffac;
                #pragma omp atomic
                sys->fy[i] += ry*ffac;
                #pragma omp atomic
                sys->fz[i] += rz*ffac;
                #pragma omp atomic
                sys->fx[j] -= rx*ffac;
                #pragma omp atomic
                sys->fy[j] -= ry*ffac;
                #pragma omp atomic
                sys->fz[j] -= rz*ffac;
            }
        }
    } 
    sys->epot = pot;
}

/* compute LJ forces */
void force_LJ(mdsys_t *sys) 
{
    double pot;
    double b = 0.5*sys->box;
    double s = sys->sigma*sys->sigma;
    double p = s*s*s;
    double csq = sys->rcut*sys->rcut;            
    int i;
    int natoms=sys->natoms;

    /* zero energy and forces */
    pot=0.0;
    azzero(sys->fx,natoms);
    azzero(sys->fy,natoms);
    azzero(sys->fz,natoms);

    #pragma omp parallel for default(shared) private(i) reduction(+:pot)
    for(i = 0; i < (natoms - 1); ++i) {
    /* private versions of sys->rx[i] */
    int j;
    double rxi,ryi,rzi;
    rxi = sys->rx[i];
    ryi = sys->ry[i];
    rzi = sys->rz[i];
        for(j = (1 + i); j < natoms; ++j) {
        double rx,ry,rz,rsq;
            /* particles have no interactions with themselves */
            if (i==j) continue;

            /* get distance between particle i and j */
            rx=pbc(rxi - sys->rx[j], b);
            ry=pbc(ryi - sys->ry[j], b);
            rz=pbc(rzi - sys->rz[j], b);
            rsq = rx*rx + ry*ry + rz*rz;
      
            /* compute force and energy if within cutoff */
            if (rsq < csq) {
		double rsqinv = 1/rsq;
		double r6 = rsqinv*rsqinv*rsqinv;
		double pr = p*r6;
		double p6 = 4.0*sys->epsilon*pr;
                double p12 = p6*pr; 
                double ffac;
		
                ffac = rsqinv*(12*p12-6*p6);
                
                pot += 0.5*(p12-p6);
                
                #pragma omp atomic
                sys->fx[i] += rx*ffac;
                #pragma omp atomic
                sys->fy[i] += ry*ffac;
                #pragma omp atomic
                sys->fz[i] += rz*ffac;
                #pragma omp atomic
                sys->fx[j] -= rx*ffac;
                #pragma omp atomic
                sys->fy[j] -= ry*ffac;
                #pragma omp atomic
                sys->fz[j] -= rz*ffac;
            }
        }
    } 
    sys->epot = pot;
}
#else

/*Compute morse forces*/

void force_morse(mdsys_t *sys) 
{
    double rsq,ffac;
    double rx,ry,rz;
    double b = 0.5*sys->box;
    double s = sys->sigma*sys->sigma;
    double p = s*s*s;
    double csq = sys->rcut*sys->rcut;            
    int i,j;
    int natoms=sys->natoms;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,natoms);
    azzero(sys->fy,natoms);
    azzero(sys->fz,natoms);

    for(i=0; i < natoms-1; ++i) {
        for(j=i+1; j < natoms; ++j) {

            /* particles have no interactions with themselves */
            if (i==j) continue;

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], b);
            ry=pbc(sys->ry[i] - sys->ry[j], b);
            rz=pbc(sys->rz[i] - sys->rz[j], b);
            rsq = rx*rx + ry*ry + rz*rz;
			rs = sqrt(rsqi);
			rsq = rx*rx + ry*ry + rz*rz;

	       /* compute force and energy if within cutoff */
	       if (rsq < csq) {

			   esp = exp((-1)*sys->a*(rs-sys->rzero));
			   espq = exp((-1)*sys->a*(rs-sys->rzero)*(rs-sys->rzero));			 
			   
			   ffac = 2*sys->De * (1-esp)*sys->a*esp;
			   sys->epot += sys->De*(1-espq)


                sys->fx[i] += rx*ffac;
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;
                sys->fx[j] -= rx*ffac;
                sys->fy[j] -= ry*ffac;
                sys->fz[j] -= rz*ffac;

            }
        }
    }
}

/*Compute LJ forces*/

void force_LJ(mdsys_t *sys) 
{
    double rsq,ffac;
    double rx,ry,rz;
    double b = 0.5*sys->box;
    double s = sys->sigma*sys->sigma;
    double p = s*s*s;
    double csq = sys->rcut*sys->rcut;            
    int i,j;
    int natoms=sys->natoms;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,natoms);
    azzero(sys->fy,natoms);
    azzero(sys->fz,natoms);

    for(i=0; i < natoms-1; ++i) {
        for(j=i+1; j < natoms; ++j) {

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
		double pr = p*r6;
		double p6 = 4.0*sys->epsilon*pr;
                double p12 = p6*pr;

		ffac = rsqinv*(12*p12-6*p6);
                
                sys->epot += 0.5*(p12-p6);

                sys->fx[i] += rx*ffac;
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;
                sys->fx[j] -= rx*ffac;
                sys->fy[j] -= ry*ffac;
                sys->fz[j] -= rz*ffac;

            }
        }
    }
}
#endif
