/* 
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

/* generic file- or pathname buffer length */
#define BLEN 200

/* a few physical constants */
const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

/* structure to hold the complete information 
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
};
typedef struct _mdsys mdsys_t;

/* helper function: read a line and then return
   the first string with whitespace stripped off */
static int get_a_line(FILE *fp, char *buf)
{
    char tmp[BLEN], *ptr;

    /* read a line and cut of comments and blanks */
    if (fgets(tmp,BLEN,fp)) {
        int i;

        ptr=strchr(tmp,'#');
        if (ptr) *ptr= '\0';
        i=strlen(tmp); --i;
        while(isspace(tmp[i])) {
            tmp[i]='\0';
            --i;
        }
        ptr=tmp;
        while(isspace(*ptr)) {++ptr;}
        i=strlen(ptr);
        strcpy(buf,tmp);
        return 0;
    } else {
        perror("problem reading input");
        return -1;
    }
    return 0;
}
 
/* helper function: zero out an array */
__global__ void azzero(double *d, const int n)
{
    int threads = gridDim.x*blockDim.x;
    int threadId  = blockDim.x*blockIdx.x + threadIdx.x;

    int chunk = n/threads;
    int start = threadId * chunk;
    int end = (threadId + 1) * chunk;

    int i;
    for (i=start; i<end; ++i) {
        d[i]=0.0;
    }
}

/* helper function: apply minimum image convention */
__kernel__ double pbc(double x, const double boxby2)
{
    int threadId  = blockDim.x*blockIdx.x + threadIdx.x;
if(
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

/* compute kinetic energy */
static void ekin(mdsys_t *sys)
{   
    int i;

    sys->ekin=0.0;
    for (i=0; i<sys->natoms; ++i) {
        sys->ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}

/* compute forces */
__global__ void force(double *d_rx, double *d_ry, double *d_rz, double *d_fx, double *d_fy, double *d_fz, double d_epot, double d_sigma, double d_epsilon, double d_box, double d_rcut, int natoms) 
{
    int threads = gridDim.x*blockDim.x;
    int threadId  = blockDim.x*blockIdx.x + threadIdx.x;

    if (threads > natoms) {
        threads = natoms);
        if (threadId >= threads) {
            return;
        }
    }


    double r,ffac;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    d_epot=0.0;
    azzero(d_fx,natoms);
    azzero(d_fy,natoms);
    azzero(d_fz,natoms);
    int chunk = (natoms)/threads;
    int start = threadId * chunk;
    int end = (threadId + 1) * chunk;

    for(i=start; i < end; ++i) {
        for(j=0; j < natoms; ++j) {

            /* particles have no interactions with themselves */
            if (i==j) continue;
            
            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*d_box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*d_box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*d_box);
            r = sqrt(rx*rx + ry*ry + rz*rz);
      
            /* compute force and energy if within cutoff */
            if (r < d_rcut) {
                ffac = -4.0*d_epsilon*(-12.0*pow(d_sigma/r,12.0)/r
                                         +6*pow(d_sigma/r,6.0)/r);
                
                d_epot += 0.5*4.0*d_epsilon*(pow(d_sigma/r,12.0)
                                               -pow(d_sigma/r,6.0));

                d_fx[i] += rx/r*ffac;
                d_fy[i] += ry/r*ffac;
                d_fz[i] += rz/r*ffac;
            }
        }
    }
}

/* velocity verlet */
static void velverlet_1(mdsys_t *sys)
{
    int i;

    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->rx[i] += sys->dt*sys->vx[i];
        sys->ry[i] += sys->dt*sys->vy[i];
        sys->rz[i] += sys->dt*sys->vz[i];
    }
}
    /* compute forces and potential energy */
    //force(sys);
static void velverlet_2(mdsys_t *sys)
{
    int i;
    /* second part: propagate velocities by another half step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }
}

/* append data to output. */
static void output(mdsys_t *sys, FILE *erg, FILE *traj)
{
    int i;
    
    printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi, sys->ekin+sys->epot);
    for (i=0; i<sys->natoms; ++i) {
        fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i], sys->rz[i]);
    }
}


/* main */
int main(int argc, char **argv) 
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;

    /* read input file */
    if(get_a_line(stdin,line)) return 1;
    sys.natoms=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.mass=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.epsilon=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.sigma=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.rcut=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.box=atof(line);
    if(get_a_line(stdin,restfile)) return 1;
    if(get_a_line(stdin,trajfile)) return 1;
    if(get_a_line(stdin,ergfile)) return 1;
    if(get_a_line(stdin,line)) return 1;
    sys.nsteps=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.dt=atof(line);
    if(get_a_line(stdin,line)) return 1;
    nprint=atoi(line);

    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.natoms*sizeof(double));
    sys.fy=(double *)malloc(sys.natoms*sizeof(double));
    sys.fz=(double *)malloc(sys.natoms*sizeof(double));

    /* read restart */
    fp=fopen(restfile,"r");
    if(fp) {
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
        }
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
        }
        fclose(fp);
        azzero(sys.fx, sys.natoms);
        azzero(sys.fy, sys.natoms);
        azzero(sys.fz, sys.natoms);
    } else {
        perror("cannot read restart file");
        return 3;
    }

    /* initialize forces and energies.*/
    sys.nfi=0;
    cudaMalloc((void**) &d_rx, sys.natoms*sizeof(double));
    cudaMalloc((void**) &d_ry, sys.natoms*sizeof(double));
    cudaMalloc((void**) &d_rz, sys.natoms*sizeof(double));
    cudaMalloc((void**) &d_fx, sys.natoms*sizeof(double));
    cudaMalloc((void**) &d_fy, sys.natoms*sizeof(double));
    cudaMalloc((void**) &d_fz, sys.natoms*sizeof(double));
/*    cudaMalloc((void**) &d_epot, sizeof(double));
    cudaMalloc((void**) &d_sigma, sizeof(double));
    cudaMalloc((void**) &d_epsilon, sizeof(double));
    cudaMalloc((void**) &d_box, sizeof(double));
    cudaMalloc((void**) &d_rcut, sizeof(double)); */

    cudaMemcpy(d_rx, sys->rx, sys.natoms*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ry, sys->ry, sys.natoms*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rz, sys->rz, sys.natoms*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fx, sys->fx, sys.natoms*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fy, sys->fy, sys.natoms*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fz, sys->fz, sys.natoms*sizeof(double), cudaMemcpyHostToDevice);
/*    cudaMemcpy(d_epot, sys->epot, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sigma, sys->sigma, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_epsilon, sys->epsilon, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_box, sys->box, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rcut, sys->rcut, sizeof(double), cudaMemcpyHostToDevice); */

    force(&sys);


    cudaMemcpy(sys->rx, d_rx, sys.natoms*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->ry, d_ry, sys.natoms*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->rz, d_rz, sys.natoms*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->fx, d_fx, sys.natoms*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->fy, d_fy, sys.natoms*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->fz, d_fz, sys.natoms*sizeof(double), cudaMemcpyDeviceToHost);
/*    cudaMemcpy(sys->epot, d_epot, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->sigma, d_sigma, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->epsilon, d_epsilon, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->box, d_box, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->rcut, d_rcut, sizeof(double), cudaMemcpyDeviceToHost);*/


    ekin(&sys);
    
    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
        velverlet_1(&sys);
    cudaMemcpy(d_rx, sys->rx, sys.natoms*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ry, sys->ry, sys.natoms*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rz, sys->rz, sys.natoms*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fx, sys->fx, sys.natoms*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fy, sys->fy, sys.natoms*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fz, sys->fz, sys.natoms*sizeof(double), cudaMemcpyHostToDevice);
/*    cudaMemcpy(d_epot, sys->epot, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sigma, sys->sigma, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_epsilon, sys->epsilon, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_box, sys->box, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rcut, sys->rcut, sizeof(double), cudaMemcpyHostToDevice);*/
	force(&sys);
    cudaMemcpy(sys->rx, d_rx, sys.natoms*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->ry, d_ry, sys.natoms*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->rz, d_rz, sys.natoms*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->fx, d_fx, sys.natoms*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->fy, d_fy, sys.natoms*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->fz, d_fz, sys.natoms*sizeof(double), cudaMemcpyDeviceToHost);
/*    cudaMemcpy(sys->epot, d_epot, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->sigma, d_sigma, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->epsilon, d_epsilon, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->box, d_box, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sys->rcut, d_rcut, sizeof(double), cudaMemcpyDeviceToHost);*/
	velverlet_2(&sys);
        ekin(&sys);
    }
    /**************************************************/
    cudaFree(d_rx);
    cudaFree(d_ry);
    cudaFree(d_rz);
    cudaFree(d_fx);
    cudaFree(d_fy);
    cudaFree(d_fz);
/*    cudaFree(d_epot);
    cudaFree(d_sigma);
    cudaFree(d_epsilon);
    cudaFree(d_box);
    cudaFree(d_rcut);*/

    /* clean up: close files, free memory */
    printf("Simulation Done.\n");
    fclose(erg);
    fclose(traj);

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
