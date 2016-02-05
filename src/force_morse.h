#ifndef FORCE_H
#define FORCE_H

#include "mdsys.h"
#include <math.h>

void force_LJ(mdsys_t *);  
void force_morse(mdsys_t *);

void azzero(double *, const int);

double pbc(double, const double);

#endif
