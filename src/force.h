#ifndef FORCE_H
#define FORCE_H

#include "mdsys.h"
#include <math.h>

typedef struct _mdsys mdsys_t;

void force(mdsys_t *);  

void azzero(double *, const int);

double pbc(double, const double);

#endif
