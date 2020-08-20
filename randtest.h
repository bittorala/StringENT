
#ifndef MONTEN
#define MONTEN 6 
/* Bytes used as Monte Carlo co-ordinates.	
This should be no more bits than the mantissa of your "double" floating point type. */

#endif
/*  Random test function prototypes  */
extern void rt_init(int binmode);
extern void rt_add(void *buf, int bufl);
extern void rt_end(double *r_ent, double *r_chisq, double *r_mean,
                   double *r_montepicalc, double *r_scc, long *r_nruns,
                   long *r_n1, long *r_n2, double *r_lmx2, int *r_lmnblocks);
