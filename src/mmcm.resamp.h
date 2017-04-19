//**********************************************************************************/
// file name	mmcm.resamp.h
// purpose		MMCM wrap function for R shared library (header file)
// license    GPL-3
// Copyright (c) 2009-2011, Kengo NAGASHIMA and Yasunori SATO.
//**********************************************************************************/
struct mmcmdat {
	double param;
	double clsrnd;
};

void mmcm_rwrap(double *, double *, double *, int *, int *, int *, int *, double *, int *, double *, double *);
int stat_resamp(struct mmcmdat *, double *, long, int, int, int, double, int, double *, double *);
int mean_get(struct mmcmdat *, int, int, int, long *, double, double *);
int rmean_get(struct mmcmdat *, int, int, int, long *, double, double *);
int gsize_chk(struct mmcmdat *, int, int, int *, long *, double *);
int stat_calc(double *, double *, double *, int, int, double *);
int stat_denom(double *, int, int, double *);
int clsrnd_cmp(const void *_p0, const void *_p1);

struct mmcmdat *my_malloc_mmcmdat1(int);
struct mmcmdat **my_malloc_mmcmdat2(int, int);
long   *my_malloc_long1(int);
double *my_malloc_double1(int);
double **my_malloc_double2(int, int);

void my_free_mmcmdat1(struct mmcmdat *);
void my_free_mmcmdat2(struct mmcmdat **, int);
void my_free_long1(long *);
void my_free_double1(double *);
void my_free_double2(double **, int);
