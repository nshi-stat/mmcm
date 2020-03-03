/*
    Copyright (C) 2009-2011, Kengo NAGASHIMA and Yasunori SATO.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

//**********************************************************************************/
// file name    mmcm.resamp.c
// purpose      MMCM wrap function for R shared library
// license      GPL-3
// Copyright (c) 2009-2011, Kengo NAGASHIMA and Yasunori SATO.
//**********************************************************************************/
#define MYMAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MYMIN(X, Y) (((X) < (Y)) ? (X) : (Y))

#include <R.h>
#include "mmcm.resamp.h"
#ifdef _OPENMP
#include <omp.h>
#endif


//**********************************************************************************/
// function name  mmcm_rwrap
// purpose        MMCM wrap function for R
// argument       *rdat        - for rdat.clsrnd (class information; before resampling)
//                             - clsrnd[sample_size]
//                *param       - for rdat.param (measurement)
//                             - param[sample_size]
//                *ctr_mat     - coefficient matrix of contrast statistics
//                             - ctr_mat[class_dim * contr_dim]
//                *resamp_size - No. of resampling repetition
//                *class_dim   - No. of groups
//                *contr_dim   - No. of contrast
//                *sample_size - whole sample size
//                *abseps      - absolute error tolerance
// return value   *pval        - resampling P-Value (count)
//                             - p[contr_dim]
//                *error       - estimated absolute error, with 99% confidence level.
//**********************************************************************************/
void mmcm_rwrap(double *param, double *clsrnd, double *ctr_mat, int *resamp_size,
    int *class_dim, int *contr_dim, int *sample_size, double *abseps, int *side,
    int *nthread, double *pval, double *error) {

  int i;
  struct mmcmdat *rdat;

  rdat = my_malloc_mmcmdat1(*sample_size);

  for (i = 0; i < *sample_size; i++) {
    rdat[i].param  = param[i];
    rdat[i].clsrnd = clsrnd[i];
  }

  // !!!!!!!!!!
  // calculate resampling P-Value (count) for MMCM
  //   argument
  //     dataset, coefficient matrix, resampling size
  //     No. of class, No. of contrast, sample size
  stat_resamp(rdat, ctr_mat, *resamp_size, *class_dim, *contr_dim,
    *sample_size, *abseps, *side, *nthread, pval, error);

  /*
  my_free_mmcmdat1(rdat);
  */

}

//**********************************************************************************/
// function name   stat_resamp
// purpose         calculate resampling P-Value (count) for MMCM
// argument        *rdat       - dataset
//                             - rdat.clsrnd (class information; before resampling)
//                             - rdat.param (measurement)
//                             - rdat[sample_size]
//                 *ctr_mat    - coefficient matrix of contrast statistics
//                             - ctr_mat[class_dim * contr_dim]
//                 resamp_size - No. of resampling repetition
//                 class_dim   - No. of groups
//                 contr_dim   - No. of contrast
//                 sample_size - whole sample size
// return value    *pval       - resampling P-Value (count)
//                             - p[contr_dim]
// error code      always 0 yet, this version
//**********************************************************************************/
 int stat_resamp(struct mmcmdat *rdat, double *ctr_mat, long resamp_size,
    int class_dim, int contr_dim, int sample_size, double abseps, int side,
    int nthread, double *pval, double *error) {

  int i, class_max = -1;
  int nmin = 1000;
  long j, *class_size, k = 0, pcount = 0;
  double all_sum = 0.0, *class_mean, *ctr_denom, *t_d, er = 1.0e+10, pv = 0.0;

  class_size = my_malloc_long1(class_dim);
  class_mean = my_malloc_double1(class_dim);
  ctr_denom  = my_malloc_double1(contr_dim);
  t_d        = my_malloc_double1(contr_dim);

  // calculate denominator of modified contrast statistics
  stat_denom(ctr_mat, class_dim, contr_dim, ctr_denom);

  // calculate
  //   class of the maximum sample size
  //   sample_size of each class
  //   sum of all measurement
  qsort(rdat, sample_size, sizeof(struct mmcmdat), clsrnd_cmp);
  gsize_chk(rdat, class_dim, sample_size, &class_max, class_size, &all_sum);

  // calculate maximum modified contrast statistics (sample)
  //   \bar{Y} = class_mean
  //   T'_k = frac{c^t_k \bar{Y}}{\sqrt{c^t_k c_k}}
  //   T'_max = t_d[0]
  mean_get(rdat, class_dim, sample_size, class_max, class_size, all_sum,
      class_mean);
  stat_calc(ctr_mat, ctr_denom, class_mean, class_dim, contr_dim, t_d);
  for (i = 0; i < contr_dim; i++) {
    switch (side) {
      case 1:
        t_d[0] = MYMIN(t_d[0], t_d[i]);
        break;
      case 2:
        t_d[0] = MYMAX(t_d[0], t_d[i]);
        break;
      case 3:
        t_d[0] = MYMAX(fabs(t_d[0]), fabs(t_d[i]));
        break;
      default:
        t_d[0] = MYMAX(fabs(t_d[0]), fabs(t_d[i]));
        break;
    }
  }
  
  
  // if openmp parallelization
  #ifdef _OPENMP
  struct mmcmdat **MP_rdat;
  int th_id;
  double **MP_class_mean, **MP_rt_d, **MP_randseq;
  
  omp_set_num_threads(nthread);
  
  #pragma omp parallel shared(k,pv,er,pcount,MP_rdat,MP_randseq,MP_class_mean,MP_rt_d,nthread) firstprivate(nmin,abseps,sample_size,class_dim,class_max,class_size,all_sum,ctr_mat,ctr_denom,contr_dim) private(i,j,th_id)
  {

    #pragma omp single
    {

      MP_rdat       = my_malloc_mmcmdat2(nthread, sample_size);
      MP_randseq    = my_malloc_double2(resamp_size, sample_size);
      MP_class_mean = my_malloc_double2(nthread, class_dim);
      MP_rt_d       = my_malloc_double2(nthread, contr_dim);
      
      GetRNGstate();
      for (j = 0; j < resamp_size; j++) {
        for (i = 0; i < sample_size; i++) {
          MP_randseq[j][i] = unif_rand();
        }
      }
      PutRNGstate();

      for (j = 0; j < nthread; j++) {
        for (i = 0; i < sample_size; i++) {
          MP_rdat[j][i].param  = rdat[i].param;
          MP_rdat[j][i].clsrnd = rdat[i].clsrnd;
        }
      }
    }

    #pragma omp for schedule(guided)
    for (j = 0; j < resamp_size; j++) {

      th_id = omp_get_thread_num();

      if (k > nmin && er < abseps) continue;

      // generate psude-random number
      for (i = 0; i < sample_size; i++) {
        MP_rdat[th_id][i].clsrnd = MP_randseq[j][i];
      }

      // calculate modified contrast statistics (resampling)
      rmean_get(&MP_rdat[th_id][0], class_dim, sample_size, class_max, class_size, all_sum,
        &MP_class_mean[th_id][0]);
      stat_calc(ctr_mat, ctr_denom, &MP_class_mean[th_id][0], class_dim, contr_dim, &MP_rt_d[th_id][0]);

      // calculate modified maximum contrast statistics (resampling)
      //   T'_{max} = \max \{ T'_k \}
      for (i = 0; i < contr_dim; i++) {
        switch (side) {
          case 1:
            MP_rt_d[th_id][0] = MYMIN(MP_rt_d[th_id][0], MP_rt_d[th_id][i]);
            break;
          case 2:
            MP_rt_d[th_id][0] = MYMAX(MP_rt_d[th_id][0], MP_rt_d[th_id][i]);
            break;
          case 3:
            MP_rt_d[th_id][0] = MYMAX(fabs(MP_rt_d[th_id][0]), fabs(MP_rt_d[th_id][i]));
            break;
          default:
            MP_rt_d[th_id][0] = MYMAX(fabs(MP_rt_d[th_id][0]), fabs(MP_rt_d[th_id][i]));
            break;
        }
      }

      // compare "resampling" T'_{max} > "sample" T'_i
      //   return value is counts of "resampling" T'_{max} > "sample" T'_{max}
      #pragma omp critical
      {
        switch (side) {
          case 1:
            if (MP_rt_d[th_id][0] < t_d[0]) {
              pcount++;
            }
            break;
          case 2:
            if (MP_rt_d[th_id][0] > t_d[0]) {
              pcount++;
            }
            break;
          case 3:
            if (MP_rt_d[th_id][0] > t_d[0]) {
              pcount++;
            }
            break;
          default:
            if (MP_rt_d[th_id][0] > t_d[0]) {
              pcount++;
            }
            break;
        }
        k++;
        pv = (double) pcount / (double) k;
        er = 3.5 * sqrt(pv * (1 - pv) / (double) k);
      }

    }

    /*
    #pragma omp barrier
    #pragma omp single
    {
      my_free_mmcmdat2(MP_rdat, nthread);
      my_free_double2(MP_randseq, resamp_size);
      my_free_double2(MP_class_mean, nthread);
      my_free_double2(MP_rt_d, nthread);
    }
    */

  }

  #else

  double *rt_d;

  rt_d       = my_malloc_double1(contr_dim);

  GetRNGstate();
  for (j = 0; j < resamp_size; j++) {

    if (k > nmin && er < abseps) continue;

    for (i = 0; i < sample_size; i++) {
      rdat[i].clsrnd = unif_rand();
    }

    // calculate modified contrast statistics (resampling)
    //   approach 1 (this program)
    //     rmean_get -> stat_calc
    //   approach 2 (common method)
    //     sort by rdat[i].clsrnd -> mean_get -> stat_calc
    rmean_get(rdat, class_dim, sample_size, class_max, class_size, all_sum,
      class_mean);
    stat_calc(ctr_mat, ctr_denom, class_mean, class_dim, contr_dim, rt_d);

    // calculate modified maximum contrast statistics (resampling)
    //   T'_{max} = \max \{ T'_k \}
    for (i = 0; i < contr_dim; i++) {
      switch (side) {
        case 1:
          rt_d[0] = MYMIN(rt_d[0], rt_d[i]);
          break;
        case 2:
          rt_d[0] = MYMAX(rt_d[0], rt_d[i]);
          break;
        case 3:
          rt_d[0] = MYMAX(fabs(rt_d[0]), fabs(rt_d[i]));
          break;
        default:
          rt_d[0] = MYMAX(fabs(rt_d[0]), fabs(rt_d[i]));
          break;
      }
    }

    // compare "resampling" T'_{max} > "sample" T'_i
    //   return value is counts of "resampling" T'_{max} > "sample" T'_{max}
    switch (side) {
      case 1:
        if (rt_d[0] < t_d[0]) {
          pcount++;
        }
        break;
      case 2:
        if (rt_d[0] > t_d[0]) {
          pcount++;
        }
        break;
      case 3:
        if (rt_d[0] > t_d[0]) {
          pcount++;
        }
        break;
      default:
        if (rt_d[0] > t_d[0]) {
          pcount++;
        }
        break;
    }
    k++;
    pv = (double) pcount / (double) k;
    er = 3.5 * sqrt(pv * (1 - pv) / (double) k);

  }
  PutRNGstate();

  #endif

  *pval = (double) pcount / (double) k;
  *error = 3.5 * sqrt(*pval * (1 - *pval) / (double) k);

  /*
  my_free_long1(class_size);
  my_free_double1(class_mean);
  my_free_double1(ctr_denom);
  my_free_double1(t_d);
  my_free_double1(rt_d);
  */

  return 0;

}

//**********************************************************************************/
// function name   stat_calc
// purpose         calculate "modified contrast statistics"
// argument        *ctr_mat    - coefficient matrix of contrast statistics
//                             - ctr_mat[class_dim * contr_dim]
//                 *ctr_denom  - denominator of modified contrast statistics
//                             - ctr_denom[contr_dim]
//                 *class_mean - "sample" or "resampling" mean of each class
//                             - class_mean[class_dim]
//                 class_dim   - No. of groups
//                 contr_dim   - No. of contrast
// return value    *t_d        - modified contrast statistics
//                             - t_d[contr_dim]
// error code      always 0 yet, this version
//**********************************************************************************/
int stat_calc(double *ctr_mat, double *ctr_denom, double *class_mean,
    int class_dim, int contr_dim, double *t_d) {

  int i, j;
  double sum;

  // modified contrast statistics
  // T'_k = frac{c^t_k \bar{Y}}{\sqrt{c^t_k c_k}}
  for (i = 0; i < contr_dim; i++) {
    sum = 0;
    for (j = 0; j < class_dim; j++) {
      sum += class_mean[j] * ctr_mat[class_dim * i + j];
    }

    t_d[i] = sum / ctr_denom[i];

  }

  return 0;

}

//**********************************************************************************/
// function name   stat_denom
// purpose         calculate denominator of "modified contrast statistics"
// argument        *ctr_mat   - coefficient matrix of contrast statistics
//                            - ctr_mat[class_dim * contr_dim]
//                 class_dim  - No. of groups
//                 contr_dim  - No. of contrast
// return value    *ctr_denom - denominator of modified contrast statistics
//                            - ctr_denom[contr_dim]
// error code      always 0 yet, this version
//**********************************************************************************/
int stat_denom(double *ctr_mat, int class_dim, int contr_dim, double *ctr_denom) {

  int i, j;
  double sum;

  // denominator of modified contrast statistics
  // \sqrt{c^t_k c_k} = \sum_{i=1}^{contr_dim} c_{i}^2
  for (i = 0; i < contr_dim; i++) {
    sum = 0;
    for (j = 0; j < class_dim; j ++) {
      sum += pow(ctr_mat[class_dim * i + j], 2);
    }
    ctr_denom[i] = sqrt(sum);
  }

  return 0;

}

//**********************************************************************************/
// function name   gsize_chk
// purpose         calculate sample size of each class
//                 calculate class has the maximum sample size
//                 calculate sum of all measurements
// argument        *rdat      - dataset
//                            - rdat.clsrnd (class information; before resampling)
//                            - rdat.param (measurement)
//                            - rdat[sample_size]
//                class_dim   - No. of groups
//                sample_size - whole sample size
// return value   *class_max  - class of the maximum sample size
//                *class_size - sample_size of each class
//                            - class_size[class_dim]
//                *all_sum    - sum of all measurement
// error code     always 0 yet, this version
//**********************************************************************************/
int gsize_chk(struct mmcmdat *rdat, int class_dim, int sample_size,
    int *class_max, long *class_size, double *all_sum) {

  long i, max = 0, tmp = 1;
  double sum = 0;

  // calculate class_size[i] & all_sum
  for (i = 0; i < sample_size; i++) {
    sum += rdat[i].param;
    if (tmp - rdat[i].clsrnd < 0.0) {
      class_size[tmp-1] = i;
      tmp++;
    }
  }
  class_size[class_dim-1] = sample_size;

  *all_sum = sum;

  // calculate class_max & cnvert class_size[i]
  for (i = class_dim - 1; i > 0; i--) {
    class_size[i] = class_size[i] - class_size[i-1];
    if (max < class_size[i]) {
      max = class_size[i];
      *class_max = (int)i;
    }
  }
  if (max < class_size[0]) {
    *class_max = 0;
  }

  return 0;

}

//**********************************************************************************/
// function name   mean_get
// purpose         calculate "sample" mean of each class (class_mean[i])
// argument        *rdat       - dataset
//                             - rdat.param (measurement)
//                             - rdat[sample_size]
//                 class_dim   - No. of groups
//                 sample_size - whole sample size
//                 class_max   - class of the maximum sample size
//                 *class_size - sample_size of each class
//                             - class_size[class_dim]
//                 all_sum     - sum of all measurement
// return value    *class_mean - "sample" mean of each class
//                             - class_mean[class_dim]
// error code      always 0 yet, this version
//**********************************************************************************/
int mean_get(struct mmcmdat *rdat, int class_dim, int sample_size,
    int class_max, long *class_size, double all_sum, double *class_mean) {

  int i;
  long j, nsum = sample_size, start = 0, end = 0;
  double sum = 0, dsum = all_sum;

  // calculate [i] class mean by class_size[i]
  for (i = 0; i < class_dim; i++) {
    end += class_size[i];
    sum = 0;
    if (i != class_max) {
      for (j = start; j < end; j++) {
        sum += rdat[j].param;
      }
      class_mean[i] = sum / class_size[i];
      dsum -= sum;
      nsum -= class_size[i];
    }
    start += class_size[i];
  }

  class_mean[class_max] = dsum / nsum;

  return 0;

}

//**********************************************************************************/
// function name   rmean_get
// purpose         calculate "resampling" mean of each class (class_mean[i])
// argument        *rdat       - dataset
//                             - rdat.clsrnd (pseudo-random number for resampling)
//                             - rdat.param (measurement)
//                             - rdat[sample_size]
//                 class_dim   - No. of groups
//                 sample_size - whole sample size
//                 class_max   - class of the maximum sample size
//                 *class_size - sample_size of each class
//                             - class_size[class_dim]
//                 all_sum     - sum of all measurement
// return value    *class_mean - "resampling" mean of each class
//                             - class_mean[class_dim]
// error code      always 0 yet, this version
//**********************************************************************************/
int rmean_get(struct mmcmdat *rdat, int class_dim, int sample_size,
    int class_max, long *class_size, double all_sum, double *class_mean) {

  int i;
  long j, k, pre_param = -1, now_param = 0, nsum = sample_size;
  double sum, dsum = all_sum;

  for (i = 0; i < class_dim; i++) {

    sum = 0;

    if (i != class_max) {

      for (j = 0; j < class_size[i]; j++) {

        // only first sample
        // select maximum pseudo-random number
        if (pre_param == -1) {

          for (k = 0; k < sample_size; k++) {
            if (rdat[now_param].clsrnd < rdat[k].clsrnd)
              now_param = k;
          }
        }
        // other sample
        // select next largest pseudo-random number
        else {
          for (k = 0; k < sample_size; k++) {
            if (rdat[pre_param].clsrnd > rdat[k].clsrnd) {
              now_param = k;
              break;
            }
          }
          for (k = 0; k < sample_size; k++) {
            if (rdat[now_param].clsrnd < rdat[k].clsrnd
                && rdat[pre_param].clsrnd > rdat[k].clsrnd)
              now_param = k;
          }
        }
        sum += rdat[now_param].param;
        pre_param = now_param;
      }

      // calculate [i] class mean
      class_mean[i] = sum / class_size[i];

      dsum -= sum;
      nsum -= class_size[i];

    }
  }

  // calculat last class mean (class has the maximum sample size)
  // reduce calculation amount
  class_mean[class_max] = dsum / nsum;

  return 0;

}

//**********************************************************************************/
// function name  clsrnd_cmp
// purpose        comparison function in qsort()
//                sort by mmcmdat.clsrnd value in ascending order
//**********************************************************************************/
int clsrnd_cmp(const void *_p0, const void *_p1) {
  struct mmcmdat *p0 = (struct mmcmdat *)_p0;
  struct mmcmdat *p1 = (struct mmcmdat *)_p1;
  if (p0->clsrnd < p1->clsrnd)
    return -1; // in ascending order
  else if (p0->clsrnd > p1->clsrnd)
    return 1; // in ascending order
  else
    return 0;
}

//**********************************************************************************/
// function name  my_malloc_mmcmdat1
// purpose        memory allocation function for 1 dimensional mmcmdat array
//**********************************************************************************/
struct mmcmdat *my_malloc_mmcmdat1(int m) {

  struct mmcmdat *p;
  if ((p = (struct mmcmdat *) R_alloc(m, sizeof(struct mmcmdat))) == NULL) {
    error("Error: my_malloc_mmcmdat1: failed to allocate size %d mmcmdat array\n", m);
  }
  return p;

}

//**********************************************************************************/
// function name  my_malloc_mmcmdat2
// purpose        memory allocation function for 2 dimensional mmcmdat array
//**********************************************************************************/
struct mmcmdat **my_malloc_mmcmdat2(int m, int n) {

  struct mmcmdat **p;

  if ((p = (struct mmcmdat **) R_alloc(m, sizeof(struct mmcmdat *))) == NULL) {
    error("Error: my_malloc_mmcmdat2: failed to allocate size %d mmcmdat * array\n", m);
  }
	for (int i = 0; i < m; i++) {
    if ((p[i] = (struct mmcmdat *) R_alloc(n, sizeof(struct mmcmdat))) == NULL) {
      error("Error: my_malloc_mmcmdat2: failed to allocate size %d mmcmdat array\n", n);
    }
	}
  return p;

}

//**********************************************************************************/
// function name  my_malloc_long1
// purpose        memory allocation function for 1 dimensional long array
//**********************************************************************************/
long *my_malloc_long1(int m) {

  long *p;
  if ((p = (long *) R_alloc(m, sizeof(long))) == NULL) {
    error("Error: my_malloc_long1: failed to allocate size %d long array\n", m);
  }
  return p;

}

//**********************************************************************************/
// function name  my_malloc_double1
// purpose        memory allocation function for 1 dimensional double array
//**********************************************************************************/
double *my_malloc_double1(int m) {

  double *p;
  if ((p = (double *) R_alloc(m, sizeof(double))) == NULL) {
    error("Error: my_malloc_double1: failed to allocate size %d double array\n", m);
  }
  return p;

}

//**********************************************************************************/
// function name  my_malloc_double2
// purpose        memory allocation function for 2 dimensional double array
//**********************************************************************************/
double **my_malloc_double2(int m, int n) {

  double **p;

  if ((p = (double **) R_alloc(m, sizeof(double *))) == NULL) {
    error("Error: my_malloc_double2: failed to allocate size %d double * array\n", m);
  }
	for (int i = 0; i < m; i++) {
    if ((p[i] = (double *) R_alloc(n, sizeof(double))) == NULL) {
      error("Error: my_malloc_double2: failed to allocate size %d double array\n", n);
    }
	}
  return p;

}

//**********************************************************************************/
// function name  my_free_mmcmdat1
// purpose        memory free function for 1 dimensional mmcmdat array
//**********************************************************************************/
void my_free_mmcmdat1(struct mmcmdat *p) {
  
	if (p) {
		free(p);
	}
	
}

//**********************************************************************************/
// function name  my_free_mmcmdat1
// purpose        memory free function for 2 dimensional mmcmdat array
//**********************************************************************************/
void my_free_mmcmdat2(struct mmcmdat **p, int m){
	
	if (p) {
  	for (int i = 0; i < m; i++) {
    	if (p[i]) {
    		free(p[i]);
    	}
  	}
		free(p);
	}
	
}

//**********************************************************************************/
// function name  my_free_long1
// purpose        memory free function for 1 dimensional long array
//**********************************************************************************/
void my_free_long1(long *p) {
  
	if (p) {
		free(p);
	}
	
}

//**********************************************************************************/
// function name  my_free_double1
// purpose        memory free function for 1 dimensional double array
//**********************************************************************************/
void my_free_double1(double *p) {
  
	if (p) {
		free(p);
	}
	
}

//**********************************************************************************/
// function name  my_free_double2
// purpose        memory free function for 2 dimensional double array
//**********************************************************************************/
void my_free_double2(double **p, int m) {
	
	if (p) {
  	for (int i = 0; i < m; i++) {
    	if (p[i]) {
    		free(p[i]);
    	}
  	}
		free(p);
	}
	
}
