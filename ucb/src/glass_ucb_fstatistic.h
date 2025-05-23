/*
 * Copyright 2019 Travis Robson, Tyson B. Littenberg, Neil J. Cornish
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#ifndef ucb_fstatistic_h
#define ucb_fstatistic_h

#include <stdio.h>

struct Filter
{
    double *A1_fX, *A1_fY, *A1_fZ, *A1_fA, *A1_fE;
    double *A2_fX, *A2_fY, *A2_fZ, *A2_fA, *A2_fE;
    double *A3_fX, *A3_fY, *A3_fZ, *A3_fA, *A3_fE;
    double *A4_fX, *A4_fY, *A4_fZ, *A4_fA, *A4_fE;
    
    long M_filter, N_filter;
    
    double N1_X,  N2_X,  N3_X,  N4_X;
    double N1_AE, N2_AE, N3_AE, N4_AE;
    
    double **M_inv_X, **M_inv_AE;
    
    double Fstat_X, Fstat_AE;
    
    double a1_X, a1_AE;
    double a2_X, a2_AE;
    double a3_X, a3_AE;
    double a4_X, a4_AE;
    
    double psi_X_Fstat,  A_X_Fstat,  iota_X_Fstat,  phase_X_Fstat;
    double psi_AE_Fstat, A_AE_Fstat, iota_AE_Fstat, phase_AE_Fstat;
    
    double f0;
    double fdot;
    double fddot;
    long   q;
    double theta, phi;
    
};

void initialize_XLS(long M, double *XLS, double *AA, double *EE);

void get_filters(struct Orbit *orbit, struct Data *data, int filter_id, struct Filter *F_filter);

void get_N(struct Data *data, struct Filter *F_filter);


void get_M(struct Filter *F_filter, double **M_inv_X, double **M_inv_AE, struct Data *data);

void calc_Fstat_logL(struct Filter *F_filter);

void calc_a_i(struct Filter *F_filter);

void init_A_filters(struct Orbit *orbit, struct Data *data, struct Filter *F_filter);

void init_M_matrix(struct Filter *F_filter, struct Data *data);

void free_Filter(struct Filter *F_filter);


void get_F_params(struct Filter *F_filter);

int sgn(double v);

void get_Fstat_logL(struct Orbit *orbit, struct Data *data, double f0, double fdot, double theta, double phi, double *logL_X, double *logL_AE, double *Fparams);

void get_Fstat_xmax(struct Orbit *orbit, struct Data *data, double *x, double *xmax);

double get_Fstat_logL_wavelet(struct Orbit *orbit, struct Data *data, double f0, double fdot, double theta, double phi);

#endif /* ucb_fstatistic_h */
