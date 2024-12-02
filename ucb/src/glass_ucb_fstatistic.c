/*
 *  Copyright (C) 2019 Travis Robson, Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <glass_utils.h>

#include "glass_ucb_model.h"
#include "glass_ucb_waveform.h"
#include "glass_ucb_fstatistic.h"


void initialize_XLS(long M, double *XLS, double *AA, double *EE)
{
    for (int i=0; i<2*M; i++) // initialize arrays to zeros
    {
        XLS[i] = 0.0;
        AA[i]  = 0.0;
        EE[i]  = 0.0;
    }
}

void init_A_filters(struct Orbit *orbit, struct Data *data, struct Filter *F_filter)
{
    long M = F_filter->M_filter;
    
    F_filter->A1_fX = calloc((2*M),sizeof(double));
    F_filter->A2_fX = calloc((2*M),sizeof(double));
    F_filter->A3_fX = calloc((2*M),sizeof(double));
    F_filter->A4_fX = calloc((2*M),sizeof(double));

    F_filter->A1_fY = calloc((2*M),sizeof(double));
    F_filter->A2_fY = calloc((2*M),sizeof(double));
    F_filter->A3_fY = calloc((2*M),sizeof(double));
    F_filter->A4_fY = calloc((2*M),sizeof(double));

    F_filter->A1_fZ = calloc((2*M),sizeof(double));
    F_filter->A2_fZ = calloc((2*M),sizeof(double));
    F_filter->A3_fZ = calloc((2*M),sizeof(double));
    F_filter->A4_fZ = calloc((2*M),sizeof(double));

    F_filter->A1_fA = calloc((2*M),sizeof(double));
    F_filter->A2_fA = calloc((2*M),sizeof(double));
    F_filter->A3_fA = calloc((2*M),sizeof(double));
    F_filter->A4_fA = calloc((2*M),sizeof(double));
    
    F_filter->A1_fE = calloc((2*M),sizeof(double));
    F_filter->A2_fE = calloc((2*M),sizeof(double));
    F_filter->A3_fE = calloc((2*M),sizeof(double));
    F_filter->A4_fE = calloc((2*M),sizeof(double));
    
    initialize_XLS(M, F_filter->A1_fX, F_filter->A1_fA, F_filter->A1_fE);
    initialize_XLS(M, F_filter->A2_fX, F_filter->A2_fA, F_filter->A2_fE);
    initialize_XLS(M, F_filter->A3_fX, F_filter->A3_fA, F_filter->A3_fE);
    initialize_XLS(M, F_filter->A4_fX, F_filter->A4_fA, F_filter->A4_fE);
    
       get_filters(orbit, data, 1, F_filter);
       get_filters(orbit, data, 2, F_filter);
       get_filters(orbit, data, 3, F_filter);
       get_filters(orbit, data, 4, F_filter);
    
    // Make use of a phase shift to quickly generate other filters
//    get_filters(orbit, data, 3, F_filter);    
//    // copy  F_filter->A3_fX into F_filter->A1_fX
//    for (int i=0; i<M; i++)
//    {
//        F_filter->A1_fX[2*i+1]   = -F_filter->A3_fX[2*i];
//        F_filter->A1_fX[2*i] =  F_filter->A3_fX[2*i+1];
//        
//        F_filter->A1_fA[2*i+1]   = -F_filter->A3_fA[2*i];
//        F_filter->A1_fA[2*i] =  F_filter->A3_fA[2*i+1];
//        
//        F_filter->A1_fE[2*i+1]   = -F_filter->A3_fE[2*i];
//        F_filter->A1_fE[2*i] =  F_filter->A3_fE[2*i+1];
//    }
//    
//    get_filters(orbit, data, 4, F_filter);
//    for (int i=0; i<M; i++)
//    {
//        F_filter->A2_fX[2*i+1]    = -F_filter->A4_fX[2*i];
//        F_filter->A2_fX[2*i] =  F_filter->A4_fX[2*i+1];
//        
//        F_filter->A2_fA[2*i+1]    = -F_filter->A4_fA[2*i];
//        F_filter->A2_fA[2*i] =  F_filter->A4_fA[2*i+1];
//        
//        F_filter->A2_fE[2*i+1]    = -F_filter->A4_fE[2*i];
//        F_filter->A2_fE[2*i] =  F_filter->A4_fE[2*i+1];
//    }
//    
//    char filename[128];
//    FILE *fptr;
//    
//    sprintf(filename,"data_dft/filter_0.dat");
//    fptr=fopen(filename,"w");
//    for(int n=0; n<data->NFFT; n++)
//    {
//        int re = 2*n;
//        int im = re+1;
//        double f = data->fmin + (double)n/data->T;
//
//        fprintf(fptr,"%.12g ",f);
//        fprintf(fptr,"%.12g ",F_filter->A1_fA[re]);
//        fprintf(fptr,"%.12g ",F_filter->A1_fA[im]);
//        fprintf(fptr,"%.12g ",F_filter->A1_fE[re]);
//        fprintf(fptr,"%.12g\n",F_filter->A1_fE[im]);
//    }
//    fclose(fptr);
//
//    
//    sprintf(filename,"data_dft/filter_1.dat");
//    fptr=fopen(filename,"w");
//    for(int n=0; n<data->NFFT; n++)
//    {
//        int re = 2*n;
//        int im = re+1;
//        double f = data->fmin + (double)n/data->T;
//
//        fprintf(fptr,"%.12g ",f);
//        fprintf(fptr,"%.12g ",F_filter->A2_fA[re]);
//        fprintf(fptr,"%.12g ",F_filter->A2_fA[im]);
//        fprintf(fptr,"%.12g ",F_filter->A2_fE[re]);
//        fprintf(fptr,"%.12g\n",F_filter->A2_fE[im]);
//    }
//    fclose(fptr);
//
//    
//    sprintf(filename,"data_dft/filter_2.dat");
//    fptr=fopen(filename,"w");
//    for(int n=0; n<data->NFFT; n++)
//    {
//        int re = 2*n;
//        int im = re+1;
//        double f = data->fmin + (double)n/data->T;
//
//        fprintf(fptr,"%.12g ",f);
//        fprintf(fptr,"%.12g ",F_filter->A3_fA[re]);
//        fprintf(fptr,"%.12g ",F_filter->A3_fA[im]);
//        fprintf(fptr,"%.12g ",F_filter->A3_fE[re]);
//        fprintf(fptr,"%.12g\n",F_filter->A3_fE[im]);
//    }
//    fclose(fptr);
//
//    
//    sprintf(filename,"data_dft/filter_3.dat");
//    fptr=fopen(filename,"w");
//    for(int n=0; n<data->NFFT; n++)
//    {
//        int re = 2*n;
//        int im = re+1;
//        double f = data->fmin + (double)n/data->T;
//
//        fprintf(fptr,"%.12g ",f);
//        fprintf(fptr,"%.12g ",F_filter->A4_fA[re]);
//        fprintf(fptr,"%.12g ",F_filter->A4_fA[im]);
//        fprintf(fptr,"%.12g ",F_filter->A4_fE[re]);
//        fprintf(fptr,"%.12g\n",F_filter->A4_fE[im]);
//    }
//    fclose(fptr);
//
    //////////////////////
}

void init_M_matrix(struct Filter *F_filter, struct Data *data)
{
    long i,j;
    
    F_filter->M_inv_X  = malloc(4*sizeof(double *));//dmatrix(0,3,0,3);
    F_filter->M_inv_AE = malloc(4*sizeof(double *));//dmatrix(0,3,0,3);
    for(i=0; i<4; i++)
    {
        F_filter->M_inv_X[i]  = calloc(4,sizeof(double));//dmatrix(0,3,0,3);
        F_filter->M_inv_AE[i] = calloc(4,sizeof(double));//dmatrix(0,3,0,3);
        
        for (j=0;j<4;j++)
        {
            F_filter->M_inv_X[i][j]  = 0.0;
            F_filter->M_inv_AE[i][j] = 0.0;
        }
    }
    
    get_M(F_filter, F_filter->M_inv_X, F_filter->M_inv_AE, data);
    
}

void free_Filter(struct Filter *F_filter)
{
    free(F_filter->A1_fX);
    free(F_filter->A2_fX);
    free(F_filter->A3_fX);
    free(F_filter->A4_fX);

    free(F_filter->A1_fY);
    free(F_filter->A2_fY);
    free(F_filter->A3_fY);
    free(F_filter->A4_fY);

    free(F_filter->A1_fZ);
    free(F_filter->A2_fZ);
    free(F_filter->A3_fZ);
    free(F_filter->A4_fZ);

    free(F_filter->A1_fA);
    free(F_filter->A2_fA);
    free(F_filter->A3_fA);
    free(F_filter->A4_fA);
    
    free(F_filter->A1_fE);
    free(F_filter->A2_fE);
    free(F_filter->A3_fE);
    free(F_filter->A4_fE);
    
    for(int i=0; i<4; i++)
    {
        free(F_filter->M_inv_X[i]);
        free(F_filter->M_inv_AE[i]);
    }
    free(F_filter->M_inv_X);
    free(F_filter->M_inv_AE);
    free(F_filter);
}

void get_filters(struct Orbit *orbit, struct Data *data, int filter_id, struct Filter *F_filter)
{
    int i;
    int d=9;
    double A_f,iota_f,psi_f,phase_f;
    double *params;
    long M_filter;
    
    M_filter = F_filter->M_filter;
    
    params = calloc(d,sizeof(double));  // allocate memory for filter parameters
    for (i=0;i<d;i++)  		 // initialize the array to zeros
    {
        params[i] = 0.0;
    }
    
    // set the parameters
    params[0] = F_filter->f0;				 // f0
    params[1] = F_filter->theta;    	     // theta
    params[2] = F_filter->phi; 			     // phi
    
    params[7] = F_filter->fdot;				 // fdot
    params[8] = F_filter->fddot;			 // fddot
    
    // map to conventions for waveform generator
    params[0]*=data->T;
    params[1]=cos(params[1]);
    params[7]*=(data->T*data->T);
    params[8]*=(data->T*data->T*data->T);
    
    // Now need to calculate the Filters A^{i} (Cornish & Crowder '05)
    // A^{1} -> iota = pi/2, psi = 0,    A = 2, phase = 0
    // A^{2} -> iota = pi/2, psi = pi/4, A = 2, phase = pi
    // A^{3} -> iota = pi/2, psi = 0,    A = 2, phase = 3*pi/2
    // A^{4} -> iota = pi/2, psi = pi/4, A = 2, phase = pi/2
    
    
    // determine which parameters to pass into FAST_LISA
    if (filter_id == 1)
    {
        A_f      = 2.0;
        iota_f   = M_PI_2;
        psi_f    = 0.0;
        phase_f  = 0.0;
        
        params[3] = A_f;
        params[4] = iota_f;
        params[5] = psi_f;
        params[6] = phase_f;
        
        
        
        //FAST_LISA(params, N_filter, M_filter, F_filter->A1_fX, F_filter->A1_fA, F_filter->A1_fE);
        // map to conventions for waveform generator
        params[3]=log(params[3]);
        params[4]=cos(params[4]);
        ucb_waveform(orbit, data->format, data->T, data->t0, params, UCB_MODEL_NP, F_filter->A1_fX, F_filter->A1_fY, F_filter->A1_fZ, F_filter->A1_fA, F_filter->A1_fE, (int)M_filter, 2);
        
    } else if (filter_id == 2){
        
        A_f      = 2.0;
        iota_f   = M_PI_2;
        psi_f    = M_PI_4;
        phase_f  = M_PI;
        
        params[3] = A_f;
        params[4] = iota_f;
        params[5] = psi_f;
        params[6] = phase_f;
        
        
        //FAST_LISA(params, N_filter, M_filter, F_filter->A1_fX, F_filter->A1_fA, F_filter->A1_fE);
        // map to conventions for waveform generator
        params[3]=log(params[3]);
        params[4]=cos(params[4]);
        ucb_waveform(orbit, data->format, data->T, data->t0, params, UCB_MODEL_NP, F_filter->A2_fX, F_filter->A2_fY, F_filter->A2_fZ, F_filter->A2_fA, F_filter->A2_fE, (int)M_filter, 2);
        
    } else if (filter_id == 3){
        
        A_f      = 2.0;
        iota_f   = M_PI_2;
        psi_f    = 0.0;
        phase_f  = 3.0*M_PI_2;
        
        params[3] = A_f;
        params[4] = iota_f;
        params[5] = psi_f;
        params[6] = phase_f;
        
        
        //FAST_LISA(params, N_filter, M_filter, F_filter->A1_fX, F_filter->A1_fA, F_filter->A1_fE);
        // map to conventions for waveform generator
        params[3]=log(params[3]);
        params[4]=cos(params[4]);
        ucb_waveform(orbit, data->format, data->T, data->t0, params, UCB_MODEL_NP, F_filter->A3_fX, F_filter->A3_fY, F_filter->A3_fZ, F_filter->A3_fA, F_filter->A3_fE, (int)M_filter, 2);
        
    } else {
        
        A_f      = 2.0;
        iota_f   = M_PI_2;
        psi_f    = M_PI_4;
        phase_f  = M_PI_2;
        
        params[3] = A_f;
        params[4] = iota_f;
        params[5] = psi_f;
        params[6] = phase_f;
        
        
        //FAST_LISA(params, N_filter, M_filter, F_filter->A1_fX, F_filter->A1_fA, F_filter->A1_fE);
        // map to conventions for waveform generator
        params[3]=log(params[3]);
        params[4]=cos(params[4]);
        ucb_waveform(orbit, data->format, data->T, data->t0, params, 9, F_filter->A4_fX, F_filter->A4_fY, F_filter->A4_fZ, F_filter->A4_fA, F_filter->A4_fE, (int)M_filter, 2);
    }
    
    free(params);
}

void get_N(struct Data *data, struct Filter *F_filter)
{
    long     i,k;
    long     M_filter;
    long q;
    double *XfLS, *AALS, *EELS;
    
    /////////
    //
    // Calculate N^{i} (Cornish & Crowder '05)
    // N^{i} = (s|A^{i})
    //
    /////////
    XfLS = data->tdi->X;
    AALS = data->tdi->A;
    EELS = data->tdi->E;
    
    q  = F_filter->q - data->qmin;
    
    M_filter = F_filter->M_filter;
    
    F_filter->N1_X = 0.; F_filter->N1_AE = 0.;
    F_filter->N2_X = 0.; F_filter->N2_AE = 0.;
    F_filter->N3_X = 0.; F_filter->N3_AE = 0.;
    F_filter->N4_X = 0.; F_filter->N4_AE = 0.;
    
    for (i=0; i<M_filter; i++)
    {
        k = (q + i - M_filter/2);
        
        if(k>0 && k<data->NFFT)
        {
            
            F_filter->N1_X  += (XfLS[2*k]*F_filter->A1_fX[2*i] + XfLS[2*k+1]*F_filter->A1_fX[2*i+1])/data->noise->C[0][0][k];
            F_filter->N2_X  += (XfLS[2*k]*F_filter->A2_fX[2*i] + XfLS[2*k+1]*F_filter->A2_fX[2*i+1])/data->noise->C[0][0][k];
            F_filter->N3_X  += (XfLS[2*k]*F_filter->A3_fX[2*i] + XfLS[2*k+1]*F_filter->A3_fX[2*i+1])/data->noise->C[0][0][k];
            F_filter->N4_X  += (XfLS[2*k]*F_filter->A4_fX[2*i] + XfLS[2*k+1]*F_filter->A4_fX[2*i+1])/data->noise->C[0][0][k];
            
            F_filter->N1_AE += (AALS[2*k]*F_filter->A1_fA[2*i] + AALS[2*k+1]*F_filter->A1_fA[2*i+1]
                                +EELS[2*k]*F_filter->A1_fE[2*i] + EELS[2*k+1]*F_filter->A1_fE[2*i+1])/data->noise->C[0][0][k];
            F_filter->N2_AE += (AALS[2*k]*F_filter->A2_fA[2*i] + AALS[2*k+1]*F_filter->A2_fA[2*i+1]
                                +EELS[2*k]*F_filter->A2_fE[2*i] + EELS[2*k+1]*F_filter->A2_fE[2*i+1])/data->noise->C[0][0][k];
            F_filter->N3_AE += (AALS[2*k]*F_filter->A3_fA[2*i] + AALS[2*k+1]*F_filter->A3_fA[2*i+1]
                                +EELS[2*k]*F_filter->A3_fE[2*i] + EELS[2*k+1]*F_filter->A3_fE[2*i+1])/data->noise->C[0][0][k];
            F_filter->N4_AE += (AALS[2*k]*F_filter->A4_fA[2*i] + AALS[2*k+1]*F_filter->A4_fA[2*i+1]
                                +EELS[2*k]*F_filter->A4_fE[2*i] + EELS[2*k+1]*F_filter->A4_fE[2*i+1])/data->noise->C[0][0][k];
        }
    }
    
    F_filter->N1_X  *= 4.0;
    F_filter->N2_X  *= 4.0;
    F_filter->N3_X  *= 4.0;
    F_filter->N4_X  *= 4.0;
    
    F_filter->N1_AE *= 4.0;
    F_filter->N2_AE *= 4.0;
    F_filter->N3_AE *= 4.0;
    F_filter->N4_AE *= 4.0;
    
}

void get_M(struct Filter *F_filter, double **M_inv_X, double **M_inv_AE, struct Data *data)
{
    int i,j,k;
    long M;
    long q, M_filter;
    
    M = F_filter->M_filter;
    q = F_filter->q-data->qmin;
    M_filter = F_filter->M_filter;
    
    for (i=0; i<M; i++)
    {
        k = ((int)q + i - (int)M_filter/2);
        
        if(k>0 && k<data->NFFT)
        {
            M_inv_X[0][0] += 4.0*(F_filter->A1_fX[2*i]  *F_filter->A1_fX[2*i]
                                  +F_filter->A1_fX[2*i+1]*F_filter->A1_fX[2*i+1])/data->noise->C[0][0][k];
            M_inv_X[0][1] += 4.0*(F_filter->A1_fX[2*i]  *F_filter->A2_fX[2*i]
                                  +F_filter->A1_fX[2*i+1]*F_filter->A2_fX[2*i+1])/data->noise->C[0][0][k];
            M_inv_X[0][2] += 4.0*(F_filter->A1_fX[2*i]  *F_filter->A3_fX[2*i]
                                  +F_filter->A1_fX[2*i+1]*F_filter->A3_fX[2*i+1])/data->noise->C[0][0][k];
            M_inv_X[0][3] += 4.0*(F_filter->A1_fX[2*i]  *F_filter->A4_fX[2*i]
                                  +F_filter->A1_fX[2*i+1]*F_filter->A4_fX[2*i+1])/data->noise->C[0][0][k];
            
            M_inv_X[1][1] += 4.0*(F_filter->A2_fX[2*i]  *F_filter->A2_fX[2*i]
                                  +F_filter->A2_fX[2*i+1]*F_filter->A2_fX[2*i+1])/data->noise->C[0][0][k];
            M_inv_X[1][2] += 4.0*(F_filter->A2_fX[2*i]  *F_filter->A3_fX[2*i]
                                  +F_filter->A2_fX[2*i+1]*F_filter->A3_fX[2*i+1])/data->noise->C[0][0][k];
            M_inv_X[1][3] += 4.0*(F_filter->A2_fX[2*i]  *F_filter->A4_fX[2*i]
                                  +F_filter->A2_fX[2*i+1]*F_filter->A4_fX[2*i+1])/data->noise->C[0][0][k];
            
            M_inv_X[2][2] += 4.0*(F_filter->A3_fX[2*i]  *F_filter->A3_fX[2*i]
                                  +F_filter->A3_fX[2*i+1]*F_filter->A3_fX[2*i+1])/data->noise->C[0][0][k];
            M_inv_X[2][3] += 4.0*(F_filter->A3_fX[2*i]  *F_filter->A4_fX[2*i]
                                  +F_filter->A3_fX[2*i+1]*F_filter->A4_fX[2*i+1])/data->noise->C[0][0][k];
            
            M_inv_X[3][3] += 4.0*(F_filter->A4_fX[2*i]  *F_filter->A4_fX[2*i]
                                  +F_filter->A4_fX[2*i+1]*F_filter->A4_fX[2*i+1])/data->noise->C[0][0][k];
            
            
            
            
            M_inv_AE[0][0] += 4.0*(F_filter->A1_fA[2*i]  *F_filter->A1_fA[2*i]
                                   +F_filter->A1_fA[2*i+1]*F_filter->A1_fA[2*i+1])/data->noise->C[0][0][k]
            +4.0*(F_filter->A1_fE[2*i]  *F_filter->A1_fE[2*i]
                  +F_filter->A1_fE[2*i+1]*F_filter->A1_fE[2*i+1])/data->noise->C[0][0][k];
            M_inv_AE[0][1] += 4.0*(F_filter->A1_fA[2*i]  *F_filter->A2_fA[2*i]
                                   +F_filter->A1_fA[2*i+1]*F_filter->A2_fA[2*i+1])/data->noise->C[0][0][k]
            +4.0*(F_filter->A1_fE[2*i]  *F_filter->A2_fE[2*i]
                  +F_filter->A1_fE[2*i+1]*F_filter->A2_fE[2*i+1])/data->noise->C[0][0][k];
            M_inv_AE[0][2] += 4.0*(F_filter->A1_fA[2*i]  *F_filter->A3_fA[2*i]
                                   +F_filter->A1_fA[2*i+1]*F_filter->A3_fA[2*i+1])/data->noise->C[0][0][k]
            +4.0*(F_filter->A1_fE[2*i]  *F_filter->A3_fE[2*i]
                  +F_filter->A1_fE[2*i+1]*F_filter->A3_fE[2*i+1])/data->noise->C[0][0][k];
            M_inv_AE[0][3] += 4.0*(F_filter->A1_fA[2*i]  *F_filter->A4_fA[2*i]
                                   +F_filter->A1_fA[2*i+1]*F_filter->A4_fA[2*i+1])/data->noise->C[0][0][k]
            +4.0*(F_filter->A1_fE[2*i]  *F_filter->A4_fE[2*i]
                  +F_filter->A1_fE[2*i+1]*F_filter->A4_fE[2*i+1])/data->noise->C[0][0][k];
            
            
            M_inv_AE[1][1] += 4.0*(F_filter->A2_fA[2*i]  *F_filter->A2_fA[2*i]
                                   +F_filter->A2_fA[2*i+1]*F_filter->A2_fA[2*i+1])/data->noise->C[0][0][k]
            +4.0*(F_filter->A2_fE[2*i]  *F_filter->A2_fE[2*i]
                  +F_filter->A2_fE[2*i+1]*F_filter->A2_fE[2*i+1])/data->noise->C[0][0][k];
            M_inv_AE[1][2] += 4.0*(F_filter->A2_fA[2*i]  *F_filter->A3_fA[2*i]
                                   +F_filter->A2_fA[2*i+1]*F_filter->A3_fA[2*i+1])/data->noise->C[0][0][k]
            +4.0*(F_filter->A2_fE[2*i]  *F_filter->A3_fE[2*i]
                  +F_filter->A2_fE[2*i+1]*F_filter->A3_fE[2*i+1])/data->noise->C[0][0][k];
            M_inv_AE[1][3] += 4.0*(F_filter->A2_fA[2*i]  *F_filter->A4_fA[2*i]
                                   +F_filter->A2_fA[2*i+1]*F_filter->A4_fA[2*i+1])/data->noise->C[0][0][k]
            +4.0*(F_filter->A2_fE[2*i]  *F_filter->A4_fE[2*i]
                  +F_filter->A2_fE[2*i+1]*F_filter->A4_fE[2*i+1])/data->noise->C[0][0][k];
            
            
            M_inv_AE[2][2] += 4.0*(F_filter->A3_fA[2*i]  *F_filter->A3_fA[2*i]
                                   +F_filter->A3_fA[2*i+1]*F_filter->A3_fA[2*i+1])/data->noise->C[0][0][k]
            +4.0*(F_filter->A3_fE[2*i]  *F_filter->A3_fE[2*i]
                  +F_filter->A3_fE[2*i+1]*F_filter->A3_fE[2*i+1])/data->noise->C[0][0][k];
            M_inv_AE[2][3] += 4.0*(F_filter->A3_fA[2*i]  *F_filter->A4_fA[2*i]
                                   +F_filter->A3_fA[2*i+1]*F_filter->A4_fA[2*i+1])/data->noise->C[0][0][k]
            +4.0*(F_filter->A3_fE[2*i]  *F_filter->A4_fE[2*i]
                  +F_filter->A3_fE[2*i+1]*F_filter->A4_fE[2*i+1])/data->noise->C[0][0][k];
            
            
            M_inv_AE[3][3] += 4.0*(F_filter->A4_fA[2*i]  *F_filter->A4_fA[2*i]
                                   +F_filter->A4_fA[2*i+1]*F_filter->A4_fA[2*i+1])/data->noise->C[0][0][k]
            +4.0*(F_filter->A4_fE[2*i]  *F_filter->A4_fE[2*i]
                  +F_filter->A4_fE[2*i+1]*F_filter->A4_fE[2*i+1])/data->noise->C[0][0][k];
        }
    }
    
    // make use of symmetry to populate the other pieces
    for (j=0; j<4; j++)
    {
        for (k=j+1; k<4; k++)
        {
            M_inv_X[k][j]  = M_inv_X[j][k];
            M_inv_AE[k][j] = M_inv_AE[j][k];
        }
    }
    
    invert_matrix(M_inv_X,4);
    invert_matrix(M_inv_AE,4);
    
}


int sgn(double v) {
    if (v < 0) return -1;
    if (v > 0) return 1;
    return 0;
}

void get_F_params(struct Filter *F_filter)
{
    
    double Aplus_X,Across_X;
    double Aplus_AE,Across_AE;
    
    double psi_X_Fstat,  c,   A_X_Fstat,  iota_X_Fstat,  phase_X_Fstat;
    double psi_AE_Fstat, cAE, A_AE_Fstat, iota_AE_Fstat, phase_AE_Fstat;
    
    // TODO: redundant terms here, save time break it up
    double a1p4X = F_filter->a1_X + F_filter->a4_X;
    double a2m3X = F_filter->a2_X - F_filter->a3_X;
    double a1m4X = F_filter->a1_X - F_filter->a4_X;
    double a2p3X = F_filter->a2_X + F_filter->a3_X;
    
    Aplus_X  = sqrt( a1p4X*a1p4X + a2m3X*a2m3X  ) + sqrt( a1m4X*a1m4X + a2p3X*a2p3X );
    Across_X = sqrt( a1p4X*a1p4X + a2m3X*a2m3X  ) - sqrt( a1m4X*a1m4X + a2p3X*a2p3X );
    
    psi_X_Fstat   = 0.5*atan2(Aplus_X*F_filter->a4_X - Across_X*F_filter->a1_X,
                              -(Across_X*F_filter->a2_X + Aplus_X*F_filter->a3_X) );
    
    c             = (double)sgn(sin(2.0*psi_X_Fstat));
    A_X_Fstat     = 0.5*(   Aplus_X + sqrt(Aplus_X*Aplus_X - Across_X*Across_X)  );
    iota_X_Fstat  = acos(  -Across_X/(Aplus_X + sqrt(Aplus_X*Aplus_X - Across_X*Across_X))  );
    phase_X_Fstat = atan2(  c*(Aplus_X*F_filter->a4_X - Across_X*F_filter->a1_X),
                          -c*(Across_X*F_filter->a3_X + Aplus_X*F_filter->a2_X)  );
    
    double a1p4AE = F_filter->a1_AE + F_filter->a4_AE;
    double a2m3AE = F_filter->a2_AE - F_filter->a3_AE;
    double a1m4AE = F_filter->a1_AE - F_filter->a4_AE;
    double a2p3AE = F_filter->a2_AE + F_filter->a3_AE;
    
    Aplus_AE  = sqrt( a1p4AE*a1p4AE + a2m3AE*a2m3AE  ) + sqrt( a1m4AE*a1m4AE + a2p3AE*a2p3AE );
    Across_AE = sqrt( a1p4AE*a1p4AE + a2m3AE*a2m3AE  ) - sqrt( a1m4AE*a1m4AE + a2p3AE*a2p3AE );
    
    
    psi_AE_Fstat   = atan2(Aplus_AE*F_filter->a4_AE - Across_AE*F_filter->a1_AE,
                               -(Across_AE*F_filter->a2_AE + Aplus_AE*F_filter->a3_AE) );
    
    psi_AE_Fstat *= 0.5;
    
    cAE            = (double)sgn(sin(2.0*psi_AE_Fstat));
    A_AE_Fstat     = 0.5*(  Aplus_AE + sqrt(Aplus_AE*Aplus_AE - Across_AE*Across_AE)  );
    iota_AE_Fstat  = acos(  -Across_AE/(Aplus_AE + sqrt(Aplus_AE*Aplus_AE - Across_AE*Across_AE))  );
    phase_AE_Fstat = atan2(  cAE*(Aplus_AE*F_filter->a4_AE - Across_AE*F_filter->a1_AE),
                           -cAE*(Across_AE*F_filter->a3_AE + Aplus_AE*F_filter->a2_AE)  );

    // 	printf("\nX F-stat ML Values                         AE F-stat ML Values\n");
    // 	printf("-------------------------------            -------------------------------\n");
    // 	printf("psi X:   	%e               psi AE:   	%e\n",  tan(psi_X_Fstat),   tan(psi_AE_Fstat));
    // 	printf("A X:     	%e               A AE:     	%e\n",  A_X_Fstat,          A_AE_Fstat);
    // 	printf("iota X:  	%e               iota AE:  	%e \n", cos(iota_X_Fstat),  cos(iota_AE_Fstat));
    // 	printf("phase X: 	%e               phase AE: 	%e\n",  tan(phase_X_Fstat), tan(phase_AE_Fstat));
    //
    // 	printf("\nDifference                                 Difference\n");
    // 	printf("-------------------------------            -------------------------------\n");
    // 	printf("psi X:   	%e               psi AE:   	%e\n", tan(psi_X_Fstat)-tan(src->params[5]),   tan(psi_AE_Fstat)-tan(src->params[5]));
    // 	printf("A X:     	%e               A AE:     	%e\n", A_X_Fstat - src->params[3],             A_AE_Fstat - src->params[3]);
    // 	printf("iota X:  	%e               iota AE:  	%e\n", cos(iota_X_Fstat)-cos(src->params[4]),  cos(iota_AE_Fstat)-cos(src->params[4]));
    // 	printf("phase X: 	%e              phase AE: 	%e\n", tan(phase_X_Fstat)-tan(src->params[6]), tan(phase_AE_Fstat)-tan(src->params[6]));
    //
    // 	printf("\n================================================================================\n\n");
    
    F_filter->psi_X_Fstat   = psi_X_Fstat;
    F_filter->A_X_Fstat     = A_X_Fstat;
    F_filter->iota_X_Fstat  = iota_X_Fstat;
    F_filter->phase_X_Fstat = phase_X_Fstat;
    
    F_filter->psi_AE_Fstat   = psi_AE_Fstat;
    F_filter->A_AE_Fstat     = A_AE_Fstat;
    F_filter->iota_AE_Fstat  = iota_AE_Fstat;
    F_filter->phase_AE_Fstat = phase_AE_Fstat;
}

void calc_Fstat_logL(struct Filter *F_filter)
{
    F_filter->Fstat_X  = 0.5*(F_filter->M_inv_X[0][0]*F_filter->N1_X*F_filter->N1_X
                              + 2*F_filter->M_inv_X[0][1]*F_filter->N1_X*F_filter->N2_X
                              + 2*F_filter->M_inv_X[0][2]*F_filter->N1_X*F_filter->N3_X
                              + 2*F_filter->M_inv_X[0][3]*F_filter->N1_X*F_filter->N4_X
                              +   F_filter->M_inv_X[1][1]*F_filter->N2_X*F_filter->N2_X
                              + 2*F_filter->M_inv_X[1][2]*F_filter->N2_X*F_filter->N3_X
                              + 2*F_filter->M_inv_X[1][3]*F_filter->N2_X*F_filter->N4_X
                              +   F_filter->M_inv_X[2][2]*F_filter->N3_X*F_filter->N3_X
                              + 2*F_filter->M_inv_X[2][3]*F_filter->N3_X*F_filter->N4_X
                              +   F_filter->M_inv_X[3][3]*F_filter->N4_X*F_filter->N4_X);
    
    F_filter->Fstat_AE = 0.5*(F_filter->M_inv_AE[0][0]*F_filter->N1_AE*F_filter->N1_AE
                              + 2*F_filter->M_inv_AE[0][1]*F_filter->N1_AE*F_filter->N2_AE
                              + 2*F_filter->M_inv_AE[0][2]*F_filter->N1_AE*F_filter->N3_AE
                              + 2*F_filter->M_inv_AE[0][3]*F_filter->N1_AE*F_filter->N4_AE
                              +   F_filter->M_inv_AE[1][1]*F_filter->N2_AE*F_filter->N2_AE
                              + 2*F_filter->M_inv_AE[1][2]*F_filter->N2_AE*F_filter->N3_AE
                              + 2*F_filter->M_inv_AE[1][3]*F_filter->N2_AE*F_filter->N4_AE
                              +   F_filter->M_inv_AE[2][2]*F_filter->N3_AE*F_filter->N3_AE
                              + 2*F_filter->M_inv_AE[2][3]*F_filter->N3_AE*F_filter->N4_AE
                              +   F_filter->M_inv_AE[3][3]*F_filter->N4_AE*F_filter->N4_AE);
}

void calc_a_i(struct Filter *F_filter)
{
    
    F_filter->a1_X = F_filter->M_inv_X[0][0]*F_filter->N1_X + F_filter->M_inv_X[0][1]*F_filter->N2_X
    +F_filter->M_inv_X[0][2]*F_filter->N3_X + F_filter->M_inv_X[0][3]*F_filter->N4_X;
    F_filter->a2_X = F_filter->M_inv_X[1][0]*F_filter->N1_X + F_filter->M_inv_X[1][1]*F_filter->N2_X
    +F_filter->M_inv_X[1][2]*F_filter->N3_X + F_filter->M_inv_X[1][3]*F_filter->N4_X;
    F_filter->a3_X = F_filter->M_inv_X[2][0]*F_filter->N1_X + F_filter->M_inv_X[2][1]*F_filter->N2_X
    +F_filter->M_inv_X[2][2]*F_filter->N3_X + F_filter->M_inv_X[2][3]*F_filter->N4_X;
    F_filter->a4_X = F_filter->M_inv_X[3][0]*F_filter->N1_X + F_filter->M_inv_X[3][1]*F_filter->N2_X
    +F_filter->M_inv_X[3][2]*F_filter->N3_X + F_filter->M_inv_X[3][3]*F_filter->N4_X;
    
    F_filter->a1_AE = F_filter->M_inv_AE[0][0]*F_filter->N1_AE + F_filter->M_inv_AE[0][1]*F_filter->N2_AE
    +F_filter->M_inv_AE[0][2]*F_filter->N3_AE + F_filter->M_inv_AE[0][3]*F_filter->N4_AE;
    F_filter->a2_AE = F_filter->M_inv_AE[1][0]*F_filter->N1_AE + F_filter->M_inv_AE[1][1]*F_filter->N2_AE
    +F_filter->M_inv_AE[1][2]*F_filter->N3_AE + F_filter->M_inv_AE[1][3]*F_filter->N4_AE;
    F_filter->a3_AE = F_filter->M_inv_AE[2][0]*F_filter->N1_AE + F_filter->M_inv_AE[2][1]*F_filter->N2_AE
    +F_filter->M_inv_AE[2][2]*F_filter->N3_AE + F_filter->M_inv_AE[2][3]*F_filter->N4_AE;
    F_filter->a4_AE = F_filter->M_inv_AE[3][0]*F_filter->N1_AE + F_filter->M_inv_AE[3][1]*F_filter->N2_AE
    +F_filter->M_inv_AE[3][2]*F_filter->N3_AE + F_filter->M_inv_AE[3][3]*F_filter->N4_AE;
}

void get_Fstat_logL(struct Orbit *orbit, struct Data *data, double f0, double fdot, double theta, double phi, double *logL_X, double *logL_AE, double *Fparams)
{
    long M_filter, N_filter;
    long q;
    
    q = (long)(f0*data->T); 	// carrier frequency bin
    
    struct Filter *F_filter = malloc(sizeof(struct Filter));
    
    
    F_filter->f0     = f0;
    F_filter->fdot   = fdot;
    F_filter->fddot  = 0.;    //11.0/3.0*fdot*fdot/f0 
    F_filter->q      = q;
    F_filter->theta  = theta;
    F_filter->phi    = phi;
    

    int BW = ucb_bandwidth(orbit->L, orbit->fstar, f0, fdot, cos(theta), 1.e-22, data->T, data->NFFT);
    M_filter = BW;
    N_filter = BW;

    F_filter->M_filter = M_filter;
    F_filter->N_filter = N_filter;
    init_A_filters(orbit, data, F_filter);
    
    /////////
    //
    // Calculate N^{i} (Cornish & Crowder '05)
    // N^{i} = (s|A^{i})
    //
    /////////	
    get_N(data, F_filter);
    
    /////////
    //
    // Calculate M^{ij} inverse (Cornish & Crowder '05)
    // M^{ij} = (A^{i}|A^{j})
    // It is symmetric
    //
    /////////
    init_M_matrix(F_filter, data);	
    
    calc_Fstat_logL(F_filter);
    
    *logL_X  = F_filter->Fstat_X;
    *logL_AE = F_filter->Fstat_AE;
    
    calc_a_i(F_filter);           // calculate the a_{i}'s associated with filters to get F-stat ML params 
    get_F_params(F_filter);  // get the F-stat ML parameters	
    
    Fparams[0] = F_filter->A_AE_Fstat;
    Fparams[1] = F_filter->iota_AE_Fstat;
    Fparams[2] = F_filter->psi_AE_Fstat;
    Fparams[3] = F_filter->phase_AE_Fstat;
    
    free_Filter(F_filter);
}

void get_Fstat_xmax(struct Orbit *orbit, struct Data *data, double *x, double *xmax)
{
    struct Filter *F_filter = malloc(sizeof(struct Filter));
        
    double f0 = x[0]/data->T;
    double theta = acos(x[1]);
    double phi = x[2];
    double fdot = x[7]/(data->T*data->T);

    F_filter->f0     = f0;
    F_filter->fdot   = fdot;
    F_filter->fddot  = 0.;    //11.0/3.0*fdot*fdot/f0
    F_filter->q      = (long)x[0];
    F_filter->theta  = theta;
    F_filter->phi    = phi;
    
 
    long M_filter, N_filter;
    
    int BW = ucb_bandwidth(orbit->L, orbit->fstar, f0, fdot, cos(theta), 1.e-22, data->T, data->NFFT);
    M_filter = BW;
    N_filter = BW;

    F_filter->M_filter = M_filter;
    F_filter->N_filter = N_filter;

    
    init_A_filters(orbit, data, F_filter);
    
    /////////
    //
    // Calculate N^{i} (Cornish & Crowder '05)
    // N^{i} = (s|A^{i})
    //
    /////////
    get_N(data, F_filter);
    
    /////////
    //
    // Calculate M^{ij} inverse (Cornish & Crowder '05)
    // M^{ij} = (A^{i}|A^{j})
    // It is symmetric
    //
    /////////
    init_M_matrix(F_filter, data);
    
    calc_a_i(F_filter);           // calculate the a_{i}'s associated with filters to get F-stat ML params
    get_F_params(F_filter);  // get the F-stat ML parameters
    
    xmax[3] = log(F_filter->A_AE_Fstat);
    xmax[4] = cos(F_filter->iota_AE_Fstat);
    xmax[5] = F_filter->psi_AE_Fstat;
    xmax[6] = F_filter->phase_AE_Fstat;
    
    free_Filter(F_filter);
}


double get_Fstat_logL_wavelet(struct Orbit *orbit, struct Data *data, double f0, double fdot, double theta, double phi)
{

    /* compute filters A_i */
    
    int Nfilter = 4;
    
    struct Source **A = malloc(Nfilter*(sizeof(struct Source *)));
    for(int i=0; i<Nfilter; i++) 
    {
        A[i] = malloc(sizeof(struct Source));
        alloc_source(A[i],data->N,data->Nchannel);
    }
    
    //set parameters for each filter
    for(int i=0; i<Nfilter; i++)
    {
        A[i]->params[0] = f0*data->T;  //frequency
        A[i]->params[1] = cos(theta);  //ecliptic lat
        A[i]->params[2] = phi;         //ecliptic lon
        A[i]->params[3] = log(2.0);    //log amplitude
        A[i]->params[4] = cos(M_PI_2); //cos inclination

        switch(i)
        {
            case 0:
                A[i]->params[5] = 0;      //polarization
                A[i]->params[6] = 0;      //initial phase
                break;
            case 1:
                A[i]->params[5] = M_PI_4; //polarization
                A[i]->params[6] = M_PI;   //initial phase
                break;
            case 2:
                A[i]->params[5] = 0;       //polarization
                A[i]->params[6] = 3*M_PI_2;//initial phase
                break;
            case 3:
                A[i]->params[5] = M_PI_4;  //polarization
                A[i]->params[6] = M_PI_2;  //initial phase
                break;
            default:
                break;
        }

        A[i]->params[7] = fdot*data->T*data->T;
        if(UCB_MODEL_NP>8) A[i]->params[8] = 0.0;

        ucb_waveform_wavelet(orbit,data->wdm,data->T, data->t0, A[i]->params, A[i]->list, &A[i]->Nlist, A[i]->tdi->X, A[i]->tdi->Y, A[i]->tdi->Z);

        //catch waveforms that are out of band
        if(A[i]->Nlist==0) return 1.0;
        
    }
    
    /* compute vectors N_i = (s|A_i) */
    double *N = double_vector(Nfilter);

    for(int i=0; i<Nfilter; i++)
    {
        N[i] += wavelet_nwip(data->tdi->X, A[i]->tdi->X, data->noise->invC[0][0], A[i]->list, A[i]->Nlist);
        N[i] += wavelet_nwip(data->tdi->Y, A[i]->tdi->Y, data->noise->invC[1][1], A[i]->list, A[i]->Nlist);
        N[i] += wavelet_nwip(data->tdi->Z, A[i]->tdi->Z, data->noise->invC[2][2], A[i]->list, A[i]->Nlist);
        N[i] += wavelet_nwip(data->tdi->X, A[i]->tdi->Y, data->noise->invC[0][1], A[i]->list, A[i]->Nlist);
        N[i] += wavelet_nwip(data->tdi->X, A[i]->tdi->Z, data->noise->invC[0][2], A[i]->list, A[i]->Nlist);
        N[i] += wavelet_nwip(data->tdi->Y, A[i]->tdi->Z, data->noise->invC[1][2], A[i]->list, A[i]->Nlist);
        N[i] += wavelet_nwip(data->tdi->Y, A[i]->tdi->X, data->noise->invC[1][0], A[i]->list, A[i]->Nlist);
        N[i] += wavelet_nwip(data->tdi->Z, A[i]->tdi->X, data->noise->invC[2][0], A[i]->list, A[i]->Nlist);
        N[i] += wavelet_nwip(data->tdi->Z, A[i]->tdi->Y, data->noise->invC[2][1], A[i]->list, A[i]->Nlist);
    }

    /* compute matrix M_ij = (A_i|A_j) */
    double **M = double_matrix(Nfilter,Nfilter);
    for(int i=0; i<Nfilter; i++)
    {
        for(int j=0; j<Nfilter; j++) //Mij symmetric so this could be faster
        {
            int *list = int_vector(data->N);
            int Nlist;
            list_union(A[i]->list, A[j]->list, A[i]->Nlist, A[j]->Nlist, list, &Nlist);

            M[i][j] += wavelet_nwip(A[i]->tdi->X, A[j]->tdi->X, data->noise->invC[0][0], list, Nlist);
            M[i][j] += wavelet_nwip(A[i]->tdi->Y, A[j]->tdi->Y, data->noise->invC[1][1], list, Nlist);
            M[i][j] += wavelet_nwip(A[i]->tdi->Z, A[j]->tdi->Z, data->noise->invC[2][2], list, Nlist);
            M[i][j] += wavelet_nwip(A[i]->tdi->X, A[j]->tdi->Y, data->noise->invC[0][1], list, Nlist);
            M[i][j] += wavelet_nwip(A[i]->tdi->X, A[j]->tdi->Z, data->noise->invC[0][2], list, Nlist);
            M[i][j] += wavelet_nwip(A[i]->tdi->Y, A[j]->tdi->Z, data->noise->invC[1][2], list, Nlist);
            M[i][j] += wavelet_nwip(A[i]->tdi->Y, A[j]->tdi->X, data->noise->invC[1][0], list, Nlist);
            M[i][j] += wavelet_nwip(A[i]->tdi->Z, A[j]->tdi->X, data->noise->invC[2][0], list, Nlist);
            M[i][j] += wavelet_nwip(A[i]->tdi->Z, A[j]->tdi->Y, data->noise->invC[2][1], list, Nlist);

            free(list);
        }
    }
    
    //invert M-matrix
    invert_matrix(M,Nfilter);

    /* get logL = N * M^-1 * N */
    double logL = 0.0;
    for(int i=0; i<Nfilter; i++) 
        for(int j=0; j<Nfilter; j++) 
            logL += N[i]*M[i][j]*N[j]; //also ignoring symmetry


    for(int i=0; i<Nfilter; i++) free_source(A[i]);
    free(A);

    free_double_vector(N);
    free_double_matrix(M,Nfilter);

    return logL;
}

