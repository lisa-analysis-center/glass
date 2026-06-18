/*
 * Copyright 2023 Tyson B. Littenberg
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


#include <glass_utils.h>

#include "glass_noise.h"


void print_spline_state(struct SplineModel *model, FILE *fptr, int step)
{
    fprintf(fptr,"%i %.12g %i\n",step,model->logL,model->spline->N);
}

void print_instrument_state(struct InstrumentModel *model, FILE *fptr)
{
    for(int i=0; i<model->Nlink; i++)fprintf(fptr,"%.12g ", model->sacc[i]);
    for(int i=0; i<model->Nlink; i++)fprintf(fptr,"%.12g ", model->soms[i]);
    //time-varying noise: append the per-link linear drift parameters
    if(model->drift_enabled)
    {
        for(int i=0; i<model->Nlink; i++)fprintf(fptr,"%.12g ", model->sacc_drift[i]);
        for(int i=0; i<model->Nlink; i++)fprintf(fptr,"%.12g ", model->soms_drift[i]);
    }
}

inline void print_foreground_state(struct ForegroundModel *model, FILE *fptr)
{
    for(int i=0; i<model->Nparams; i++)fprintf(fptr,"%.12g ", model->sgal[i]);
}

inline void print_sgwb_state(struct SGWBModel *model, FILE *fptr)
{
    for(int i=0; i<model->Nparams; i++)fprintf(fptr,"%.12g ", model->params[i]);
}

void print_noise_model_dynamic(struct Data *data, struct Noise *noise, char filename[])
{
    struct Wavelets* wdm = data->wdm;
    FILE *fptr = fopen(filename,"w");
    if (!fptr) {
        fprintf(stderr, "Unable to open file for writing! %s\n", filename);
        exit(-1);
    }
    int k;
    int jmin=(int)round(noise->f[0]/wdm->df);
    int jmax=(int)round(noise->f[noise->N-1]/wdm->df)+1; 
    for(int i=0; i<wdm->NT; i++)
    {
        double t = i*wdm->dt;
        for(int j=jmin; j<jmax; j++)
        {
            double f = j*wdm->df;
            wavelet_pixel_to_index(wdm,i,j,&k);
            k-=wdm->kmin;
            fprintf(fptr,"%lg ",t);
            fprintf(fptr,"%lg ",f);
            for(size_t m=0; m<noise->Nchannel; m++)
                fprintf(fptr,"%lg ",noise->C[m][m][k]);
            fprintf(fptr,"%lg ",noise->C[0][1][k]);
            fprintf(fptr,"%lg ",noise->C[0][2][k]);
            fprintf(fptr,"%lg ",noise->C[1][2][k]);
            fprintf(fptr,"\n");
        } //loop over frequency layers
    } //loop over time slices
    fclose(fptr);
}

void print_noise_model_dynamic_coarse(struct Data *data, struct Noise *coarse, int Q, char filename[])
{
    struct Wavelets* wdm = data->wdm;
    FILE *fptr = fopen(filename,"w");
    if (!fptr) {
        fprintf(stderr, "Unable to open file for writing! %s\n", filename);
        exit(-1);
    }
    int Ncoarse = wdm->NT / Q;
    int jmin=(int)round(coarse->f[0]/wdm->df);
    int jmax=(int)round(coarse->f[coarse->N-1]/wdm->df)+1;
    for(int q=0; q<Ncoarse; q++)
    {
        double t = (double)q*wdm->dt*Q;
        for(int j=jmin; j<jmax; j++)
        {
            double f = j*wdm->df;
            int k = q + (j-jmin)*Ncoarse;
            fprintf(fptr,"%lg ",t);
            fprintf(fptr,"%lg ",f);
            for(size_t m=0; m<coarse->Nchannel; m++)
                fprintf(fptr,"%lg ",coarse->C[m][m][k]);
            fprintf(fptr,"%lg ",coarse->C[0][1][k]);
            fprintf(fptr,"%lg ",coarse->C[0][2][k]);
            fprintf(fptr,"%lg ",coarse->C[1][2][k]);
            fprintf(fptr,"\n");
        }
    }
    fclose(fptr);
}

void print_noise_model(struct Noise *noise, char filename[])
{
    FILE *fptr = fopen(filename,"w");
    if (fptr==NULL) {
        fprintf(stderr, "filesystem error, couldn't open %s for writing.", filename);
        exit(-1);
    }
    for(int i=0; i<noise->N; i++)
    {
        fprintf(fptr,"%lg ",noise->f[i]);
        for(int j=0; j<noise->Nchannel; j++)
            fprintf(fptr,"%lg ",noise->C[j][j][i]);
        fprintf(fptr,"%lg ",noise->C[0][1][i]);
        fprintf(fptr,"%lg ",noise->C[0][2][i]);
        fprintf(fptr,"%lg ",noise->C[1][2][i]);
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}

void print_whitened_data(struct Data *data, struct Noise *noise, char filename[])
{
    FILE *fptr = fopen(filename,"w");
    if (!strcmp(data->basis,"wavelet"))
    {
        struct Wavelets *wdm = data->wdm;
        int k;
        for (int i=0; i<wdm->NT; i++)
        {
            double t = i*wdm->dt;
            for (int j=data->lmin; j<data->lmax; j++)
            {
                double f = j*wdm->df;
                wavelet_pixel_to_index(wdm, i, j, &k);
                k -= wdm->kmin;
                double wX = data->dwt->X[k] / sqrt(noise->C[0][0][k]);
                double wY = data->dwt->Y[k] / sqrt(noise->C[1][1][k]);
                double wZ = data->dwt->Z[k] / sqrt(noise->C[2][2][k]);
                fprintf(fptr,"%lg %lg %lg %lg %lg\n", t, f, wX, wY, wZ);
            }
            fprintf(fptr,"\n");
        }
    }
    else
    {
        for(int i=0; i<noise->N; i++)
        {
            fprintf(fptr,"%lg ",noise->f[i]);
            fprintf(fptr,"%lg %lg ",data->tdi->X[2*i]/sqrt(noise->C[0][0][i]),data->tdi->X[2*i+1]/sqrt(noise->C[0][0][i]));
            fprintf(fptr,"%lg %lg ",data->tdi->Y[2*i]/sqrt(noise->C[1][1][i]),data->tdi->Y[2*i+1]/sqrt(noise->C[1][1][i]));
            fprintf(fptr,"%lg %lg ",data->tdi->Z[2*i]/sqrt(noise->C[2][2][i]),data->tdi->Z[2*i+1]/sqrt(noise->C[2][2][i]));
            fprintf(fptr,"\n");
        }
    }
    fclose(fptr);
}

void print_noise_reconstruction(struct Data *data, struct Flags *flags)
{
    FILE *fptr_Snf;
    char filename[MAXSTRINGSIZE];
    
    sprintf(filename,"%s/power_noise_reconstruction.dat",data->dataDir);
    fptr_Snf=fopen(filename,"w");
    
    for(int i=0; i<data->NFFT; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        fprintf(fptr_Snf,"%.12g ",f);
        
        for(int j=0; j<data->Nchannel; j++)
        {
            double_sort(data->S_pow[i][j],data->Nwave);

            double S_med   = get_quantile_from_sorted_data(data->S_pow[i][j], data->Nwave, 0.50);
            double S_lo_50 = get_quantile_from_sorted_data(data->S_pow[i][j], data->Nwave, 0.25);
            double S_hi_50 = get_quantile_from_sorted_data(data->S_pow[i][j], data->Nwave, 0.75);
            double S_lo_90 = get_quantile_from_sorted_data(data->S_pow[i][j], data->Nwave, 0.05);
            double S_hi_90 = get_quantile_from_sorted_data(data->S_pow[i][j], data->Nwave, 0.95);
            
            fprintf(fptr_Snf,"%lg ",S_med);
            fprintf(fptr_Snf,"%lg ",S_lo_50);
            fprintf(fptr_Snf,"%lg ",S_hi_50);
            fprintf(fptr_Snf,"%lg ",S_lo_90);
            fprintf(fptr_Snf,"%lg ",S_hi_90);
        }
        fprintf(fptr_Snf,"\n");
    }
    fclose(fptr_Snf);
    
}
