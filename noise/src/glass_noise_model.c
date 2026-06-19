/*
 * Copyright 2023 Tyson B. Littenberg, Neil J. Cornish, Robert Rosati
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
#include <glass_math.h>

#include "glass_noise.h"


void map_array_to_noise_params(struct InstrumentModel *model)
{
    model->sacc12 = model->sacc[0];
    model->sacc21 = model->sacc[1];

    model->sacc23 = model->sacc[2];
    model->sacc32 = model->sacc[3];

    model->sacc13 = model->sacc[4];
    model->sacc31 = model->sacc[5];

    model->soms12 = model->soms[0];
    model->soms21 = model->soms[1];
        
    model->soms23 = model->soms[2];
    model->soms32 = model->soms[3];
    
    model->soms13 = model->soms[4];
    model->soms31 = model->soms[5];

}

void map_noise_params_to_array(struct InstrumentModel *model)
{
    model->sacc[0] = model->sacc12;
    model->sacc[1] = model->sacc21;
        
    model->sacc[2] = model->sacc23;
    model->sacc[3] = model->sacc32;
    
    model->sacc[4] = model->sacc13;
    model->sacc[5] = model->sacc31;

    model->soms[0] = model->soms12;
    model->soms[1] = model->soms21;
    
    model->soms[2] = model->soms23;
    model->soms[3] = model->soms32;

    model->soms[4] = model->soms13;
    model->soms[5] = model->soms31;
}

void map_array_to_foreground_params(struct ForegroundModel *model)
{
    model->Amp   = exp(model->sgal[0]);
    model->f1    = exp(model->sgal[1]);
    model->alpha = model->sgal[2];
    model->fk    = exp(model->sgal[3]);
    model->f2    = exp(model->sgal[4]);
}

void map_foreground_params_to_array(struct ForegroundModel *model)
{
    model->sgal[0] = log(model->Amp);
    model->sgal[1] = log(model->f1);
    model->sgal[2] = model->alpha;
    model->sgal[3] = log(model->fk);
    model->sgal[4] = log(model->f2);
}

void alloc_spline_model(struct SplineModel *model, int Ndata, int Nlayer, int Nchannel, int Nspline)
{
    model->Nchannel = Nchannel;
    
    model->psd = malloc(sizeof(struct Noise));
    model->spline=malloc(sizeof(struct Noise));
    
    alloc_noise(model->psd, Ndata, Nlayer, Nchannel);
    alloc_noise(model->spline, Nspline, Nlayer, Nchannel);
}

void alloc_instrument_model(struct InstrumentModel *model, int Ndata, int Nlayer, int Nchannel)
{
    model->Nlink=6;
    model->soms = malloc(model->Nlink*sizeof(double));
    model->sacc = malloc(model->Nlink*sizeof(double));
    model->psd = malloc(sizeof(struct Noise));
    alloc_noise(model->psd, Ndata, Nlayer, Nchannel);
    model->grid_cache = NULL;

    // Time-varying (drift) state. Disabled and zeroed by default so all
    // existing call sites behave exactly as before; the app enables it.
    model->drift_enabled = 0;
    model->soms_drift = calloc(model->Nlink, sizeof(double));
    model->sacc_drift = calloc(model->Nlink, sizeof(double));
    model->psd_slope = malloc(sizeof(struct Noise));
    alloc_noise(model->psd_slope, Ndata, Nlayer, Nchannel);
}

void alloc_foreground_model(struct ForegroundModel *model, int Ndata, int Nlayer, int Nchannel)
{
    model->Nparams=5;
    model->sgal = malloc(model->Nparams*sizeof(double));
    model->psd = malloc(sizeof(struct Noise));
    alloc_noise(model->psd, Ndata, Nlayer, Nchannel);
    model->grid_cache = NULL;
}

struct SGWBResponse* global_SGWBResponse = NULL;

void alloc_sgwb_model(struct SGWBModel *model, int Ndata, int Nlayer, int Nchannel, SGWB_t SGWB_type)
{
    model->SGWB_type = SGWB_type;
    model->Nparams = SGWB_TEMPLATE_NPARAMS[SGWB_type];
    model->params = malloc(model->Nparams*sizeof(double));
    model->psd = malloc(sizeof(struct Noise));
    alloc_noise(model->psd, Ndata, Nlayer, Nchannel);
    model->grid_cache = NULL;
    if (global_SGWBResponse == NULL) {
        model->R = malloc(sizeof(struct SGWBResponse)); // memory leak, only once per run
        alloc_pop_sgwb_response(model->R,"./sgwb_response_xyz2.dat");
        global_SGWBResponse = model->R;
    } else {
        model->R = global_SGWBResponse;
    }
}

void alloc_pop_sgwb_response(struct SGWBResponse* sgwbr, char* fname) {
    FILE* ff = fopen(fname,"r");
    if (ff == NULL) {
        printf("Could not open file %s ! SGWBResponse generation failed.",fname);
        exit(-2);
    }
    int *N = &sgwbr->N;
    if (fscanf(ff,"%d\n",N) != 1) {
        fclose(ff);
        printf("Failed to read first line of %s ! SGWBResponse generation failed.",fname);
        exit(-2);
    }
    if (*N > 3000 || *N < 10) {
        fclose(ff);
        printf("That seems wrong? I read N=%d . SGWBResponse generation failed.",*N);
        exit(-2);
    }

    // memory leak
    sgwbr->logf = malloc(sizeof(double) * *N); 
    sgwbr->logXX = malloc(sizeof(double) * *N);
    sgwbr->logYY = malloc(sizeof(double) * *N);
    sgwbr->logZZ = malloc(sizeof(double) * *N);
    sgwbr->asinhXY = malloc(sizeof(double) * *N);
    sgwbr->asinhXZ = malloc(sizeof(double) * *N);
    sgwbr->asinhYZ = malloc(sizeof(double) * *N);
    sgwbr->asinh_scale = 1e-7;
    for (int i=0;i<*N;i++) {
        double f,XX,YY,ZZ,XY,XZ,YZ;
        fscanf(ff,"%lf %lf %lf %lf %lf %lf %lf\n",&f,&XX,&YY,&ZZ,&XY,&XZ,&YZ);
        sgwbr->logf[i] = log(f);
        sgwbr->logXX[i] = log(XX);
        sgwbr->logYY[i] = log(YY);
        sgwbr->logZZ[i] = log(ZZ);
        sgwbr->asinhXY[i] = asinh(XY / sgwbr->asinh_scale);
        sgwbr->asinhXZ[i] = asinh(XZ / sgwbr->asinh_scale);
        sgwbr->asinhYZ[i] = asinh(YZ / sgwbr->asinh_scale);
    }
    fclose(ff);


    // memory leak
    sgwbr->spline_logRXX   = alloc_cubic_spline(*N);
    sgwbr->spline_logRYY   = alloc_cubic_spline(*N);
    sgwbr->spline_logRZZ   = alloc_cubic_spline(*N);
    sgwbr->spline_asinhRXY = alloc_cubic_spline(*N);
    sgwbr->spline_asinhRXZ = alloc_cubic_spline(*N);
    sgwbr->spline_asinhRYZ = alloc_cubic_spline(*N);
    initialize_cubic_spline(sgwbr->spline_logRXX, sgwbr->logf, sgwbr->logXX, SPLINE_EVEN_SAMPLED);
    initialize_cubic_spline(sgwbr->spline_logRYY, sgwbr->logf, sgwbr->logYY, SPLINE_EVEN_SAMPLED);
    initialize_cubic_spline(sgwbr->spline_logRZZ, sgwbr->logf, sgwbr->logZZ, SPLINE_EVEN_SAMPLED);
    initialize_cubic_spline(sgwbr->spline_asinhRXY, sgwbr->logf, sgwbr->asinhXY, SPLINE_EVEN_SAMPLED);
    initialize_cubic_spline(sgwbr->spline_asinhRXZ, sgwbr->logf, sgwbr->asinhXZ, SPLINE_EVEN_SAMPLED);
    initialize_cubic_spline(sgwbr->spline_asinhRYZ, sgwbr->logf, sgwbr->asinhYZ, SPLINE_EVEN_SAMPLED);

}

void free_sgwb_response(struct SGWBResponse *sgwbr) {
    free(sgwbr->logf);
    free(sgwbr->logXX);
    free(sgwbr->logYY);
    free(sgwbr->logZZ);
    free(sgwbr->asinhXY);
    free(sgwbr->asinhXZ);
    free(sgwbr->asinhYZ);
    free_cubic_spline(sgwbr->spline_logRXX);
    free_cubic_spline(sgwbr->spline_logRYY);
    free_cubic_spline(sgwbr->spline_logRZZ);
    free_cubic_spline(sgwbr->spline_asinhRXY);
    free_cubic_spline(sgwbr->spline_asinhRXZ);
    free_cubic_spline(sgwbr->spline_asinhRYZ);
    free(sgwbr);
}

void free_spline_model(struct SplineModel *model)
{
    free_noise(model->psd);
    free_noise(model->spline);
    free(model);
}

void free_instrument_model(struct InstrumentModel *model)
{
    if (model->grid_cache) free_instrument_model(model->grid_cache);
    free(model->soms);
    free(model->sacc);
    free(model->soms_drift);
    free(model->sacc_drift);
    free_noise(model->psd);
    free_noise(model->psd_slope);
    free(model);
}

void free_galaxy_modulation(struct GalaxyModulation *gm) {
    if(gm->XX_spline) free_cubic_spline(gm->XX_spline);
    if(gm->YY_spline) free_cubic_spline(gm->YY_spline);
    if(gm->ZZ_spline) free_cubic_spline(gm->ZZ_spline);
    if(gm->XY_spline) free_cubic_spline(gm->XY_spline);
    if(gm->XZ_spline) free_cubic_spline(gm->XZ_spline);
    if(gm->YZ_spline) free_cubic_spline(gm->YZ_spline);
    if(gm->cache_XX) free(gm->cache_XX);
    if(gm->cache_YY) free(gm->cache_YY);
    if(gm->cache_ZZ) free(gm->cache_ZZ);
    if(gm->cache_XY) free(gm->cache_XY);
    if(gm->cache_XZ) free(gm->cache_XZ);
    if(gm->cache_YZ) free(gm->cache_YZ);
}

void free_foreground_model(struct ForegroundModel *model)
{
    if (model->grid_cache) free_foreground_model(model->grid_cache);
    if (model->modulation) free_galaxy_modulation(model->modulation);
    free_noise(model->psd);
    free(model->sgal);
    free(model);
}

void free_sgwb_model(struct SGWBModel *model)
{
    if (model->grid_cache) free_sgwb_model(model->grid_cache);
    free_noise(model->psd);
    // TODO should actually only do this if we're the last SGWBModel instance.
    // for now, leak this memory
    // free_sgwb_response(model->R);
    free(model->params);
    free(model);
}

void copy_spline_model(struct SplineModel *origin, struct SplineModel *copy)
{
    copy->Nchannel = origin->Nchannel;

    //Spline parameters
    copy_noise(origin->spline,copy->spline);
    
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);
    
    //Model likelihood
    copy->logL = origin->logL;
    
    //Priors
    copy->Nmin = origin->Nmin;
    copy->Nmax = origin->Nmax;
}

void copy_instrument_model(struct InstrumentModel *origin, struct InstrumentModel *copy)
{
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);

    //Model likelihood
    copy->logL = origin->logL;

    //Instrument model Parameters
    copy->Nlink = origin->Nlink;
    memcpy(copy->soms, origin->soms, origin->Nlink*sizeof(double));
    memcpy(copy->sacc, origin->sacc, origin->Nlink*sizeof(double));

    //Time-varying (drift) state
    copy->drift_enabled = origin->drift_enabled;
    memcpy(copy->soms_drift, origin->soms_drift, origin->Nlink*sizeof(double));
    memcpy(copy->sacc_drift, origin->sacc_drift, origin->Nlink*sizeof(double));
    if (origin->drift_enabled)
        copy_noise(origin->psd_slope, copy->psd_slope);
}

void copy_foreground_model(struct ForegroundModel *origin, struct ForegroundModel *copy)
{
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);
    
    //Model likelihood
    copy->logL = origin->logL;

    //Foreground model Parameters
    copy->Tobs = origin->Tobs;
    copy->Nparams = origin->Nparams;
    memcpy(copy->sgal, origin->sgal, origin->Nparams*sizeof(double));
    
}
void copy_sgwb_model(struct SGWBModel *origin, struct SGWBModel *copy)
{
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);
    
    //Model likelihood
    copy->logL = origin->logL;

    //SGWB model Parameters
    copy->Tobs = origin->Tobs;
    copy->Nparams = origin->Nparams;
    copy->SGWB_type = origin->SGWB_type;
    memcpy(copy->params, origin->params, origin->Nparams*sizeof(double));
    
}

void update_spline_noise_model(struct SplineModel *model, int new_knot, int min_knot, int max_knot)
{
    struct Noise *psd = model->psd;
    struct Noise *spline = model->spline;
    
    /* find location in data vector of knot */
    double T = 1./(psd->f[1] - psd->f[0]);

    struct CubicSpline **cspline = malloc(model->Nchannel*sizeof(struct CubicSpline *));
    
    /* have to recompute the spline everywhere (derivatives on boundary) */
    for(int n=0; n<model->Nchannel; n++)
    {
        cspline[n] = alloc_cubic_spline(spline->N);
        initialize_cubic_spline(cspline[n], spline->f, spline->C[n][n], SPLINE_BINARY_SEARCH);
    }
    
    int imin = (int)((spline->f[min_knot]-psd->f[0])*T);
    int imax = (int)((spline->f[max_knot]-psd->f[0])*T);
    
    
    for(int i=imin; i<imax; i++)
    {
        for(int n=0; n<model->Nchannel; n++)
        {
            psd->C[n][n][i]=spline_interpolation(cspline[n],psd->f[i]);
            
            /*
             apply transfer function
             -this catches the sharp features in the spectrum from f/fstar
             -without needing to interpolate
             */
            psd->C[n][n][i]+=psd->transfer[i];
        }
    }
    invert_noise_covariance_matrix(psd);

    for(int n=0; n<model->Nchannel; n++)
    {
        free_cubic_spline(cspline[n]);
    }
    free(cspline);
}


void generate_spline_noise_model(struct SplineModel *model)
{
    struct Noise *psd = model->psd;
    struct Noise *spline = model->spline;
    
    for(int i=0; i<model->Nchannel; i++)
        CubicSplineGLASS(spline->N, spline->f, spline->C[i][i], psd->N, psd->f, psd->C[i][i]);

    //set a floor on Sn so likelihood doesn't go crazy
    for(int n=0; n<psd->N; n++)
    {
        /*
         apply transfer function
         -this catches the sharp features in the spectrum from f/fstar
         -without needing to interpolate
         */
        for(int i=0; i<model->Nchannel; i++)
            psd->C[i][i][n]+=psd->transfer[n];
        
    }
    invert_noise_covariance_matrix(psd);
}

void generate_instrument_noise_model(struct Orbit *orbit, struct InstrumentModel *model)
{
    double f,f2;
    double x;
    double cosx;
    double oms_transfer_function;
    double acc_transfer_function;
    double tdi_transfer_function;
    double Sacc;
    double Soms;
    double Sacc12,Sacc21,Sacc13,Sacc31,Sacc23,Sacc32;
    double Soms12,Soms21,Soms13,Soms31,Soms23,Soms32;

    map_array_to_noise_params(model);

    double acc_units = 2.81837551648e-19;//1./(PI2*CLIGHT)/(PI2*CLIGHT)
    double oms_units = 4.39256635604e-16;//PI2*PI2/CLIGHT/CLIGHT
    double acc_fmax  = 0.01; //Hz (maximum frequency for computing acc contribution
    
    for(int n=0; n<model->psd->N; n++)
    {
        f = model->psd->f[n];
        if(f == 0.0) {
            model->psd->C[0][0][n] = 0.0;
            model->psd->C[0][1][n] = 0.0;
            model->psd->C[0][2][n] = 0.0;
            model->psd->C[1][2][n] = 0.0;
            model->psd->C[2][2][n] = 0.0;
            continue;
        }
        x = f/orbit->fstar;
        f2= f*f;
        
        /*
        //at low frequency use linear approximation for trig functions
        if(x<0.1)
        {
            cosx = 1.0;
            tdi_transfer_function = x*x;
        }
        else
        {
            cosx = cos(x);
            tdi_transfer_function = noise_transfer_function(x);
        }
        */

        cosx = cos(x);
        tdi_transfer_function = noise_transfer_function(x);

        if(f<acc_fmax) //at f<~.1mHz acceleration noise is < 1% of the total noise budget
            acc_transfer_function = acc_units / f2 * (1.0 + pow(0.4e-3/f,2)) * (1.0 + pow(f/8.0e-3,4));
        else
            acc_transfer_function = 0.0;

        oms_transfer_function = oms_units * f2 * (1.0 + pow(2.0e-3/f,4));

        //absorb factor of -8 into tdi transfer function to cut down on multiplications
        tdi_transfer_function *= -8.0;
        
        switch(model->psd->Nchannel)
        {
            case 1:
                Sacc = model->sacc12 * acc_transfer_function;
                Soms = model->soms12 * oms_transfer_function;
                model->psd->C[0][0][n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Sacc, Soms);
                break;
            case 2:
                Sacc = model->sacc12 * acc_transfer_function;
                Soms = model->soms12 * oms_transfer_function;
                model->psd->C[0][0][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Sacc, Soms);
                model->psd->C[1][1][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Sacc, Soms);
                model->psd->C[0][1][n] = 0.0;
                break;
            case 3:
                
                if(f<acc_fmax) //at f<~.1mHz acceleration noise is < 1% of the total noise budget
                {
                    Sacc12 = model->sacc12 * acc_transfer_function;
                    Sacc21 = model->sacc21 * acc_transfer_function;
                    Sacc23 = model->sacc23 * acc_transfer_function;
                    Sacc32 = model->sacc32 * acc_transfer_function;
                    Sacc13 = model->sacc13 * acc_transfer_function;
                    Sacc31 = model->sacc31 * acc_transfer_function;
                }
                
                Soms12 = model->soms12 * oms_transfer_function;
                Soms21 = model->soms21 * oms_transfer_function;
                Soms23 = model->soms23 * oms_transfer_function;
                Soms32 = model->soms32 * oms_transfer_function;
                Soms13 = model->soms13 * oms_transfer_function;
                Soms31 = model->soms31 * oms_transfer_function;
            
                //Start from scratch
                model->psd->C[0][0][n] = 0.;
                model->psd->C[1][1][n] = 0.;
                model->psd->C[2][2][n] = 0.;

                model->psd->C[0][1][n] = 0.;
                model->psd->C[0][2][n] = 0.;
                model->psd->C[1][2][n] = 0.;

                //OMS noise
                model->psd->C[0][0][n] += (Soms12+Soms21+Soms13+Soms31)*.25;
                model->psd->C[1][1][n] += (Soms23+Soms32+Soms21+Soms12)*.25;
                model->psd->C[2][2][n] += (Soms31+Soms13+Soms32+Soms23)*.25;

                model->psd->C[0][1][n] += (Soms12 + Soms21)*.5;
                model->psd->C[0][2][n] += (Soms13 + Soms31)*.5;
                model->psd->C[1][2][n] += (Soms32 + Soms23)*.5;
                
                //Acceleration noise
                if(f<acc_fmax) //at f<~.1mHz acceleration noise is < 1% of the total noise budget
                {
                    //distribute through factors of 2
                    /*
                    model->psd->C[0][0][n] += 2. * ( (Sacc12 + Sacc13)*.5 + ((Sacc21 + Sacc31)*.5)*cosx*cosx );
                    model->psd->C[1][1][n] += 2. * ( (Sacc23 + Sacc21)*.5 + ((Sacc32 + Sacc12)*.5)*cosx*cosx );
                    model->psd->C[2][2][n] += 2. * ( (Sacc31 + Sacc32)*.5 + ((Sacc13 + Sacc23)*.5)*cosx*cosx );

                    model->psd->C[0][1][n] += 4. * ( (Sacc12 + Sacc21)*.5 );
                    model->psd->C[0][2][n] += 4. * ( (Sacc13 + Sacc31)*.5 );
                    model->psd->C[1][2][n] += 4. * ( (Sacc32 + Sacc23)*.5 );
                    */
                    model->psd->C[0][0][n] += ( (Sacc12 + Sacc13) + ((Sacc21 + Sacc31))*cosx*cosx );
                    model->psd->C[1][1][n] += ( (Sacc23 + Sacc21) + ((Sacc32 + Sacc12))*cosx*cosx );
                    model->psd->C[2][2][n] += ( (Sacc31 + Sacc32) + ((Sacc13 + Sacc23))*cosx*cosx );

                    model->psd->C[0][1][n] += 2.*(Sacc12 + Sacc21);
                    model->psd->C[0][2][n] += 2.*(Sacc13 + Sacc31);
                    model->psd->C[1][2][n] += 2.*(Sacc32 + Sacc23);
                }
                

                //TDI transfer functions (note tdi_transfer_function has abosrbed a factor of -8)
                model->psd->C[0][0][n] *= -2. * tdi_transfer_function;
                model->psd->C[1][1][n] *= -2. * tdi_transfer_function;
                model->psd->C[2][2][n] *= -2. * tdi_transfer_function;

                model->psd->C[0][1][n] *= tdi_transfer_function * cosx ;
                model->psd->C[0][2][n] *= tdi_transfer_function * cosx ;
                model->psd->C[1][2][n] *= tdi_transfer_function * cosx ;
                
                //Symmetry
                model->psd->C[1][0][n] = model->psd->C[0][1][n];
                model->psd->C[2][0][n] = model->psd->C[0][2][n];
                model->psd->C[2][1][n] = model->psd->C[1][2][n];
                
                break;
        }
        
        //Normalization
        for(int i=0; i<model->psd->Nchannel; i++)
            for(int j=0; j<model->psd->Nchannel; j++)
                model->psd->C[i][j][n] /= 4.0;
    }
}
//TODO: is it worth oversampling these very smooth functional forms...?

void generate_instrument_noise_model_wavelet(struct Wavelets *wdm, struct Orbit *orbit, struct InstrumentModel *model)
{
    int NF = wdm->NF;
    int NT = wdm->NT;
    int ND = NF * NT;

    // active layers (model->psd->C[i][j][k] holds layer imin+k)
    int imin = (int)round(model->psd->f[0]/wdm->df);
    int imax = (int)round(model->psd->f[model->psd->N-1]/wdm->df)+1;
    int Nactive = imax - imin;

    struct InstrumentModel *grid = model->grid_cache;
    if (grid && (grid->psd->N != Nactive ||
                 grid->psd->f[0] != (double)imin * wdm->df)) {
        // band or wdm changed: drop and re-alloc
        free_instrument_model(grid);
        grid = NULL;
    }
    if (!grid) {
        grid = malloc(sizeof(struct InstrumentModel));
        alloc_instrument_model(grid, Nactive, NF, 3);
        for(int n=0; n<Nactive; n++)
            grid->psd->f[n] = (double)(imin + n) * wdm->df;
        model->grid_cache = grid;
    }
    for(int i=0; i<grid->Nlink; i++)
    {
        grid->soms[i] = model->soms[i];
        grid->sacc[i] = model->sacc[i];
    }
    generate_instrument_noise_model(orbit, grid);

    // Layer variance = PSD(f_m) / ND, matching stationary_dft_psd_to_wdm_layer_var_approx
    double inv_ND = 1.0 / (double)ND;
    for(int n=0; n<3; n++)
        for(int m=n; m<3; m++)
            for(int k=0; k<Nactive; k++) {
                double v = grid->psd->C[n][m][k] * inv_ND;
                model->psd->C[n][m][k] = v;
                model->psd->C[m][n][k] = v;
            }

    // Time-varying instrument noise: linear drift per link coefficient.
    // soms_i(t) = soms_i*(1 + soms_drift_i*u), likewise sacc_i, with u in [-1,1].
    // The PSD is linear in each coefficient, so the covariance separates as
    //   C(t,f) = C(soms,sacc; f) + u * C(soms*drift, sacc*drift; f).
    // Reuse the same generator on the drift-scaled coefficients to get the slope.
    if (model->drift_enabled)
    {
        for(int i=0; i<grid->Nlink; i++)
        {
            grid->soms[i] = model->soms[i] * model->soms_drift[i];
            grid->sacc[i] = model->sacc[i] * model->sacc_drift[i];
        }
        generate_instrument_noise_model(orbit, grid);
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++)
                for(int k=0; k<Nactive; k++) {
                    double v = grid->psd->C[n][m][k] * inv_ND;
                    model->psd_slope->C[n][m][k] = v;
                    model->psd_slope->C[m][n][k] = v;
                }
    }
}

void generate_galactic_foreground_model(struct ForegroundModel *model)
{
    double f;
    double Sgal;
    
    map_array_to_foreground_params(model);
    
    for(int n=0; n<model->psd->N; n++)
    {
        f = model->psd->f[n];
        double fstar=CLIGHT/(PI2*LARM);
        double x = f/fstar;
        double t = 4.*x*x*sin(x)*sin(x); //transfer function

        //skip confusion noise calculation at high frequencies where contribution is negligible
        int high_f_check = 0;
        if(f>model->f1*10) high_f_check = 1;
        
        if(high_f_check) Sgal = 0.0;
        else Sgal = galaxy_foreground(f,model->Amp, model->f1, model->alpha, model->fk, model->f2);
        
        switch(model->psd->Nchannel)
        {
            case 1:
                model->psd->C[0][0][n] = Sgal;
                break;
            case 2:
                model->psd->C[0][0][n] = model->psd->C[1][1][n] = 1.5*t*Sgal;
                model->psd->C[0][1][n] = model->psd->C[1][0][n] = 0;
                break;
            case 3:
                model->psd->C[0][0][n] = model->psd->C[1][1][n] = model->psd->C[2][2][n] = t*Sgal;
                model->psd->C[0][1][n] = model->psd->C[0][2][n] = model->psd->C[1][2][n] = -0.5*t*Sgal;
                model->psd->C[1][0][n] = model->psd->C[2][0][n] = model->psd->C[2][1][n] = -0.5*t*Sgal;
                break;
        }
        
        //Normalization
        for(int i=0; i<model->psd->Nchannel; i++)
            for(int j=0; j<model->psd->Nchannel; j++)
                model->psd->C[i][j][n] /= 4.0;
    }
}

void generate_galactic_foreground_model_wavelet(struct Wavelets *wdm, struct ForegroundModel *model)
{
    /*
     * Sample the foreground PSD directly at the active-band layer-center
     * frequencies. The approx FFT->WDM map only reads layer-center bins,
     * so evaluating outside the active band was wasted work.
     */
    int NF = wdm->NF;
    int NT = wdm->NT;
    int ND = NF * NT;

    int imin = (int)round(model->psd->f[0]/wdm->df);
    int imax = (int)round(model->psd->f[model->psd->N-1]/wdm->df)+1;
    int Nactive = imax - imin;

    struct ForegroundModel *grid = model->grid_cache;
    if (grid && (grid->psd->N != Nactive ||
                 grid->psd->f[0] != (double)imin * wdm->df)) {
        free_foreground_model(grid);
        grid = NULL;
    }
    if (!grid) {
        grid = malloc(sizeof(struct ForegroundModel));
        alloc_foreground_model(grid, Nactive, NF, 3);
        for(int n=0; n<Nactive; n++)
            grid->psd->f[n] = (double)(imin + n) * wdm->df;
        model->grid_cache = grid;
    }

    // copy foreground parameters
    map_array_to_foreground_params(model);
    grid->Tobs  = model->Tobs;
    grid->Amp   = model->Amp;
    grid->f1    = model->f1;
    grid->alpha = model->alpha;
    grid->fk    = model->fk;
    grid->f2    = model->f2;
    map_foreground_params_to_array(grid);

    generate_galactic_foreground_model(grid);

    // galaxy_foreground diverges at f=0; only present in the grid if imin==0.
    if (imin == 0) {
        for(int n=0; n<3; n++)
            for(int m=0; m<3; m++)
                grid->psd->C[n][m][0] = 0.0;
    }

    // Undo the isotropic -1/2 baked into generate_galactic_foreground_model
    // off-diagonals (modulation(t) is applied later).
    for(int n=0; n<3; n++)
        for(int m=n+1; m<3; m++)
            for(int k=0; k<Nactive; k++)
                grid->psd->C[n][m][k] *= -2.0;

    // Layer variance = PSD(f_m) / ND, matching stationary_dft_psd_to_wdm_layer_var_approx
    double inv_ND = 1.0 / (double)ND;
    for(int n=0; n<3; n++)
        for(int m=n; m<3; m++)
            for(int k=0; k<Nactive; k++) {
                double v = grid->psd->C[n][m][k] * inv_ND;
                model->psd->C[n][m][k] = v;
                model->psd->C[m][n][k] = v;
            }
}


_Static_assert(SGWB_TEMPLATE_COUNT == 4, "Did you add an SGWB template? Implement its form in frequency in a new function here");
inline double sgwb_powerlaw(double f, const double* params) {
    static double fref = 25.0;
    double A     = pow(10.0,params[0]);
    double alpha = params[1];
    double prefactor = Hscale / (f*f*f);
    return prefactor*A*pow(f/fref,alpha);
}

// Pi & Sasaki, JCAP 2020 (arXiv:2005.12306), wide-Delta limit, eq. (3.29).
// Falls back to the asymptotic numerical value (3.33) for D >= 9 where the
// closed-form expression suffers from catastrophic cancellation.
// note that this is also used for SGWB_TEMPLATE_LOGNORMAL_FIXEDF
inline double sgwb_lognormal(double f, const double* params) {
    double A     = pow(10.0,params[0]);
    double fstar = pow(10.0,params[1]);
    double D     = pow(10.0,params[2]);
    double prefactor = Hscale / (f*f*f);
    double ft = f/fstar;
    if (D < 9.0) {
        double logft = log(ft);
        double logK  = logft + 1.5*D*D;
        double sqrtpi = sqrt(M_PI);
        double half_log32 = 0.5*log(1.5);
        double D2 = D*D;
        double t1 = 4.0/5.0/sqrtpi * ft*ft*ft * exp(9.0*D2/4.0) / D
                    * ((logK*logK + 0.5*D2) * glass_erfc((logK + half_log32)/D)
                       - D/sqrtpi * exp(-(logK + half_log32)*(logK + half_log32)/D2)
                         * (logK - half_log32));
        double t2 = 0.0659/D2 * ft*ft * exp(D2)
                    * exp(-(logft + D2 - 0.5*log(4.0/3.0))*(logft + D2 - 0.5*log(4.0/3.0))/D2);
        double t3 = (1.0/3.0) * sqrt(2.0)/sqrtpi * pow(ft,-4) * exp(8.0*D2) / D
                    * exp(-logft*logft/(2.0*D2))
                    * glass_erfc((4.0*D2 - logft + log(4.0))/(D*sqrt(2.0)));
        return prefactor * cgOr0 * A*A * (t1 + t2 + t3);
    } else {
        // Numerical asymptote from Pi & Sasaki eq. (3.33).
        return prefactor * cgOr0 * A*A * 0.783/1e3;
    }
}

inline double sgwb_phase_transition(double f, const double* params) {
    double rb = params[0];
    double b  = params[1];
    double Ap = pow(10.0,params[2]);
    double fp = pow(10.0,params[3]);
    double rb4 = rb*rb*rb*rb;
    double m = (9*rb4+b)/(rb4+1);
    double s = f/fp;
    double s4 = s*s*s*s;
    double s9 = s4*s4*s;
    double M = s9 * pow((1+rb4)/(rb4+s4),(9-b)/4) * pow((b+4)/(b+4-m + m*s*s),(b+4)/2);
    double prefactor = Hscale / (f*f*f);
    return Ap*M*prefactor;
}

void generate_sgwb_model(struct SGWBModel *model)
{
    double f;
    double Sgw;

    
    for(int n=0; n<model->psd->N; n++)
    {
        f = model->psd->f[n];
        
        _Static_assert(SGWB_TEMPLATE_COUNT == 4, "Did you add an SGWB template? Edit this switch case, it needs to be exhaustive.");
        switch(model->SGWB_type) {
            case SGWB_TEMPLATE_POWERLAW:
                Sgw = sgwb_powerlaw(f,model->params);
                break;
            case SGWB_TEMPLATE_LOGNORMAL:
                Sgw = sgwb_lognormal(f,model->params);
                break;
            case SGWB_TEMPLATE_PHASE_TRANSITION:
                Sgw = sgwb_phase_transition(f,model->params);
                break;
            case SGWB_TEMPLATE_LOGNORMAL_FIXEDF:
		double params[] = {model->params[0], log(1e-3)/log(10.0), 0.0};
                Sgw = sgwb_lognormal(f,params);
                break;
            default:
                fprintf(stderr,"Unimplemented SGWB type?\n\tTemplate index: %d\n\tTemplate name: %s\n",model->SGWB_type,SGWB_TEMPLATE_NAMES[model->SGWB_type]);
                exit(1);
                break;
        }
        
        switch(model->psd->Nchannel)
        {
            double Rxx, Ryy, Rzz, Rxy, Rxz, Ryz;
            double logf;
            case 1:
            case 2:
                fprintf(stderr, "Unimplemented error! SGWB covariance generation with Nchannels = %d\n",model->psd->Nchannel);
                exit(-3);
                break;
            case 3:
                if (f == 0.0) {
                    // TODO better zero check
                    model->psd->C[0][0][n] = 0.0;
                    model->psd->C[1][1][n] = 0.0;
                    model->psd->C[2][2][n] = 0.0;
                    model->psd->C[0][1][n] = model->psd->C[1][0][n] = 0.0;
                    model->psd->C[0][2][n] = model->psd->C[2][0][n] = 0.0;
                    model->psd->C[1][2][n] = model->psd->C[2][1][n] = 0.0;
                    continue;
                }
                logf = log(f);
                Rxx =  exp(spline_interpolation(model->R->spline_logRXX, logf));
                Ryy =  exp(spline_interpolation(model->R->spline_logRYY, logf));
                Rzz =  exp(spline_interpolation(model->R->spline_logRZZ, logf));
                Rxy = model->R->asinh_scale * sinh(spline_interpolation(model->R->spline_asinhRXY, logf));
                Rxz = model->R->asinh_scale * sinh(spline_interpolation(model->R->spline_asinhRXZ, logf));
                Ryz = model->R->asinh_scale * sinh(spline_interpolation(model->R->spline_asinhRYZ, logf));
                // note that, in equal-arm LISA, XYZ: Rxx = Ryy = Rzz, and Rxy = Rxz = Ryz = -0.5*Rxx
                model->psd->C[0][0][n] = Rxx*Sgw;
                model->psd->C[1][1][n] = Ryy*Sgw;
                model->psd->C[2][2][n] = Rzz*Sgw;
                model->psd->C[0][1][n] = model->psd->C[1][0][n] = Rxy*Sgw;
                model->psd->C[0][2][n] = model->psd->C[2][0][n] = Rxz*Sgw;
                model->psd->C[1][2][n] = model->psd->C[2][1][n] = Ryz*Sgw;
                break;
        }
    }
}
void generate_sgwb_model_wavelet(struct Wavelets* wdm, struct SGWBModel *model)
{
    // Sample directly at active-band layer-center frequencies (see comment in
    // generate_instrument_noise_model_wavelet).
    int NF = wdm->NF;
    int NT = wdm->NT;
    int ND = NF * NT;

    int imin = (int)round(model->psd->f[0]/wdm->df);
    int imax = (int)round(model->psd->f[model->psd->N-1]/wdm->df)+1;
    int Nactive = imax - imin;

    struct SGWBModel *grid = model->grid_cache;
    if (grid && (grid->psd->N != Nactive ||
                 grid->psd->f[0] != (double)imin * wdm->df)) {
        free_sgwb_model(grid);
        grid = NULL;
    }
    if (!grid) {
        grid = malloc(sizeof(struct SGWBModel));
        alloc_sgwb_model(grid, Nactive, NF, 3, model->SGWB_type);
        for(int n=0; n<Nactive; n++)
            grid->psd->f[n] = (double)(imin + n) * wdm->df;
        model->grid_cache = grid;
    }

    grid->Nparams = model->Nparams;
    grid->Tobs = model->Tobs;
    for (int i=0; i<model->Nparams; i++)
        grid->params[i] = model->params[i];

    generate_sgwb_model(grid);

    // DC layer is only in the grid if the active band starts at 0.
    if (imin == 0) {
        for(int n=0; n<3; n++)
            for(int m=0; m<3; m++)
                grid->psd->C[n][m][0] = 0.0;
    }

    // Layer variance = PSD(f_m) / ND, matching stationary_dft_psd_to_wdm_layer_var_approx
    double inv_ND = 1.0 / (double)ND;
    for(int n=0; n<3; n++)
        for(int m=n; m<3; m++)
            for(int k=0; k<Nactive; k++) {
                double v = grid->psd->C[n][m][k] * inv_ND;
                model->psd->C[n][m][k] = v;
                model->psd->C[m][n][k] = v;
            }
}

void generate_full_covariance_matrix(struct Noise *full, struct Noise *component, int Nchannel)
{
    for(int n=0; n<full->N; n++)
    {
        for(int i=0; i<Nchannel; i++)
        {
            for(int j=0; j<Nchannel; j++)
            {
                full->C[i][j][n] += component->C[i][j][n];
            }
        }
    }
}


// Q=1 is the degenerate (full-resolution) case of the coarse path: Ncoarse=NT,
// k = q + jrel*NT, which equals wavelet_pixel_to_index(i, lmin+jrel) - kmin.
// The full-res LL (my_noise_log_likelihood_wavelet) reads noise->C[k] via that
// same linearization, so this delegation is layout-equivalent.
void generate_full_dynamic_covariance_matrix(struct Wavelets *wdm, struct InstrumentModel *inst, struct ForegroundModel *conf, struct SGWBModel *sgwb, struct Noise *full)
{
    generate_full_dynamic_covariance_matrix_coarse(wdm, 1, inst, conf, sgwb, full);
}

void generate_full_stationary_covariance_matrix(struct Wavelets *wdm, struct InstrumentModel *inst, struct ForegroundModel *conf, struct SGWBModel *sgwb, struct Noise *full)
{
    generate_full_stationary_covariance_matrix_coarse(wdm, 1, inst, conf, sgwb, full);
}

void generate_full_dynamic_covariance_matrix_coarse(struct Wavelets *wdm, int Q, struct InstrumentModel *inst, struct ForegroundModel *conf, struct SGWBModel *sgwb, struct Noise *coarse)
{
    int jmin=(int)round(inst->psd->f[0]/wdm->df);
    int jmax=(int)round(inst->psd->f[inst->psd->N-1]/wdm->df)+1;
    int Ncoarse = wdm->NT / Q;

    if (conf) galaxy_modulation_cache_for_Q(conf->modulation, wdm, Q);

    // Time-varying instrument noise: C_inst(t,f) = psd(f) + u(t)*psd_slope(f),
    // with u in [-1,1] mapped to the cell-center fraction of the run. u is
    // anchored to absolute time, so it is independent of the coarse factor Q.
    int drift = inst->drift_enabled;

    for(int j=jmin; j<jmax; j++)
    {
        int jrel = j - jmin;
        for(int q=0; q<Ncoarse; q++)
        {
            double mxx=0, myy=0, mzz=0, mxy=0, mxz=0, myz=0;
            if (conf) {
                mxx = conf->modulation->cache_XX[q];
                myy = conf->modulation->cache_YY[q];
                mzz = conf->modulation->cache_ZZ[q];
                mxy = conf->modulation->cache_XY[q];
                mxz = conf->modulation->cache_XZ[q];
                myz = conf->modulation->cache_YZ[q];
            }

            int k = q + jrel*Ncoarse;

            //instrument noise (linear drift in time when enabled)
            coarse->C[0][0][k] = inst->psd->C[0][0][jrel];
            coarse->C[1][1][k] = inst->psd->C[1][1][jrel];
            coarse->C[2][2][k] = inst->psd->C[2][2][jrel];
            coarse->C[0][1][k] = inst->psd->C[0][1][jrel];
            coarse->C[0][2][k] = inst->psd->C[0][2][jrel];
            coarse->C[1][2][k] = inst->psd->C[1][2][jrel];

            if (drift) {
                double u = (2.0*q + 1.0)/(double)Ncoarse - 1.0;
                coarse->C[0][0][k] += u*inst->psd_slope->C[0][0][jrel];
                coarse->C[1][1][k] += u*inst->psd_slope->C[1][1][jrel];
                coarse->C[2][2][k] += u*inst->psd_slope->C[2][2][jrel];
                coarse->C[0][1][k] += u*inst->psd_slope->C[0][1][jrel];
                coarse->C[0][2][k] += u*inst->psd_slope->C[0][2][jrel];
                coarse->C[1][2][k] += u*inst->psd_slope->C[1][2][jrel];
            }

            //modulated galactic foreground (sampled at coarse-cell midpoint)
            if (conf) {
                coarse->C[0][0][k] += conf->psd->C[0][0][jrel]*mxx;
                coarse->C[1][1][k] += conf->psd->C[1][1][jrel]*myy;
                coarse->C[2][2][k] += conf->psd->C[2][2][jrel]*mzz;
                coarse->C[0][1][k] += conf->psd->C[0][1][jrel]*mxy;
                coarse->C[0][2][k] += conf->psd->C[0][2][jrel]*mxz;
                coarse->C[1][2][k] += conf->psd->C[1][2][jrel]*myz;
            }

            //stationary stochastic background
            if (sgwb) {
                coarse->C[0][0][k] += sgwb->psd->C[0][0][jrel];
                coarse->C[1][1][k] += sgwb->psd->C[1][1][jrel];
                coarse->C[2][2][k] += sgwb->psd->C[2][2][jrel];
                coarse->C[0][1][k] += sgwb->psd->C[0][1][jrel];
                coarse->C[0][2][k] += sgwb->psd->C[0][2][jrel];
                coarse->C[1][2][k] += sgwb->psd->C[1][2][jrel];
            }

            coarse->C[1][0][k] = coarse->C[0][1][k];
            coarse->C[2][0][k] = coarse->C[0][2][k];
            coarse->C[2][1][k] = coarse->C[1][2][k];
        }
    }
}

void generate_full_stationary_covariance_matrix_coarse(struct Wavelets *wdm, int Q, struct InstrumentModel *inst, struct ForegroundModel *conf, struct SGWBModel *sgwb, struct Noise *coarse)
{
    int jmin=(int)round(inst->psd->f[0]/wdm->df);
    int jmax=(int)round(inst->psd->f[inst->psd->N-1]/wdm->df)+1;
    int Ncoarse = wdm->NT / Q;

    // j outer, q inner for contiguous writes. See dynamic variant above.
    for(int j=jmin; j<jmax; j++)
    {
        int jrel = j - jmin;
        for(int q=0; q<Ncoarse; q++)
        {
            int k = q + jrel*Ncoarse;

            coarse->C[0][0][k] = inst->psd->C[0][0][jrel];
            coarse->C[1][1][k] = inst->psd->C[1][1][jrel];
            coarse->C[2][2][k] = inst->psd->C[2][2][jrel];
            coarse->C[0][1][k] = inst->psd->C[0][1][jrel];
            coarse->C[0][2][k] = inst->psd->C[0][2][jrel];
            coarse->C[1][2][k] = inst->psd->C[1][2][jrel];

            // mirror the stationary-path foreground convention from
            // generate_full_stationary_covariance_matrix
            if (conf) {
                coarse->C[0][0][k] += conf->psd->C[0][0][jrel];
                coarse->C[1][1][k] += conf->psd->C[1][1][jrel];
                coarse->C[2][2][k] += conf->psd->C[2][2][jrel];
                coarse->C[0][1][k] -= conf->psd->C[0][1][jrel]/2.;
                coarse->C[0][2][k] -= conf->psd->C[0][2][jrel]/2.;
                coarse->C[1][2][k] -= conf->psd->C[1][2][jrel]/2.;
            }

            if (sgwb) {
                coarse->C[0][0][k] += sgwb->psd->C[0][0][jrel];
                coarse->C[1][1][k] += sgwb->psd->C[1][1][jrel];
                coarse->C[2][2][k] += sgwb->psd->C[2][2][jrel];
                coarse->C[0][1][k] += sgwb->psd->C[0][1][jrel];
                coarse->C[0][2][k] += sgwb->psd->C[0][2][jrel];
                coarse->C[1][2][k] += sgwb->psd->C[1][2][jrel];
            }

            coarse->C[1][0][k] = coarse->C[0][1][k];
            coarse->C[2][0][k] = coarse->C[0][2][k];
            coarse->C[2][1][k] = coarse->C[1][2][k];
        }
    }
}

void alloc_coarse_stats(struct CoarseStats *s, int Nlayer, int Q, int NT)
{
    int Ncoarse = NT / Q;
    size_t Pdim = (size_t)Nlayer * (size_t)Ncoarse;
    // Allocate via non-const handles, expose as const through the struct.
    double *Pxx = malloc(Pdim*sizeof(double));
    double *Pyy = malloc(Pdim*sizeof(double));
    double *Pzz = malloc(Pdim*sizeof(double));
    double *Pxy = malloc(Pdim*sizeof(double));
    double *Pxz = malloc(Pdim*sizeof(double));
    double *Pyz = malloc(Pdim*sizeof(double));
    // Default Qeff to the integer Q in every cell, so the likelihood matches
    // the pre-ws_approx behavior unless precompute_coarse_Qeff overwrites it.
    double *Qeff = malloc(Pdim*sizeof(double));
    for(size_t i=0; i<Pdim; i++) Qeff[i] = (double)Q;
    s->Pxx = Pxx; s->Pyy = Pyy; s->Pzz = Pzz;
    s->Pxy = Pxy; s->Pxz = Pxz; s->Pyz = Pyz;
    s->Qeff = Qeff;
    s->Q = Q;
    s->Ncoarse = Ncoarse;
    s->Nlayer = Nlayer;
    // Default: sum over every coarse cell. Narrow [qlo,qhi) to crop edge cells.
    s->qlo = 0;
    s->qhi = Ncoarse;
}

void free_coarse_stats(struct CoarseStats *s)
{
    free((void*)s->Pxx); free((void*)s->Pyy); free((void*)s->Pzz);
    free((void*)s->Pxy); free((void*)s->Pxz); free((void*)s->Pyz);
    free((void*)s->Qeff);
    s->Pxx = s->Pyy = s->Pzz = NULL;
    s->Pxy = s->Pxz = s->Pyz = NULL;
    s->Qeff = NULL;
}

void coarse_grain_wavelet_data(struct Data *data, struct CoarseStats *stats)
{
    struct Wavelets *wdm = data->wdm;
    struct TDI *tdi = data->tdi;
    int Q = stats->Q;
    int Ncoarse = stats->Ncoarse;
    double invQ = 1.0/(double)Q;
    // Populator owns mutable views of the arrays it allocated.
    double *Pxx = (double*)stats->Pxx;
    double *Pyy = (double*)stats->Pyy;
    double *Pzz = (double*)stats->Pzz;
    double *Pxy = (double*)stats->Pxy;
    double *Pxz = (double*)stats->Pxz;
    double *Pyz = (double*)stats->Pyz;

    for(int j=data->lmin; j<data->lmax; j++)
    {
        int jrel = j - data->lmin;
        for(int q=0; q<Ncoarse; q++)
        {
            double sxx=0, syy=0, szz=0, sxy=0, sxz=0, syz=0;
            for(int i=q*Q; i<(q+1)*Q; i++)
            {
                int k_full;
                wavelet_pixel_to_index(wdm, i, j, &k_full);
                k_full -= wdm->kmin;
                double wx = tdi->X[k_full];
                double wy = tdi->Y[k_full];
                double wz = tdi->Z[k_full];
                sxx += wx*wx;
                syy += wy*wy;
                szz += wz*wz;
                sxy += wx*wy;
                sxz += wx*wz;
                syz += wy*wz;
            }
            int k_coarse = q + jrel*Ncoarse;
            Pxx[k_coarse] = sxx*invQ;
            Pyy[k_coarse] = syy*invQ;
            Pzz[k_coarse] = szz*invQ;
            Pxy[k_coarse] = sxy*invQ;
            Pxz[k_coarse] = sxz*invQ;
            Pyz[k_coarse] = syz*invQ;
        }
    }
}

void precompute_coarse_Qeff(struct Data *data, struct CoarseStats *stats)
{
    struct Wavelets *wdm = data->wdm;
    struct Noise *noise = data->noise;   // injected, full-resolution (Q=1) covariance
    int Q = stats->Q;
    int Ncoarse = stats->Ncoarse;
    // Populator owns a mutable view of the array it allocated.
    double *Qeff = (double*)stats->Qeff;

    if (Q <= 1) return; // each cell is one fine pixel: Qeff == Q == 1 already.

    for(int j=data->lmin; j<data->lmax; j++)
    {
        int jrel = j - data->lmin;
        for(int q=0; q<Ncoarse; q++)
        {
            // First two moments of the injected per-fine-pixel variance within
            // the coarse cell, per diagonal channel. With C_a(t_n) the injected
            // variance at fine pixel n, the heterogeneous-variance sum
            // sum_n C_a(t_n) z_n^2 (z_n ~ N(0,1)) is moment-matched to a scaled
            // chi^2 with Welch-Satterthwaite dof (sum_n C_a)^2 / sum_n C_a^2.
            double s1x=0, s2x=0, s1y=0, s2y=0, s1z=0, s2z=0;
            for(int i=q*Q; i<(q+1)*Q; i++)
            {
                int k_full;
                wavelet_pixel_to_index(wdm, i, j, &k_full);
                k_full -= wdm->kmin;
                double cxx = noise->C[0][0][k_full];
                double cyy = noise->C[1][1][k_full];
                double czz = noise->C[2][2][k_full];
                s1x += cxx; s2x += cxx*cxx;
                s1y += cyy; s2y += cyy*cyy;
                s1z += czz; s2z += czz*czz;
            }
            double qxx = s1x*s1x/s2x;
            double qyy = s1y*s1y/s2y;
            double qzz = s1z*s1z/s2z;
            // Single scalar dof per cell (the likelihood applies one Qd to all
            // channels and cross terms), taken as the channel average.
            Qeff[q + jrel*Ncoarse] = (qxx + qyy + qzz)/3.0;
        }
    }
}

double noise_log_likelihood(struct Data *data, struct Noise *noise)
{
    double logL = 0.0;

    struct TDI *tdi = data->tdi;

    int N = data->NFFT;
    
    switch(data->Nchannel)
    {
        case 1:
            logL += -0.5*fourier_nwip(tdi->X, tdi->X, noise->invC[0][0], N);
            break;
        case 2:
            logL += -0.5*fourier_nwip(tdi->A, tdi->A, noise->invC[0][0], N);
            logL += -0.5*fourier_nwip(tdi->E, tdi->E, noise->invC[1][1], N);
            break;
        case 3:
            logL += -0.5*fourier_nwip(tdi->X, tdi->X, noise->invC[0][0], N);
            logL += -0.5*fourier_nwip(tdi->Y, tdi->Y, noise->invC[1][1], N);
            logL += -0.5*fourier_nwip(tdi->Z, tdi->Z, noise->invC[2][2], N);
            logL += -fourier_nwip(tdi->X, tdi->Y, noise->invC[0][1], N);
            logL += -fourier_nwip(tdi->X, tdi->Z, noise->invC[0][2], N);
            logL += -fourier_nwip(tdi->Y, tdi->Z, noise->invC[1][2], N);
            break;
    }
    for(int n=0; n<N; n++)
        logL -= noise->logdetC[n];
    
    return logL;
}

static inline double wavelet_nwip_linear(const double* __restrict a, const double* __restrict b, const double* __restrict invC, size_t N)
{
    // TODO: use lapack/blas here? this is just a.b.invC
    double arg = 0.0;
    #pragma omp simd reduction(+:arg)
    for(size_t k=0; k<N; k++)
    {
        arg += a[k]*b[k]*invC[k];
    }
    return arg;
}

// Per-pixel multivariate Gaussian log-likelihood contribution. Q absorbs the
// Welch sufficient-statistic weighting: Q=1 is the standard Gaussian (one fine
// pixel = one degree of freedom); Q>1 is the Wishart sufficient-statistic form
// with the in-cell sample (cross-)covariances P_ab as observed data.
//
// Note about this likelihood:
// When Qd is 1, this is just the normal Gaussian / Whittle likelihood
// When Qd  > 1, this is mostly equivalent to a Gamma likelihood
// (see appendix of WDM note!)
// But essentially, this only differs from the Gamma by dropping terms that only depend
// on Q or the data, but are noise-model-independent. These don't affect the posterior
//
// TODO: in practice, the non-stationarity changes this distribution from a Gamma
// see WDM note. It's an infinite sum of Gammas
//
// The dropped pieces are the normalization of the Wishart/Gamma sufficient-
// statistic density. Treating the cell statistic as A = Qd*P ~ Wishart_p(Cbar, Qd)
// with p=3 channels, the full log-density is
//   ((Qd-p-1)/2) log|A| - 1/2 tr(Cbar^-1 A) - (Qd p / 2) log2 - logGamma_p(Qd/2),
// and only the two terms kept below (-1/2 tr(Cbar^-1 A) and the -Qd/2 log|Cbar|
// hidden in logdetC) depend on the *model* when Qd is a fixed constant. The
// commented block restores the remaining Qd-dependent normalization; see the
// NOTE there for when it is actually needed.

static inline double per_pixel_logL_contribution(
        double Qd,
        double pxx, double pyy, double pzz,
        double pxy, double pxz, double pyz,
        double inv_xx, double inv_yy, double inv_zz,
        double inv_xy, double inv_xz, double inv_yz,
        double logdetC)
{
    double s = 0.0;
    s -= 0.5*Qd*inv_xx*pxx;
    s -= 0.5*Qd*inv_yy*pyy;
    s -= 0.5*Qd*inv_zz*pzz;
    s -=     Qd*inv_xy*pxy;
    s -=     Qd*inv_xz*pxz;
    s -=     Qd*inv_yz*pyz;
    s -= 0.5*Qd*logdetC;

    /* ---------------------------------------------------------------------
     * Qd-dependent normalization of the Wishart (p=3) sufficient-statistic
     * density. These terms are CONSTANT across the chain when Qd is fixed --
     * which it is both for plain coarse-graining (Qd = integer Q) and for the
     * --ws-approx path (Qeff frozen at the injected covariance) -- so they
     * cancel in every MCMC acceptance ratio and are omitted by default.
     *
     * NOTE: restore this block if Qd ever becomes model-dependent, i.e. if Qeff
     * is recomputed each step from the *current* noise parameters rather than
     * frozen at injection. It is also required whenever an absolute, correctly
     * normalized logL is needed (e.g. thermodynamic-integration evidence, or
     * comparing models/segments that carry different effective dof). Valid only
     * for Qd > p-1 = 2 (multivariate gamma domain); the channel-averaged Qeff
     * can dip below 2 in strongly confusion-dominated cells, so guard before use.
     *
     * const double p = 3.0;
     * double detP = pxx*(pyy*pzz - pyz*pyz)
     *             - pxy*(pxy*pzz - pyz*pxz)
     *             + pxz*(pxy*pyz - pyy*pxz);
     * double x = 0.5*Qd;
     * double logGamma_p = 1.5*log(M_PI)             // pi^{p(p-1)/4}, p=3 -> pi^{3/2}
     *                   + lgamma(x)
     *                   + lgamma(x - 0.5)
     *                   + lgamma(x - 1.0);
     * double logdetA = p*log(Qd) + log(detP);       // |A| = Qd^p |P|
     * s += 0.5*(Qd - p - 1.0)*logdetA;
     * s -= 0.5*Qd*p*log(2.0);
     * s -= logGamma_p;
     * --------------------------------------------------------------------- */

    return s;
}

double my_noise_log_likelihood_wavelet(struct Data *data, struct Noise *noise) {
    struct TDI *tdi = data->tdi;
    const double log2pi = log(2*M_PI);
    int NT = data->wdm->NT;
    int Nlayer = data->Nlayer;
    double logL = 0.0;
    // j outer / i inner so k = i + j*NT increments by 1 in the inner loop
    // (matches the column-major wavelet layout, equivalent to
    // wavelet_pixel_to_index(i, lmin+j) - kmin).
    #pragma omp simd reduction(+:logL)
    for (int j=0; j<Nlayer; j++) {
        for (int i=0; i<NT; i++) {
            int k = i + j*NT;
            double wx = tdi->X[k];
            double wy = tdi->Y[k];
            double wz = tdi->Z[k];
            logL += per_pixel_logL_contribution(1.0,
                    wx*wx, wy*wy, wz*wz,
                    wx*wy, wx*wz, wy*wz,
                    noise->invC[0][0][k], noise->invC[1][1][k], noise->invC[2][2][k],
                    noise->invC[0][1][k], noise->invC[0][2][k], noise->invC[1][2][k],
                    noise->logdetC[k]);
        }
    }
    logL -= 0.5 * 3 * data->N * log2pi;

    return logL;
}

double my_noise_log_likelihood_wavelet_coarse(struct Data *data, struct Noise *coarse_noise, const struct CoarseStats *stats)
{
    const double log2pi = log(2*M_PI);
    int Nlayer = data->Nlayer;
    int Ncoarse = stats->Ncoarse;
    // Per-cell effective dof. Without --ws-approx this array is uniformly the
    // integer Q (filled in alloc_coarse_stats), so the result is identical to
    // the old scalar-Q path; with --ws-approx it holds the frozen Qeff <= Q.
    const double *Qeff = stats->Qeff;
    const double *Pxx = stats->Pxx;
    const double *Pyy = stats->Pyy;
    const double *Pzz = stats->Pzz;
    const double *Pxy = stats->Pxy;
    const double *Pxz = stats->Pxz;
    const double *Pyz = stats->Pyz;
    // Time-pixel crop: sum only coarse cells in [qlo, qhi). Defaults to the full
    // [0, Ncoarse) set in alloc_coarse_stats.
    int qlo = stats->qlo;
    int qhi = stats->qhi;
    double logL = 0.0;
    #pragma omp simd reduction(+:logL)
    for (int j=0; j<Nlayer; j++) {
        for (int q=qlo; q<qhi; q++) {
            int k = q + j*Ncoarse;
            logL += per_pixel_logL_contribution(Qeff[k],
                    Pxx[k], Pyy[k], Pzz[k],
                    Pxy[k], Pxz[k], Pyz[k],
                    coarse_noise->invC[0][0][k], coarse_noise->invC[1][1][k], coarse_noise->invC[2][2][k],
                    coarse_noise->invC[0][1][k], coarse_noise->invC[0][2][k], coarse_noise->invC[1][2][k],
                    coarse_noise->logdetC[k]);
        }
    }
    // constant uses the fine-pixel count actually included: each of the (qhi-qlo)
    // summed cells per layer carries Q fine Gaussian dof, over 3 channels. With
    // the full [0,Ncoarse) range this equals 3 * Nlayer * Ncoarse * Q == 3*data->N.
    logL -= 0.5 * 3.0 * (double)Nlayer * (double)(qhi - qlo) * (double)stats->Q * log2pi;

    return logL;
}

double my_noise_log_likelihood(struct Data *data, struct Noise *noise)
{
    double logL = 0.0;
    const double log2pi = log(2*M_PI);
    
    struct TDI *tdi = data->tdi;
    
    int N = data->NFFT;
    #pragma omp simd reduction(+:logL)
    for (int n=0; n<N; n++) {
        double rex = tdi->X[2*n], imx = tdi->X[2*n+1];
        double rey = tdi->Y[2*n], imy = tdi->Y[2*n+1];
        double rez = tdi->Z[2*n], imz = tdi->Z[2*n+1];
        double xx = rex*rex + imx*imx;
        double yy = rey*rey + imy*imy;
        double zz = rez*rez + imz*imz;
        double xy = rex*rey + imx*imy; // note this is really (x*conj(y) + y*conj(x))/2
        double xz = rex*rez + imx*imz;
        double yz = rey*rez + imy*imz;
        logL -= 0.5*xx*noise->invC[0][0][n]*2;
        logL -= 0.5*yy*noise->invC[1][1][n]*2;
        logL -= 0.5*zz*noise->invC[2][2][n]*2;
        logL -=     xy*noise->invC[0][1][n]*2;
        logL -=     xz*noise->invC[0][2][n]*2;
        logL -=     yz*noise->invC[1][2][n]*2;
        logL -= noise->logdetC[n];
    }
    logL -= 3 * N * log2pi;
    
    return logL;
}

double noise_log_likelihood_wavelet(struct Data *data, struct Noise *noise)
{
    double logL = 0.0;
    
    struct TDI *tdi = data->tdi;
    
    int N = data->N;
    
    switch(data->Nchannel)
    {
        case 1:
            logL += -0.5*wavelet_nwip_linear(tdi->X, tdi->X, noise->invC[0][0], N);
            break;
        case 2:
            logL += -0.5*wavelet_nwip_linear(tdi->A, tdi->A, noise->invC[0][0], N);
            logL += -0.5*wavelet_nwip_linear(tdi->E, tdi->E, noise->invC[1][1], N);
            break;
        case 3:
            logL += -0.5*wavelet_nwip_linear(tdi->X, tdi->X, noise->invC[0][0], N);
            logL += -0.5*wavelet_nwip_linear(tdi->Y, tdi->Y, noise->invC[1][1], N);
            logL += -0.5*wavelet_nwip_linear(tdi->Z, tdi->Z, noise->invC[2][2], N);
            logL += -wavelet_nwip_linear(tdi->X, tdi->Y, noise->invC[0][1], N);
            logL += -wavelet_nwip_linear(tdi->X, tdi->Z, noise->invC[0][2], N);
            logL += -wavelet_nwip_linear(tdi->Y, tdi->Z, noise->invC[1][2], N);
            break;
    }
    for(int n=0; n<N; n++)
        logL -= 0.5*noise->logdetC[n];

    return logL;
    // this comment constitutes an offering to the deity responsible for the correct factors of 2
}

double noise_delta_log_likelihood(struct Data *data, struct SplineModel *model_x, struct SplineModel *model_y, double fmin, double fmax,int ic)
{
    double dlogL = 0.0;
    
    struct TDI *tdi = data->tdi;
    struct Noise *psd_x = model_x->psd;
    struct Noise *psd_y = model_y->psd;

    int N = (int)floor((fmax-fmin)*data->T);
    int imin = (int)floor((fmin-psd_x->f[0])*data->T);
    if(imin<0)imin=0;
    
    /* remove contribution for current state x */
    dlogL -= -0.5*fourier_nwip(tdi->A+2*imin, tdi->A+2*imin, psd_x->invC[0][0]+imin, N);
    dlogL -= -0.5*fourier_nwip(tdi->E+2*imin, tdi->E+2*imin, psd_x->invC[1][1]+imin, N);
    for(int n=imin; n<imin+N; n++)
        dlogL += psd_x->logdetC[n];

    /* add contribution for proposed state y */
    dlogL += -0.5*fourier_nwip(tdi->A+2*imin, tdi->A+2*imin, psd_y->invC[0][0]+imin, N);
    dlogL += -0.5*fourier_nwip(tdi->E+2*imin, tdi->E+2*imin, psd_y->invC[1][1]+imin, N);
    for(int n=imin; n<imin+N; n++)
        dlogL -= psd_y->logdetC[n];
    
    return dlogL;
}

void initialize_spline_model(struct Orbit *orbit, struct Data *data, struct SplineModel *model, int Nspline)
{
    
    // Initialize data models
    alloc_spline_model(model, data->NFFT, data->Nlayer, data->Nchannel, Nspline);
    
    //set max and min spline points
    model->Nmin = MIN_SPLINE_STENCIL;
    model->Nmax = Nspline;
    model->Nchannel = data->Nchannel;
    
    //set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
    {
        double f = data->fmin + (double)n/data->T;
        double Spm, Sop;
        get_noise_levels("sangria", f, &Spm, &Sop);
        model->psd->f[n] = f;
        if(model->Nchannel==2) model->psd->transfer[n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);//noise_transfer_function(f/orbit->fstar);
        else model->psd->transfer[n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
    }
    
    //divide into Nspline control points
    double logdf = (log(data->fmax) - log(data->fmin))/(Nspline-1);
    for(int i=0; i<Nspline; i++)
    {
        double f = exp(log(data->fmin) + (double)i*logdf);
        model->spline->f[i] = f;
        
        for(int n=0; n<model->Nchannel; n++)
        {
            /* initialize model to theoretical level without transfer function applied */
            for(int m=0; m<model->Nchannel; m++)
                model->spline->C[n][m][i] = 0.0;
        }
    }
    //shift first spline control point by half a bin to avoid rounding problems
    model->spline->f[0] -= 0.5/data->T;
    
    generate_spline_noise_model(model);
    model->logL = noise_log_likelihood(data, model->psd);
}

void initialize_instrument_model(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model)
{
    // initialize data models
    alloc_instrument_model(model, data->NFFT, data->Nlayer, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = data->fmin + (double)n/data->T;

    // initialize noise levels
    /* NOTE: changed to match wavelet */
    /*
    for(int i=0; i<model->Nlink; i++)
    {
        model->soms[i] = 1.28e-22;
        model->sacc[i] = 5.76e-30;
    }
    */
    for(int i=0; i<model->Nlink; i++)
    {
        model->soms[i] = 1.28e-22;
        model->sacc[i] = 5.76e-30;
    }
    
    // get noise covariance matrix for initial parameters
    generate_instrument_noise_model(orbit,model);
}

void initialize_instrument_model_wavelet(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model)
{
    // wavelet basis
    struct Wavelets *wdm = data->wdm;

    // initialize data models
    alloc_instrument_model(model, data->lmax - data->lmin, data->Nlayer, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = (data->lmin+n)*wdm->df;

    // initialize noise levels
    for(int i=0; i<model->Nlink; i++)
    {
        model->soms[i] = 1.28e-22;
        model->sacc[i] = 5.76e-30;
    }
    
    // get noise covariance matrix for initial parameters
    generate_instrument_noise_model_wavelet(wdm, orbit, model);
}

void initialize_foreground_model(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model)
{
    // initialize data models
    alloc_foreground_model(model, data->NFFT, data->Nlayer, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = data->fmin + (double)n/data->T;

    // initialize constant foreground parameters levels 
    model->Tobs  =  data->T;
    model->Amp   =  1.2826e-44;
    model->alpha =  1.629667;
    model->f2    =  4.810781e-4;

    // initialize time dependent foreground parameter elves
    double af1 = -2.235e-1;
    double bf1 = -2.7040844;
    double afk = -3.60976122e-1;
    double bfk = -2.37822436;

    model->f1  =  pow(10., af1*log10(model->Tobs/YEAR) + bf1);
    model->fk  =  pow(10., afk*log10(model->Tobs/YEAR) + bfk);

    map_foreground_params_to_array(model);
    
    // get noise covariance matrix for initial parameters
    generate_galactic_foreground_model(model);

}

void default_sgwb_injection(double* params, const SGWB_t SGWB_type) {
    _Static_assert(SGWB_TEMPLATE_COUNT == 4, "Did you add an SGWB template? Edit this switch case, it needs to be exhaustive.");
    // set default values
    switch (SGWB_type) {
        case SGWB_TEMPLATE_POWERLAW:
            params[0] = -21.0;
            //params[0] = -8.45;
            //params[0] = -16.0;
            params[1] = 2./3.;
            //params[1] = 0.0;
            break;
        case SGWB_TEMPLATE_LOGNORMAL:
            // log10 A, log10 fstar [Hz], log10 Delta
            //params[0] = -5.0;
            params[0] = -2.5;
            params[1] = log(2e-3) / log(10.0);
            params[2] = 0.0;
            break;
        case SGWB_TEMPLATE_PHASE_TRANSITION:
            // rb, b, log10 Ap, log10 fp [Hz]
            params[0] = 1.0;
            params[1] = 1.0;
            params[2] = -18.0;
            params[3] = -2.5;
            break;
        case SGWB_TEMPLATE_LOGNORMAL_FIXEDF:
            // log10 A, log10 fstar [Hz], log10 Delta
            //params[0] = -5.0;
            params[0] = -2.5;
            break;
        default:
            fprintf(stderr,"need default values for SGWB type: %s", SGWB_TEMPLATE_NAMES[SGWB_type]);
            exit(1);
            break;
    }
}


void initialize_sgwb_model(struct Orbit *orbit, struct Data *data, struct SGWBModel *model, SGWB_t SGWB_type, const double *inj_params)
{
    // initialize data models
    alloc_sgwb_model(model, data->NFFT, data->Nlayer, data->Nchannel, SGWB_type);

    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = data->fmin + (double)n/data->T;

    if(inj_params) memcpy(model->params, inj_params, SGWB_TEMPLATE_NPARAMS[SGWB_type]*sizeof(double));
    else           default_sgwb_injection(model->params, SGWB_type);

    model->Tobs  =  data->T;
    // get covariance matrix for initial parameters
    generate_sgwb_model(model);

}

void initialize_sgwb_model_wavelet(struct Orbit *orbit, struct Data *data, struct SGWBModel *model, SGWB_t SGWB_type, const double *inj_params)
{
    struct Wavelets* wdm = data->wdm;

    // initialize data models
    alloc_sgwb_model(model, data->lmax - data->lmin, data->Nlayer, data->Nchannel, SGWB_type);

    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = (data->lmin+n)*wdm->df;

    if(inj_params) memcpy(model->params, inj_params, SGWB_TEMPLATE_NPARAMS[SGWB_type]*sizeof(double));
    else           default_sgwb_injection(model->params, SGWB_type);


    model->Tobs  =  data->T;
    // get covariance matrix for initial parameters
    generate_sgwb_model_wavelet(wdm, model);

}

void initialize_foreground_model_wavelet(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model)
{
    //wavelet basis
    struct Wavelets *wdm = data->wdm;

    // initialize data models
    alloc_foreground_model(model, data->lmax - data->lmin, data->Nlayer, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = (data->lmin+n)*wdm->df;

    // initialize constant foreground parameters levels 
    model->Tobs  =  data->T;
    model->Amp   =  1.2826e-44;
    //model->Amp   =  1.2826e-53;
    model->alpha =  1.629667;
    model->f2    =  4.810781e-4;

    // initialize time dependent foreground parameter elves
    double af1 = -2.235e-1;
    double bf1 = -2.7040844;
    double afk = -3.60976122e-1;
    double bfk = -2.37822436;

    model->f1  =  pow(10., af1*log10(model->Tobs/YEAR) + bf1);
    model->fk  =  pow(10., afk*log10(model->Tobs/YEAR) + bfk);
    map_foreground_params_to_array(model);
    
    // get noise covariance matrix for initial parameters
    generate_galactic_foreground_model_wavelet(wdm,model);

    // get galaxy modulation
    model->modulation = malloc(sizeof(struct GalaxyModulation));
    initialize_galaxy_modulation(model->modulation, data->wdm, orbit, data->T, data->t0);
    
    /**************************************************
     * Compute galaxy modulation
    **************************************************/
    // TODO: use Pozzoli et al's analytic response for modulation
    
    double *galaxy_params = double_vector(6); // defines galaxy shape
    galaxy_params[0] = 0.25; // A 0.25    bulge fraction
    galaxy_params[1] = 0.8;  // Rb 0.8    bulge radius (kpc)
    galaxy_params[2] = 2.5;  // Rd 2.5    disk radius (kpc)
    galaxy_params[3] = 0.4;  // Zd 0.4    disk height (kpc)
    galaxy_params[4] = 7.2;  // RGC 7.2   distance from solar BC to GC (kpc)
    galaxy_params[5] = 3.5;  // Rcut 3.5  radius out to which all sources are found (kpc)

    //computes the modulation of the confusion noise
    galaxy_modulation(model->modulation, galaxy_params);
    free_double_vector(galaxy_params);
}

void GetDynamicNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    /**************************************************
     * Compute instrument noise levels
    **************************************************/
    struct InstrumentModel *inst_noise = malloc(sizeof(struct InstrumentModel));
    initialize_instrument_model_wavelet(orbit, data, inst_noise);

    /**************************************************
     * Compute galactic foreground noise levels
    **************************************************/
    struct ForegroundModel *conf_noise = malloc(sizeof(struct ForegroundModel));
    initialize_foreground_model_wavelet(orbit, data, conf_noise);

    /**************************************************
     * Compute SGWB noise levels
    **************************************************/
    struct SGWBModel *sgwb = malloc(sizeof(struct SGWBModel));
    initialize_sgwb_model_wavelet(orbit, data, sgwb, flags->sgwbTemplate, flags->sgwbInjN ? flags->sgwbInjParams : NULL);

    /**************************************************
     * Combine noise components
    **************************************************/
    generate_full_dynamic_covariance_matrix(data->wdm, inst_noise, conf_noise, sgwb, data->noise);
    invert_noise_covariance_matrix(data->noise);

    char filename[128];
    sprintf(filename,"%s/power_noise.dat",data->dataDir);
    FILE *fptr=fopen(filename,"w");
    int k;
    for(int j=data->lmin; j<data->lmax; j++)
    {
        double f = j*data->wdm->df;
        for(int i=0; i<data->wdm->NT; i++)
        {
            double t = i*data->wdm->dt;
            wavelet_pixel_to_index(data->wdm,i,j,&k);
            k-=data->wdm->kmin;
            fprintf(fptr,"%lg %lg %.14e\n", t, f, data->noise->C[0][0][k]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);

    free_instrument_model(inst_noise);
    free_foreground_model(conf_noise);
    free_sgwb_model(sgwb);
}

void GetStationaryNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Noise *noise)
{
    /**************************************************
     * Compute instrument noise levels
    **************************************************/
    struct InstrumentModel *inst_noise = malloc(sizeof(struct InstrumentModel));
    initialize_instrument_model_wavelet(orbit, data, inst_noise);

    /**************************************************
     * Compute galactic foreground noise levels
    **************************************************/
    struct ForegroundModel *conf_noise = malloc(sizeof(struct ForegroundModel));
    initialize_foreground_model_wavelet(orbit, data, conf_noise);

    /**************************************************
     * Compute SGWB noise levels
    **************************************************/
    struct SGWBModel *sgwb = malloc(sizeof(struct SGWBModel));
    initialize_sgwb_model_wavelet(orbit, data, sgwb, flags->sgwbTemplate, flags->sgwbInjN ? flags->sgwbInjParams : NULL);


    /**************************************************
     * Combine noise components
    **************************************************/
    generate_full_stationary_covariance_matrix(data->wdm, inst_noise, conf_noise, sgwb, noise);
    invert_noise_covariance_matrix(noise);

    char filename[128];
    sprintf(filename,"%s/power_stationary_noise.dat",data->dataDir);
    FILE *fptr=fopen(filename,"w");
    int k;
    // TODO: check bounds here
    for(int j=data->lmin; j<data->lmax; j++)
    {
        double f = j*data->wdm->df;
        for(int i=0; i<data->wdm->NT; i++)
        {
            double t = i*data->wdm->dt;
            wavelet_pixel_to_index(data->wdm,i,j,&k);
            k-=data->wdm->kmin;
            fprintf(fptr,"%lg %lg %.14e\n", t, f,noise->C[0][0][k]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);

    free_instrument_model(inst_noise);
    free_foreground_model(conf_noise);
    free_sgwb_model(sgwb);
}
 
