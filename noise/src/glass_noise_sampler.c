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

void spline_ptmcmc(struct SplineModel **model, struct Chain *chain, struct Flags *flags)
{
    int a, b;
    int olda, oldb;
    
    double heat1, heat2;
    double logL1, logL2;
    double dlogL;
    double H;
    double alpha;
    double beta;
    
    int NC = chain->NC;
    
    for(b=NC-1; b>0; b--)
    {
        a = b - 1;
        chain->acceptance[a]=0;
        
        olda = chain->index[a];
        oldb = chain->index[b];
        
        heat1 = chain->temperature[a];
        heat2 = chain->temperature[b];
        
        logL1 = model[olda]->logL;
        logL2 = model[oldb]->logL;
        
        //Hot chains jump more rarely
        if(rand_r_U_0_1(&chain->r[a])<1.0)
        {
            dlogL = logL2 - logL1;
            H  = (heat2 - heat1)/(heat2*heat1);
            
            alpha = exp(dlogL*H);
            beta  = rand_r_U_0_1(&chain->r[a]);
            
            if(alpha >= beta)
            {
                chain->index[a] = oldb;
                chain->index[b] = olda;
                chain->acceptance[a]=1;
            }
        }
    }
}

// TODO: compare ptmcmc methods. Some people do pairwise random swaps

void noise_ptmcmc(struct InstrumentModel **model, struct Chain *chain, struct Flags *flags)
{
    int a, b;
    int olda, oldb;
    
    double heat1, heat2;
    double logL1, logL2;
    double dlogL;
    double H;
    double alpha;
    double beta;
    
    int NC = chain->NC;
    
    for(b=NC-1; b>0; b--)
    {
        a = b - 1;
        chain->acceptance[a]=0;
        
        olda = chain->index[a];
        oldb = chain->index[b];
        
        heat1 = chain->temperature[a];
        heat2 = chain->temperature[b];
        
        logL1 = model[olda]->logL;
        logL2 = model[oldb]->logL;
        
        //Hot chains jump more rarely
        if(rand_r_U_0_1(&chain->r[a])<1.0)
        {
            dlogL = logL2 - logL1;
            H  = (heat2 - heat1)/(heat2*heat1);
            
            alpha = exp(dlogL*H);
            beta  = rand_r_U_0_1(&chain->r[a]);
            
            if(alpha >= beta)
            {
                chain->index[a] = oldb;
                chain->index[b] = olda;
                chain->acceptance[a]=1;
            }
        }
    }
}


// TODO: this is log-uniform not uniform!
static double uniform_frequency_draw(double fmin, double fmax, unsigned int *r)
{
    return exp(log(fmin) + (log(fmax) - log(fmin))*rand_r_U_0_1(r));
}

static int check_frequency_spacing(double *f, int k, double T)
{
    if (
        (f[k]-f[k-1])*T < MIN_SPLINE_SPACING ||
        (f[k+1]-f[k])*T < MIN_SPLINE_SPACING
        ) return 1;
    else return 0;
}

void noise_spline_model_mcmc(struct Orbit *orbit, struct Data *data, struct SplineModel *model, struct Chain *chain, struct Flags *flags, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct SplineModel *model_x = model;
    struct SplineModel *model_y = malloc(sizeof(struct SplineModel));
    alloc_spline_model(model_y, model_x->psd->N, data->Nlayer, data->Nchannel, model_x->spline->N);

    //alisases for pointers to frequency vectors
    double *fx = model_x->spline->f;
    double *fy = model_y->spline->f;
    
    //copy current state into trial
    copy_spline_model(model_x, model_y);
    
    //pick a point any point
    int k = (int)floor(rand_r_U_0_1(&chain->r[ic])*(double)model_y->spline->N);

    //find the minimum and maximum indecies for stencil
    //int half_stencil = (MIN_SPLINE_STENCIL-1)/2;
    //int kmin = (k-half_stencil < 0) ? 0 : k-half_stencil;
    //int kmax = (k+half_stencil > model_x->spline->N-1) ? model_x->spline->N-1 : k+half_stencil;
    
    //update frequency
    model_y->spline->f[k] = model_x->spline->f[k];
    if(k>0 && k<model_y->spline->N-1)
    {
        /* normal draw */
        if(rand_r_U_0_1(&chain->r[ic])<0.8)
        {
            
            //get shortest distance between neighboring points
            double df_left = log(fx[k]) - log(fx[k-1]);
            double df_right= log(fx[k+1]) - log(fx[k]);
            double df = (df_left < df_right) ? df_left : df_right;
            
            //set shortest distance as 3-sigma jump for Gaussian draw
            double sigma = df/3.;
            
            //draw new frequency (illegally) checking that it stays between existing points
            do fy[k] = exp(log(fx[k]) + sigma*rand_r_N_0_1(&chain->r[ic]));
            while( fy[k]<fy[k-1] || fy[k]>fy[k+1]);
        }
        else
        {
            /* uniform draw */
            fy[k] = uniform_frequency_draw(fx[k-1], fx[k+1], &chain->r[ic]);
        }

        //check frequency prior
        if(check_frequency_spacing(fy, k, data->T)) logPy = -INFINITY;
        
    }

    //update amplitude
    double Sop, Spm;
    get_noise_levels("sangria", model_y->spline->f[k], &Spm, &Sop);
    double Sn;
    if(model_y->Nchannel==2) Sn = AEnoise_FF(orbit->L, orbit->fstar, model_y->spline->f[k], Spm, Sop);
    else Sn = XYZnoise_FF(orbit->L, orbit->fstar, model_y->spline->f[k], Spm, Sop);
    double scale = pow(10., -2.0 + 2.0*rand_r_U_0_1(&chain->r[ic]));
    for(int n=0; n<model_y->Nchannel; n++)
        model_y->spline->C[n][n][k] += scale*Sn*rand_r_N_0_1(&chain->r[ic]);
    

    
    //compute spline
    if(!flags->prior)
    {
        /* compute spline model */
        generate_spline_noise_model(model_y); //full interpolation
        //update_spline_noise_model(model_y, k, kmin, kmax); //interpolation over stencil
        
        /* get spline model likelihood */
        model_y->logL = noise_log_likelihood(data, model_y->psd);
        //model_y->logL = model_x->logL + noise_delta_log_likelihood(data, model_x, model_y, model_x->spline->f[kmin] , model_x->spline->f[kmax],ic);

        
        /*
         H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
         */
        logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
    }
    logH += logPy - logPx; //priors
    
    loga = log(rand_r_U_0_1(&chain->r[ic]));
    if(logH > loga)
    {
        copy_spline_model(model_y, model_x);
    }
    
    free_spline_model(model_y);
}

void noise_spline_model_rjmcmc(struct Orbit *orbit, struct Data *data, struct SplineModel *model, struct Chain *chain, struct Flags *flags, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct SplineModel *model_x = model;
    struct SplineModel *model_y = malloc(sizeof(struct SplineModel));
    
    //decide if doing a birth or death move
    int Nspline;
    char move;
    if(rand_r_U_0_1(&chain->r[ic])>0.5)
    {
        Nspline = model_x->spline->N + 1; //birth move
        move = 'B';
    }
    else
    {
        Nspline = model_x->spline->N - 1; //death move
        move = 'D';
    }
    alloc_spline_model(model_y, model_x->psd->N, data->Nlayer, data->Nchannel, Nspline);
    model_y->Nmin = model_x->Nmin;
    model_y->Nmax = model_x->Nmax;

    //alisases for pointers to frequency vectors
    double *fx = model_x->spline->f;
    double *fy = model_y->spline->f;
    
    copy_noise(model_x->psd,model_y->psd);

    //check range
    if(Nspline < model_x->Nmin || Nspline >= model_x->Nmax)
        move = 'R'; //reject move
    
    int kmin,kmax;
    switch(move)
    {
        case 'B':
            /*
             the move is to slot a new spline point between two existing
             points [kmin,kmax]
             */
            
            // first pick the left most point (it can't be the last point)
            kmin = (int)floor(rand_r_U_0_1(&chain->r[0])*(double)(model_x->spline->N-1));
            kmax = kmin+1;
            
            //copy current state into trial
            for(int k=0; k<=kmin; k++)
            {
                model_y->spline->f[k] = model_x->spline->f[k];
                for(int n=0; n<model_y->Nchannel; n++)
                {
                    model_y->spline->C[n][n][k] = model_x->spline->C[n][n][k];
                    model_y->spline->invC[n][n][k] = model_x->spline->invC[n][n][k];
                }
                model_y->spline->logdetC[k] = model_x->spline->logdetC[k];
            }
            
            //get grid place for new point
            int birth = kmin+1;
            fy[birth] = uniform_frequency_draw(fx[kmin], fx[kmax], &chain->r[ic]);
            
            //check frequency prior
            if(check_frequency_spacing(fy, birth, data->T)) logPy = -INFINITY;

            double Spm,Sop;
            get_noise_levels("sangria", model_y->spline->f[birth], &Spm, &Sop);

            double Sn;
            if(model_y->Nchannel==2) Sn = AEnoise_FF(orbit->L, orbit->fstar, model_y->spline->f[birth], Spm, Sop);
            else Sn = XYZnoise_FF(orbit->L, orbit->fstar, model_y->spline->f[birth], Spm, Sop);

            double Snmin = -Sn*100.;
            double Snmax =  Sn*100.;
            for(int n=0; n<model_y->Nchannel; n++)
                model_y->spline->C[n][n][birth] = Snmin + (Snmax - Snmin)*rand_r_U_0_1(&chain->r[ic]);
            
            // now fill in all higher points over k index
            for(int k=kmax; k<model_x->spline->N; k++)
            {
                model_y->spline->f[k+1] = model_x->spline->f[k];
                for(int n=0; n<model_y->Nchannel; n++)
                    model_y->spline->C[n][n][k+1] = model_x->spline->C[n][n][k];
            }
            break;

        case 'D':
            /*
             the move is to remove any existing point between the end points (0,Nspline)
             */
            
            // first pick the left most point (it can't be the last point)
            kmin = 1;
            kmax = model_x->spline->N - 1;
            int kill = kmin + (int)floor( (double)(kmax-kmin)*rand_r_U_0_1(&chain->r[ic]) );

            //copy current state into trial
            for(int k=0; k<kill; k++)
            {
                model_y->spline->f[k] = model_x->spline->f[k];
                for(int n=0; n<model_y->Nchannel; n++)
                    model_y->spline->C[n][n][k] = model_x->spline->C[n][n][k];
            }
            
            // now fill in all higher points over k index
            for(int k=kill; k<model_y->spline->N; k++)
            {
                model_y->spline->f[k] = model_x->spline->f[k+1];
                for(int n=0; n<model_y->Nchannel; n++)
                    model_y->spline->C[n][n][k] = model_x->spline->C[n][n][k+1];
            }
            
            break;
            
        case 'R':
            logPy = -INFINITY;
            break;
            
        default:
            printf("Invalid case %c\n",move);
    }
    
    fflush(stdout);
    
    //compute Hasting's ratio
    if(logPy > -INFINITY && !flags->prior)
    {

        generate_spline_noise_model(model_y);
        
        //compute likelihood
        model_y->logL = noise_log_likelihood(data, model_y->psd);
        
        /*
         H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
         */
        logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
    }
    logH += logPy - logPx; //priors
    
    loga = log(rand_r_U_0_1(&chain->r[ic]));
    if(logH > loga)
    {
        //reallocate noise structure to size of new spline model
        //realloc_noise(model_x->spline, model_y->spline->N);
        
        free_noise(model_x->spline);
        model_x->spline = malloc(sizeof(struct Noise));
        alloc_noise(model_x->spline,model_y->spline->N, data->Nlayer, data->Nchannel);
         
        
        copy_spline_model(model_y, model_x);
    }
    
    free_spline_model(model_y);
}

void noise_instrument_model_mcmc(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model, struct InstrumentModel *trial, struct ForegroundModel *galaxy, struct SGWBModel *sgwb, struct Noise *psd, struct Chain *chain, struct Flags *flags, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct InstrumentModel *model_x = model;
    struct InstrumentModel *model_y = trial;
    copy_instrument_model(model_x,model_y);

    //skip initialize likelihood
    /*
    generate_instrument_noise_model(orbit,model_x);
    copy_Cij(model_x->psd->C, psd->C, psd->Nchannel, psd->N);
    if(flags->confNoise) 
        generate_full_covariance_matrix(psd,galaxy->psd, data->Nchannel);
    if(flags->sgwbTemplate>=0) 
        generate_full_covariance_matrix(psd,sgwb->psd, data->Nchannel);
    invert_noise_covariance_matrix(psd);
    
    model_x->logL = noise_log_likelihood(data, psd);
    */
    
    
    //set priors
    double Sacc = 9.00e-30;
    double Soms = 2.25e-22;
    double Sacc_min = 9.00e-30/10;
    double Sacc_max = 9.00e-30*10;
    double Soms_min = 2.25e-22/10;
    double Soms_max = 2.25e-22*10;
    
    //set correlation matrix
    double *acc_jump_vec = malloc(model_x->Nlink*sizeof(double));
    double **correlation_matrix = malloc(model_x->Nlink*sizeof(double *));
    for(int n=0; n<model_x->Nlink; n++)
    {
        correlation_matrix[n] = malloc(model_x->Nlink*sizeof(double));
        correlation_matrix[n][n] = +1.0;
    }
    correlation_matrix[0][1] = correlation_matrix[1][0] = -1.0;
    correlation_matrix[2][3] = correlation_matrix[3][2] = -1.0;
    correlation_matrix[4][5] = correlation_matrix[5][4] = -1.0;

    correlation_matrix[0][2] = correlation_matrix[2][0] = +1.0;
    correlation_matrix[0][3] = correlation_matrix[3][0] = -1.0;
    correlation_matrix[0][4] = correlation_matrix[4][0] = -1.0;
    correlation_matrix[0][5] = correlation_matrix[5][0] = +1.0;

    correlation_matrix[1][2] = correlation_matrix[2][1] = -1.0;
    correlation_matrix[1][3] = correlation_matrix[3][1] = +1.0;
    correlation_matrix[1][4] = correlation_matrix[4][1] = +1.0;
    correlation_matrix[1][5] = correlation_matrix[5][1] = -1.0;

    correlation_matrix[2][4] = correlation_matrix[4][2] = -1.0;
    correlation_matrix[2][5] = correlation_matrix[5][2] = +1.0;

    correlation_matrix[3][4] = correlation_matrix[4][3] = +1.0;
    correlation_matrix[3][5] = correlation_matrix[5][3] = -1.0;

    for(int mc=0; mc<10; mc++)
    {
        //get jump sizes
        double acc_jump,oms_jump;
        double scale;
        if(rand_r_U_0_1(&chain->r[ic])>0.75)
            scale = 1;
        else if(rand_r_U_0_1(&chain->r[ic])>0.5)
            scale = 0.1;
        else if(rand_r_U_0_1(&chain->r[ic])>0.25)
            scale = 0.01;
        else
            scale = 0.001;
        
        /* get proposed noise parameters */
        
        int type;
        
        int i;
        int j;
        type = (int)(rand_r_U_0_1(&chain->r[ic])*3.);
        
        switch(type)
        {
            case 0:
                //update one link at a time
                i = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nlink);
                model_y->sacc[i] = model_x->sacc[i] + scale * Sacc * rand_r_N_0_1(&chain->r[ic]);
                model_y->soms[i] = model_x->soms[i] + scale * Soms * rand_r_N_0_1(&chain->r[ic]);
                
                //OMS noise is degenerate on a link
                if(i%2==0) model_y->soms[i+1] = model_y->soms[i];
                else model_y->soms[i-1] = model_y->soms[i];
                
                //check priors
                if(model_y->sacc[i] < Sacc_min || model_y->sacc[i] > Sacc_max) logPy = -INFINITY;
                if(model_y->soms[i] < Soms_min || model_y->soms[i] > Soms_max) logPy = -INFINITY;
                
                break;
            case 1:
                //update pair of links according to correlations
                i = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nlink);
                do
                {
                    j = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nlink);
                }while(i==j);
                
                acc_jump = scale * Sacc * rand_r_N_0_1(&chain->r[ic]);
                oms_jump = scale * Soms * rand_r_N_0_1(&chain->r[ic]);
                
                if(abs(i-j)==1)
                {
                    model_y->sacc[i] = model_x->sacc[i] + acc_jump;
                    model_y->sacc[j] = model_x->sacc[j] - acc_jump;
                    model_y->soms[i] = model_x->soms[i] + oms_jump;
                }
                else
                {
                    model_y->sacc[i] = model_x->sacc[i] + acc_jump;
                    model_y->sacc[j] = model_x->sacc[j] + acc_jump;
                    model_y->soms[i] = model_x->soms[i] + oms_jump;
                    model_y->soms[j] = model_x->soms[j] - oms_jump;
                }
                
                //OMS noise is degenerate on a link
                if(i%2==0) model_y->soms[i+1] = model_y->soms[i];
                else model_y->soms[i-1] = model_y->soms[i];
                
                if(j%2==0) model_y->soms[j+1] = model_y->soms[j];
                else model_y->soms[j-1] = model_y->soms[j];
                
                //check priors
                if(model_y->sacc[i] < Sacc_min || model_y->sacc[i] > Sacc_max) logPy = -INFINITY;
                if(model_y->sacc[j] < Sacc_min || model_y->sacc[j] > Sacc_max) logPy = -INFINITY;
                if(model_y->soms[i] < Soms_min || model_y->soms[i] > Soms_max) logPy = -INFINITY;
                if(model_y->soms[j] < Soms_min || model_y->soms[j] > Soms_max) logPy = -INFINITY;
                
                break;
            case 2:
                //update acceleration noise using correlation matrix
                for(i=0; i<model_x->Nlink; i++)
                {
                    acc_jump_vec[i] = scale * Sacc * rand_r_N_0_1(&chain->r[ic]);
                    model_y->sacc[i] = model_x->sacc[i];
                }
                
                
                for(i=0; i<model_x->Nlink; i++)
                {
                    for(j=0; j<model_x->Nlink; j++)
                    {
                        model_y->sacc[i] += correlation_matrix[i][j] * acc_jump_vec[j];
                    }
                    //check priors
                    if(model_y->sacc[i] < Sacc_min || model_y->sacc[i] > Sacc_max) logPy = -INFINITY;
                }
                break;
        }
        
        //get noise covariance matrix for initial parameters
        if(logPy > -INFINITY && !flags->prior)
        {
            generate_instrument_noise_model(orbit,model_y);
            copy_Cij(model_y->psd->C, psd->C, psd->Nchannel, psd->N);
            
            //add foreground noise contribution
            if(flags->confNoise)
                generate_full_covariance_matrix(psd,galaxy->psd, data->Nchannel);
            //add sgwb contribution
            if(flags->sgwbTemplate>=0) 
                generate_full_covariance_matrix(psd,sgwb->psd, data->Nchannel);
            
            invert_noise_covariance_matrix(psd);
            
            model_y->logL = my_noise_log_likelihood(data, psd);
            
            logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
        }
        logH += logPy - logPx; //priors
        
        loga = log(rand_r_U_0_1(&chain->r[ic]));
        if(logH > loga)
        {
            copy_instrument_model(model_y, model_x);
            if(flags->confNoise) galaxy->logL = model_x->logL;
            if(flags->sgwbTemplate>=0) sgwb->logL = model_x->logL;
        }
    }

    free(acc_jump_vec);
    for(int n=0; n<model_x->Nlink; n++)
        free(correlation_matrix[n]);
    free(correlation_matrix);

}

void noise_foreground_model_mcmc(struct Data *data, struct InstrumentModel *noise, struct ForegroundModel *model, struct ForegroundModel *trial, struct SGWBModel *sgwb, struct Noise *psd, struct Chain *chain, struct Flags *flags, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct ForegroundModel *model_x = model;
    struct ForegroundModel *model_y = trial;
    copy_foreground_model(model_x,model_y);
    
    //initialize likelhood
    // skip initialization
    /*
    generate_galactic_foreground_model(model_x);
    copy_Cij(model_x->psd->C, psd->C, psd->Nchannel, psd->N);
    generate_full_covariance_matrix(psd, noise->psd, data->Nchannel);
    if(flags->sgwbTemplate>=0) 
        generate_full_covariance_matrix(psd,sgwb->psd, data->Nchannel);
    invert_noise_covariance_matrix(psd);
    model_x->logL = noise_log_likelihood(data, psd);
    */

    
    //set priors
    double **prior = malloc(sizeof(double *)*model_x->Nparams);
    for(int n=0; n<model_x->Nparams; n++)
        prior[n] = malloc(sizeof(double)*2);
    
    //set correlation matrix
    double *acc_jump_vec = malloc(model_x->Nparams*sizeof(double));
    double **correlation_matrix = malloc(model_x->Nparams*sizeof(double *));
    for(int n=0; n<model_x->Nparams; n++)
    {
        correlation_matrix[n] = malloc(model_x->Nparams*sizeof(double));
        correlation_matrix[n][n] = +1.0;
    }
    correlation_matrix[0][1] = correlation_matrix[1][0] = -1.0;
    correlation_matrix[0][2] = correlation_matrix[2][0] = -1.0;
    correlation_matrix[0][3] = correlation_matrix[3][0] = -1.0;
    correlation_matrix[0][4] = correlation_matrix[4][0] = +1.0;

    correlation_matrix[1][2] = correlation_matrix[2][1] = +1.0;
    correlation_matrix[1][3] = correlation_matrix[3][1] = +1.0;
    correlation_matrix[1][4] = correlation_matrix[4][1] = -1.0;

    correlation_matrix[2][3] = correlation_matrix[3][2] = +1.0;
    correlation_matrix[2][4] = correlation_matrix[4][2] = -1.0;

    correlation_matrix[3][4] = correlation_matrix[4][3] = -1.0;
    
    
    //log(A)
    prior[0][0] = -102.0;
    prior[0][1] = -99.0;
    
    //f1
    prior[1][0] = log(0.0001);
    prior[1][1] = log(0.01);
    
    //alpha
    prior[2][0] = 0.0;
    prior[2][1] = 3.0;
    
    //fk
    prior[3][0] = log(0.0001);
    prior[3][1] = log(0.1);
    
    //f2
    prior[4][0] = log(0.00001);
    prior[4][1] = log(0.01);
        
    for(int mc=0; mc<10; mc++)
    {
        
        /* get proposed noise parameters */
        
        //get jump sizes
        double scale;
        if(rand_r_U_0_1(&chain->r[ic])>0.75)
            scale = 0.1;
        else if(rand_r_U_0_1(&chain->r[ic])>0.5)
            scale = 0.01;
        else if(rand_r_U_0_1(&chain->r[ic])>0.25)
            scale = 0.001;
        else
            scale = 0.0001;
        
        if(rand_r_U_0_1(&chain->r[ic])<0.5)
        {
            //pick which parameter to update
            int i = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nparams);
            model_y->sgal[i] = model_x->sgal[i] + scale * 0.5*(prior[i][1]+prior[i][0]) * rand_r_N_0_1(&chain->r[ic]);
        }
        else
        {
            //jump along correlated directions
            double jump=rand_r_N_0_1(&chain->r[ic]);
            for(int n=0; n<model_x->Nparams; n++)
            {
                acc_jump_vec[n] = scale *  0.5*(prior[n][1]+prior[n][0]) * jump;
                model_y->sgal[n] = model_x->sgal[n];
            }
            
            //pick which vector
            int i = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nparams);
            
            //jump
            for(int j=0; j<model_x->Nparams; j++)
                model_y->sgal[j] += correlation_matrix[i][j] * acc_jump_vec[j];

        }
        
        //check priors
        for(int n=0; n<model_y->Nparams; n++)
            if(model_y->sgal[n] < prior[n][0] || model_y->sgal[n] > prior[n][1]) {
                logPy = -INFINITY;
                break;
            }
        
        //get noise covariance matrix for initial parameters
        if(logPy > -INFINITY && !flags->prior)
        {
            generate_galactic_foreground_model(model_y);
            copy_Cij(model_y->psd->C, psd->C, psd->Nchannel, psd->N);
            
            //add instrument noise contribution
            generate_full_covariance_matrix(psd, noise->psd, data->Nchannel);
            //add sgwb contribution
            if(flags->sgwbTemplate>=0) 
                generate_full_covariance_matrix(psd, sgwb->psd, data->Nchannel);
            
            invert_noise_covariance_matrix(psd);
            
            model_y->logL = my_noise_log_likelihood(data, psd);
            
            logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
        }
        logH += logPy - logPx; //priors
        
        loga = log(rand_r_U_0_1(&chain->r[ic]));
        if(logH > loga)
        {
            copy_foreground_model(model_y, model_x);
            noise->logL = model_x->logL;
            if (sgwb)
                sgwb->logL = model_x->logL;
        }
    }
 
    free(acc_jump_vec);   
    for(int n=0; n<model_x->Nparams; n++) free(prior[n]);
    free(prior);
    for(int n=0; n<model_x->Nparams; n++)
        free(correlation_matrix[n]);
    free(correlation_matrix);
}

void noise_sgwb_model_mcmc(struct Data *data, struct InstrumentModel *noise, struct ForegroundModel *galaxy, struct SGWBModel *model, struct SGWBModel *trial, struct Noise *psd, struct Chain *chain, struct Flags *flags, int ic)
{

    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct SGWBModel *model_x = model;
    struct SGWBModel *model_y = trial;
    copy_sgwb_model(model_x,model_y);
    //printf("copied first sgwb\n");
    
    //initialize likelhood
    // skip initializing
    /*
    generate_sgwb_model(model_x);
    copy_Cij(model_x->psd->C, psd->C, psd->Nchannel, psd->N);
    generate_full_covariance_matrix(psd, noise->psd, data->Nchannel);
    if(flags->confNoise) 
        generate_full_covariance_matrix(psd,galaxy->psd, data->Nchannel);
    invert_noise_covariance_matrix(psd);
    model_x->logL = noise_log_likelihood(data, psd);
    //printf("likelihood init\n");
    */

    
    //set priors
    const double (*prior)[2] = default_sgwb_priors[model_x->SGWB_type];
    
    //set correlation matrix
    double acc_jump_vec[model_x->Nparams];
    double correlation_matrix[model_x->Nparams][model_x->Nparams];
    memset(correlation_matrix, 0, sizeof(double)*model_x->Nparams*model_x->Nparams);
    for(int n=0; n<model_x->Nparams; n++)
    {
        correlation_matrix[n][n] = +1.0;
    }
    
    for(int mc=0; mc<10; mc++)
    {
        
        /* get proposed noise parameters */
        
        //get jump sizes
        // TODO: this looks statistically sus
        double scale;
        if(rand_r_U_0_1(&chain->r[ic])>0.75)
            scale = 0.1;
        else if(rand_r_U_0_1(&chain->r[ic])>0.5)
            scale = 0.01;
        else if(rand_r_U_0_1(&chain->r[ic])>0.25)
            scale = 0.001;
        else
            scale = 0.0001;
        
        if(rand_r_U_0_1(&chain->r[ic])<0.5)
        {
            //pick which parameter to update
            int i = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nparams);
            model_y->params[i] = model_x->params[i] + scale * 0.5*(prior[i][1]+prior[i][0]) * rand_r_N_0_1(&chain->r[ic]);
        }
        else
        {
            //jump along correlated directions
            double jump = rand_r_N_0_1(&chain->r[ic]);
            for(int n=0; n<model_x->Nparams; n++)
            {
                acc_jump_vec[n] = scale *  0.5*(prior[n][1]+prior[n][0]) * jump;
                model_y->params[n] = model_x->params[n];
            }
            
            //pick which vector
            int i = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nparams);
            
            //jump
            for(int j=0; j<model_x->Nparams; j++)
                model_y->params[j] += correlation_matrix[i][j] * acc_jump_vec[j];

        }
        //printf("got proposal %g %g\n", model_y->params[0], model_y->params[1]);
        
        //check priors
        for(int n=0; n<model_y->Nparams; n++)
            if(model_y->params[n] < prior[n][0] || model_y->params[n] > prior[n][1])
            {
                logPy = -INFINITY;
                break;
            }
        
        //get noise covariance matrix for initial parameters
        if(logPy > -INFINITY && !flags->prior)
        {
            generate_sgwb_model(model_y);
            copy_Cij(model_y->psd->C, psd->C, psd->Nchannel, psd->N);
            
            //add instrument noise contribution
            generate_full_covariance_matrix(psd, noise->psd, data->Nchannel);
            //add confusion noise contribution
            if(flags->confNoise) 
                generate_full_covariance_matrix(psd,galaxy->psd, data->Nchannel);
            
            invert_noise_covariance_matrix(psd);
            
            model_y->logL = my_noise_log_likelihood(data, psd);
            //printf("trial %g %g %g\n", model_y->params[0], model_y->params[1], model_y->logL);
            
            logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
        }
        logH += logPy - logPx; //priors
        
        loga = log(rand_r_U_0_1(&chain->r[ic]));
        if(logH > loga)
        {
            //printf("accepted %g %g %g\n", model_x->params[0], model_x->params[1], model_x->logL);
            copy_sgwb_model(model_y, model_x);
            noise->logL = model_x->logL;
            galaxy->logL = model_x->logL;
        }
    }
 
    free(acc_jump_vec);   
    for(int n=0; n<model_x->Nparams; n++) free(prior[n]);
    free(prior);
    for(int n=0; n<model_x->Nparams; n++)
        free(correlation_matrix[n]);
    free(correlation_matrix);
}
void noise_sgwb_model_mcmc_wavelet_dumb(struct Data *data, struct InstrumentModel *noise, struct ForegroundModel *galaxy, struct SGWBModel *model, struct SGWBModel *trial, struct Noise *psd, const struct CoarseStats *stats, struct Chain *chain, struct Flags *flags, int ic)
{
    //shorthand pointers
    struct SGWBModel *model_x = model;
    struct SGWBModel *model_y = trial;
    copy_sgwb_model(model_x, model_y);

    // skip initializing likelihood

    //set priors
    const double (*prior)[2] = default_sgwb_priors[model_x->SGWB_type];
    double u = rand_r_U_0_1(&chain->r[ic]);
    double scale = 1.0;
    if (u < 0.10) {
        scale = 4;
    }
    else if (u < 0.20) {
        scale = 1e0;
    }
    else if (u < 0.40) {
        scale = 1e-1;
    }
    else if (u < 0.50) {
        scale = 1e-2;
    }
    else if (u < 0.70) {
        scale = 1e-3;
    }
    else {
        scale = 1e-5;
    }
    for (size_t i=0; i<model_y->Nparams; i++) {
        double r = rand_r_N_0_1(&chain->r[ic]);
        //printf("%lf\n",r);
        double prior_extent = prior[i][1] - prior[i][0];
        // choose jump size
        double jump = model_y->params[i] + scale*r/12.*prior_extent;
        model_y->params[i] = jump;
        if (jump > prior[i][1] || jump < prior[i][0]) {
            // try again if out of prior
            return;
        }
    }
    generate_sgwb_model_wavelet(data->wdm, model_y);
    if (flags->stationary)
        generate_full_stationary_covariance_matrix_coarse(data->wdm, stats->Q, noise, galaxy, model_y, psd);
    else
        generate_full_dynamic_covariance_matrix_coarse(data->wdm, stats->Q, noise, galaxy, model_y, psd);
    invert_noise_covariance_matrix(psd);
    model_y->logL = my_noise_log_likelihood_wavelet_coarse(data, psd, stats);

    double logH = (model_y->logL - model_x->logL)/chain->temperature[ic];
    double loga = log(rand_r_U_0_1(&chain->r[ic]));
    if (logH > loga) {
        /*
        printf("accepted %lf, %lf\n", model_y->params[0], model_y->params[1]);
        printf("\told logL: %g\n", model_x->logL);
        printf("\tnew logL: %g\n", model_y->logL);
        printf("\tscale %g\n", scale);
        */
        copy_sgwb_model(model_y, model_x);
        noise->logL = model_x->logL;
        if (galaxy)
            galaxy->logL = model_x->logL;
    }
    else {
        /*
        printf("rejected %lf, %lf\n", model_y->params[0], model_y->params[1]);
        printf("\told logL: %g\n", model_x->logL);
        printf("\tnew logL: %g\n", model_y->logL);
        printf("\tscale %g\n", scale);
        */
    }
}

void noise_instrument_model_mcmc_wavelet(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model, struct InstrumentModel *trial, struct ForegroundModel *galaxy, struct SGWBModel *sgwb, struct Noise *psd, const struct CoarseStats *stats, struct Chain *chain, struct Flags *flags, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct InstrumentModel *model_x = model;
    struct InstrumentModel *model_y = trial;
    copy_instrument_model(model_x,model_y);

    // skip initializing likelihood

    //set priors
    double Sacc = 9.00e-30;
    double Soms = 2.25e-22;
    double Sacc_min = 9.00e-30/10;
    double Sacc_max = 9.00e-30*10;
    double Soms_min = 2.25e-22/10;
    double Soms_max = 2.25e-22*10;

    //linear-drift priors (fractional change at the run endpoints; u in [-1,1])
    double Sdrift = 0.1;        //!< typical drift jump scale
    double drift_min = -1.0;
    double drift_max =  1.0;
    int sample_drift = model_x->drift_enabled;

    //set correlation matrix
    double *acc_jump_vec = malloc(model_x->Nlink*sizeof(double));
    double **correlation_matrix = malloc(model_x->Nlink*sizeof(double *));
    for(int n=0; n<model_x->Nlink; n++)
    {
        correlation_matrix[n] = malloc(model_x->Nlink*sizeof(double));
        correlation_matrix[n][n] = +1.0;
    }
    correlation_matrix[0][1] = correlation_matrix[1][0] = -1.0;
    correlation_matrix[2][3] = correlation_matrix[3][2] = -1.0;
    correlation_matrix[4][5] = correlation_matrix[5][4] = -1.0;

    correlation_matrix[0][2] = correlation_matrix[2][0] = +1.0;
    correlation_matrix[0][3] = correlation_matrix[3][0] = -1.0;
    correlation_matrix[0][4] = correlation_matrix[4][0] = -1.0;
    correlation_matrix[0][5] = correlation_matrix[5][0] = +1.0;

    correlation_matrix[1][2] = correlation_matrix[2][1] = -1.0;
    correlation_matrix[1][3] = correlation_matrix[3][1] = +1.0;
    correlation_matrix[1][4] = correlation_matrix[4][1] = +1.0;
    correlation_matrix[1][5] = correlation_matrix[5][1] = -1.0;

    correlation_matrix[2][4] = correlation_matrix[4][2] = -1.0;
    correlation_matrix[2][5] = correlation_matrix[5][2] = +1.0;

    correlation_matrix[3][4] = correlation_matrix[4][3] = +1.0;
    correlation_matrix[3][5] = correlation_matrix[5][3] = -1.0;

    for(int mc=0; mc<10; mc++)
    {
        logH = 0.0;
        //get jump sizes
        double acc_jump,oms_jump;
        double scale;
        if(rand_r_U_0_1(&chain->r[ic])>0.75)
            scale = 1;
        else if(rand_r_U_0_1(&chain->r[ic])>0.5)
            scale = 0.1;
        else if(rand_r_U_0_1(&chain->r[ic])>0.25)
            scale = 0.01;
        else
            scale = 0.001;
        
        /* get proposed noise parameters */
        
        int type;
        
        int i;
        int j;
        //when drift is sampled, a 4th proposal type perturbs the per-link drifts
        type = (int)(rand_r_U_0_1(&chain->r[ic])*(sample_drift ? 4. : 3.));
        
        switch(type)
        {
            case 0:
                //update one link at a time
                i = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nlink);
                model_y->sacc[i] = model_x->sacc[i] + scale * Sacc * rand_r_N_0_1(&chain->r[ic]);
                model_y->soms[i] = model_x->soms[i] + scale * Soms * rand_r_N_0_1(&chain->r[ic]);
                
                //OMS noise is degenerate on a link
                if(i%2==0) model_y->soms[i+1] = model_y->soms[i];
                else model_y->soms[i-1] = model_y->soms[i];
                
                //check priors
                if(model_y->sacc[i] < Sacc_min || model_y->sacc[i] > Sacc_max) logPy = -INFINITY;
                if(model_y->soms[i] < Soms_min || model_y->soms[i] > Soms_max) logPy = -INFINITY;
                
                break;
            case 1:
                //update pair of links according to correlations
                i = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nlink);
                do
                {
                    j = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nlink);
                }while(i==j);
                
                acc_jump = scale * Sacc * rand_r_N_0_1(&chain->r[ic]);
                oms_jump = scale * Soms * rand_r_N_0_1(&chain->r[ic]);
                
                if(abs(i-j)==1)
                {
                    model_y->sacc[i] = model_x->sacc[i] + acc_jump;
                    model_y->sacc[j] = model_x->sacc[j] - acc_jump;
                    model_y->soms[i] = model_x->soms[i] + oms_jump;
                }
                else
                {
                    model_y->sacc[i] = model_x->sacc[i] + acc_jump;
                    model_y->sacc[j] = model_x->sacc[j] + acc_jump;
                    model_y->soms[i] = model_x->soms[i] + oms_jump;
                    model_y->soms[j] = model_x->soms[j] - oms_jump;
                }
                
                //OMS noise is degenerate on a link
                if(i%2==0) model_y->soms[i+1] = model_y->soms[i];
                else model_y->soms[i-1] = model_y->soms[i];
                
                if(j%2==0) model_y->soms[j+1] = model_y->soms[j];
                else model_y->soms[j-1] = model_y->soms[j];
                
                //check priors
                if(model_y->sacc[i] < Sacc_min || model_y->sacc[i] > Sacc_max) logPy = -INFINITY;
                if(model_y->sacc[j] < Sacc_min || model_y->sacc[j] > Sacc_max) logPy = -INFINITY;
                if(model_y->soms[i] < Soms_min || model_y->soms[i] > Soms_max) logPy = -INFINITY;
                if(model_y->soms[j] < Soms_min || model_y->soms[j] > Soms_max) logPy = -INFINITY;
                
                break;
            case 2:
                //update acceleration noise using correlation matrix
                for(i=0; i<model_x->Nlink; i++)
                {
                    acc_jump_vec[i] = scale * Sacc * rand_r_N_0_1(&chain->r[ic]);
                    model_y->sacc[i] = model_x->sacc[i];
                }
                
                
                for(i=0; i<model_x->Nlink; i++)
                {
                    for(j=0; j<model_x->Nlink; j++)
                    {
                        model_y->sacc[i] += correlation_matrix[i][j] * acc_jump_vec[j];
                    }
                    //check priors
                    if(model_y->sacc[i] < Sacc_min || model_y->sacc[i] > Sacc_max) logPy = -INFINITY;
                }
                break;
            case 3:
                //update the per-link linear drift (only reached when sample_drift)
                i = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nlink);
                model_y->sacc_drift[i] = model_x->sacc_drift[i] + scale * Sdrift * rand_r_N_0_1(&chain->r[ic]);
                model_y->soms_drift[i] = model_x->soms_drift[i] + scale * Sdrift * rand_r_N_0_1(&chain->r[ic]);

                //OMS noise (and hence its drift) is degenerate on a link
                if(i%2==0) model_y->soms_drift[i+1] = model_y->soms_drift[i];
                else model_y->soms_drift[i-1] = model_y->soms_drift[i];

                //check priors
                if(model_y->sacc_drift[i] < drift_min || model_y->sacc_drift[i] > drift_max) logPy = -INFINITY;
                if(model_y->soms_drift[i] < drift_min || model_y->soms_drift[i] > drift_max) logPy = -INFINITY;

                break;
        }

        //get noise covariance matrix for initial parameters
        if(logPy > -INFINITY && !flags->prior)
        {
            generate_instrument_noise_model_wavelet(data->wdm, orbit, model_y);
            if (flags->stationary)
                generate_full_stationary_covariance_matrix_coarse(data->wdm, stats->Q, model_y, galaxy, sgwb, psd);
            else
                generate_full_dynamic_covariance_matrix_coarse(data->wdm, stats->Q, model_y, galaxy, sgwb, psd);
            invert_noise_covariance_matrix(psd);

            model_y->logL = my_noise_log_likelihood_wavelet_coarse(data, psd, stats);

            logH = (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
        }
        logH += logPy - logPx; //priors

        loga = log(rand_r_U_0_1(&chain->r[ic]));
        if(logH > loga)
        {
            copy_instrument_model(model_y, model_x);
            if (sgwb)
                sgwb->logL   = model_y->logL;
            if (galaxy)
                galaxy->logL = model_y->logL;
        }
    }

    free(acc_jump_vec);
    for(int n=0; n<model_x->Nlink; n++)
        free(correlation_matrix[n]);
    free(correlation_matrix);

}

void noise_foreground_model_mcmc_wavelet(struct Data *data, struct InstrumentModel *noise, struct ForegroundModel *model, struct ForegroundModel *trial, struct SGWBModel *sgwb, struct Noise *psd, const struct CoarseStats *stats, struct Chain *chain, struct Flags *flags, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct ForegroundModel *model_x = model;
    struct ForegroundModel *model_y = trial;
    copy_foreground_model(model_x,model_y);
    
    //skip initializing likelhood
    
    //set priors
    double **prior = malloc(sizeof(double *)*model_x->Nparams);
    for(int n=0; n<model_x->Nparams; n++)
        prior[n] = malloc(sizeof(double)*2);
    
    //set correlation matrix
    double *acc_jump_vec = malloc(model_x->Nparams*sizeof(double));
    double **correlation_matrix = malloc(model_x->Nparams*sizeof(double *));
    for(int n=0; n<model_x->Nparams; n++)
    {
        correlation_matrix[n] = malloc(model_x->Nparams*sizeof(double));
        correlation_matrix[n][n] = +1.0;
    }
    correlation_matrix[0][1] = correlation_matrix[1][0] = -1.0;
    correlation_matrix[0][2] = correlation_matrix[2][0] = -1.0;
    correlation_matrix[0][3] = correlation_matrix[3][0] = -1.0;
    correlation_matrix[0][4] = correlation_matrix[4][0] = +1.0;

    correlation_matrix[1][2] = correlation_matrix[2][1] = +1.0;
    correlation_matrix[1][3] = correlation_matrix[3][1] = +1.0;
    correlation_matrix[1][4] = correlation_matrix[4][1] = -1.0;

    correlation_matrix[2][3] = correlation_matrix[3][2] = +1.0;
    correlation_matrix[2][4] = correlation_matrix[4][2] = -1.0;

    correlation_matrix[3][4] = correlation_matrix[4][3] = -1.0;
    
    
    //log(A)
    //prior[0][0] = -86.0;
    //prior[0][1] = -82.0;
    prior[0][0] = -110.0;
    prior[0][1] = -80.0;
    
    //f1
    prior[1][0] = log(0.0001);
    prior[1][1] = log(0.01);
    
    //alpha
    prior[2][0] = 0.0;
    prior[2][1] = 3.0;
    
    //fk
    prior[3][0] = log(0.0001);
    prior[3][1] = log(0.1);
    
    //f2
    prior[4][0] = log(0.00001);
    prior[4][1] = log(0.01);
        
    for(int mc=0; mc<10; mc++)
    {
        logH = 0.0;
        
        /* get proposed noise parameters */
        
        //get jump sizes
        double scale;
        if(rand_r_U_0_1(&chain->r[ic])>0.75)
            scale = 0.1;
        else if(rand_r_U_0_1(&chain->r[ic])>0.5)
            scale = 0.01;
        else if(rand_r_U_0_1(&chain->r[ic])>0.25)
            scale = 0.001;
        else
            scale = 0.0001;
        
        if(rand_r_U_0_1(&chain->r[ic])<0.5)
        {
            //pick which parameter to update
            int i = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nparams);
            model_y->sgal[i] = model_x->sgal[i] + scale * 0.5*(prior[i][1]+prior[i][0]) * rand_r_N_0_1(&chain->r[ic]);
        }
        else
        {
            //jump along correlated directions
            double jump=rand_r_N_0_1(&chain->r[ic]);
            for(int n=0; n<model_x->Nparams; n++)
            {
                acc_jump_vec[n] = scale *  0.5*(prior[n][1]+prior[n][0]) * jump;
                model_y->sgal[n] = model_x->sgal[n];
            }
            
            //pick which vector
            int i = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nparams);
            
            //jump
            for(int j=0; j<model_x->Nparams; j++)
                model_y->sgal[j] += correlation_matrix[i][j] * acc_jump_vec[j];

        }
        
        //check priors
        for(int n=0; n<model_y->Nparams; n++)
            if(model_y->sgal[n] < prior[n][0] || model_y->sgal[n] > prior[n][1])
                logPy = -INFINITY;
        
        //get noise covariance matrix for initial parameters
        if(logPy > -INFINITY && !flags->prior)
        {
            map_array_to_foreground_params(model_y);
            generate_galactic_foreground_model_wavelet(data->wdm, model_y);
            if (flags->stationary)
                generate_full_stationary_covariance_matrix_coarse(data->wdm, stats->Q, noise, model_y, sgwb, psd);
            else
                generate_full_dynamic_covariance_matrix_coarse(data->wdm, stats->Q, noise, model_y, sgwb, psd);
            invert_noise_covariance_matrix(psd);

            model_y->logL = my_noise_log_likelihood_wavelet_coarse(data, psd, stats);

            logH = (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
        }
        logH += logPy - logPx; //priors
        
        loga = log(rand_r_U_0_1(&chain->r[ic]));
        if(logH > loga)
        {
            copy_foreground_model(model_y, model_x);
            noise->logL = model_x->logL;
            if (sgwb)
                sgwb->logL = model_x->logL;
        }
    }
 
    free(acc_jump_vec);   
    for(int n=0; n<model_x->Nparams; n++) free(prior[n]);
    free(prior);
    for(int n=0; n<model_x->Nparams; n++)
        free(correlation_matrix[n]);
    free(correlation_matrix);
}

void noise_sgwb_model_mcmc_wavelet(struct Data *data, struct InstrumentModel *noise, struct ForegroundModel *galaxy, struct SGWBModel *model, struct SGWBModel *trial, struct Noise *psd, struct Chain *chain, struct Flags *flags, int ic)
{

    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct SGWBModel *model_x = model;
    struct SGWBModel *model_y = trial;
    copy_sgwb_model(model_x, model_y);
    
    //skip initializing likelhood
    
    //set priors
    const double (*prior)[2] = default_sgwb_priors[model_x->SGWB_type];
    
    //set correlation matrix
    double acc_jump_vec[model_x->Nparams];
    double correlation_matrix[model_x->Nparams][model_x->Nparams];
    memset(correlation_matrix, 0, sizeof(double)*model_x->Nparams*model_x->Nparams);
    for(int n=0; n<model_x->Nparams; n++)
    {
        correlation_matrix[n][n] = +1.0;
    }
        
    for(int mc=0; mc<10; mc++)
    {
        logH = 0.0;
        
        /* get proposed noise parameters */
        
        //get jump sizes
        double scale;
        double u = rand_r_U_0_1(&chain->r[ic]);
        if(u>0.75)
            scale = 0.1;
        else if(u>0.5)
            scale = 0.01;
        else if(u>0.25)
            scale = 0.001;
        else
            scale = 0.0001;
        
        if(rand_r_U_0_1(&chain->r[ic])<0.5)
        {
            //pick which parameter to update
            int i = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nparams);
            model_y->params[i] = model_x->params[i] + scale * 0.5*(prior[i][1]+prior[i][0]) * rand_r_N_0_1(&chain->r[ic]);
        }
        else
        {
            //jump along correlated directions
            double jump = rand_r_N_0_1(&chain->r[ic]);
            for(int n=0; n<model_x->Nparams; n++)
            {
                acc_jump_vec[n] = scale *  0.5*(prior[n][1]+prior[n][0]) * jump;
                model_y->params[n] = model_x->params[n];
            }
            
            //pick which vector
            int i = (int)(rand_r_U_0_1(&chain->r[ic])* (double)model_x->Nparams);
            
            //jump
            for(int j=0; j<model_x->Nparams; j++)
                model_y->params[j] += correlation_matrix[i][j] * acc_jump_vec[j];

        }
        
        //check priors
        for(int n=0; n<model_y->Nparams; n++)
            if(model_y->params[n] < prior[n][0] || model_y->params[n] > prior[n][1])
            {
                logPy = -INFINITY;
                break;
            }
        
        //get noise covariance matrix for initial parameters
        if(logPy > -INFINITY && !flags->prior)
        {
            generate_sgwb_model_wavelet(data->wdm, model_y);
            
            //add other stochastic contributions
            if (flags->stationary)
                generate_full_stationary_covariance_matrix(data->wdm, noise, galaxy, model_y, psd);
            else
                generate_full_dynamic_covariance_matrix(data->wdm, noise, galaxy, model_y, psd);
            invert_noise_covariance_matrix(psd);
            
            model_y->logL = my_noise_log_likelihood_wavelet(data, psd);
            //printf("trial %g %g %g\n", model_y->params[0], model_y->params[1], model_y->logL);
            
            logH = (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
        }
        logH += logPy - logPx; //priors
        
        loga = log(rand_r_U_0_1(&chain->r[ic]));
        if(logH > loga)
        {
            //printf("accepted %g %g %g\n", model_x->params[0], model_x->params[1], model_x->logL);
            copy_sgwb_model(model_y, model_x);
            noise->logL = model_x->logL;
            if (galaxy)
                galaxy->logL = model_x->logL;
        }
    }
}
