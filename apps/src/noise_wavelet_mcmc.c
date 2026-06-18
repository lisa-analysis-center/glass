/*
 *  Copyright (C) 2024 Robert Rosati MSFC ST-12 / UAH
 *
 */

/**
 @file noise_wavelet_mcmc.c
 \brief Main function for wavelet-domain noise and SGWB sampler app `noise_wavelet_mcmc` 
 */

/*  REQUIRED LIBRARIES  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <sys/stat.h>

#include <omp.h>

#include <glass_utils.h>
#include <glass_noise.h>
#include <glass_noise_sampler.h>


static void print_usage()
{
    print_glass_usage();
    fprintf(stdout,"EXAMPLE:\n");
    fprintf(stdout,"noise_wavelet_mcmc --sim-noise --conf-noise --sgwb-template 0 --duration 7864320 --fmin 1e-4 --fmax 8e-3 [--coarse-Q 4]\n");
    fprintf(stdout,"\n");
    exit(0);
}

int write_time_data(double dt, int N, double* data, char* fname) {
    FILE *fptr = NULL;
    fptr = fopen(fname,"w");
    if (!fptr) {
        fprintf(stderr, "Couldn't open %s for writing!\n", fname);
        return 1;
    }
    for (int i=0; i<N; i++) {
            double f = i*dt;
            fprintf(fptr,"%lg ",f);
            fprintf(fptr,"%lg" ,data[i]);
            fprintf(fptr,"\n");
    }
    fclose(fptr);
    return 0;
}

int write_fft_data(double df, int NFFT, double* data, char* fname) {
    FILE *fptr = NULL;
    fptr = fopen(fname,"w");
    if (!fptr) {
        fprintf(stderr, "Couldn't open %s for writing!\n", fname);
        return 1;
    }
    for (int i=0; i<NFFT; i++) {
            double f = i*df;
            fprintf(fptr,"%lg ",f);
            fprintf(fptr,"%lg ",data[2*i]);
            fprintf(fptr,"%lg" ,data[2*i+1]);
            fprintf(fptr,"\n");
    }
    fclose(fptr);
    return 0;
}

int write_wdm_data(struct Data* d, double* data, char* fname) {
    FILE *fptr = NULL;
    struct Wavelets* wdm = d->wdm;
    fptr = fopen(fname,"w");
    if (!fptr) {
        fprintf(stderr, "Couldn't open %s for writing!\n", fname);
        return 1;
    }
    int k;
    for (int i=0; i<wdm->NT; i++) {
        double t = i*wdm->dt;
        for (int j=d->lmin; j < d->lmax; j++) {
            double f = j*wdm->df;
            wavelet_pixel_to_index(wdm, i, j, &k);
            k -= wdm->kmin;
            fprintf(fptr,"%lg ",t);
            fprintf(fptr,"%lg ",f);
            fprintf(fptr,"%lg",data[k]);
            fprintf(fptr,"\n");
        }
    }
    fclose(fptr);
    return 0;
}

int write_coarse_wdm_data(struct Data* d, double* data, char* fname, int Q) {
    FILE *fptr = NULL;
    struct Wavelets* wdm = d->wdm;
    fptr = fopen(fname,"w");
    if (!fptr) {
        fprintf(stderr, "Couldn't open %s for writing!\n", fname);
        return 1;
    }
    int Ncoarse = wdm->NT / Q;
    for (int q=0; q<Ncoarse; q++) {
        double t = q*wdm->dt*Q;
        for (int j=d->lmin; j < d->lmax; j++) {
            double f = j*wdm->df;
            int k = q + (j-d->lmin)*Ncoarse;
            fprintf(fptr,"%lg ",t);
            fprintf(fptr,"%lg ",f);
            fprintf(fptr,"%lg",data[k]);
            fprintf(fptr,"\n");
        }
    }
    fclose(fptr);
    return 0;
}

int main(int argc, char *argv[])
{
    fprintf(stdout, "\n============== NOISE WAVELET MCMC ===========\n");

    time_t start, stop;
    start = time(NULL);
    char filename[MAXSTRINGSIZE];
    
    /* check arguments */
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    if(argc==1) print_usage();

    /* Allocate data structures */
    struct Data *data   = malloc(sizeof(struct Data));
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    
    parse_data_args(argc,argv,data,orbit,flags,chain,"wavelet");
    printf("old fmin=%lg, fmax=%lg\n",data->fmin,data->fmax);
    if(flags->help) print_usage();

    // Time-varying instrument noise needs the dynamic (time-resolved) likelihood.
    if(flags->instDrift && flags->stationary) {
        fprintf(stderr,"--inst-drift (time-varying instrument noise) is incompatible with --stationary\n");
        exit(-1);
    }

    /* Setup output directories for data and chain files */
    sprintf(data->dataDir,"%s/data",flags->runDir);
    sprintf(chain->chainDir,"%s/chains",flags->runDir);
    sprintf(chain->chkptDir,"%s/checkpoint",flags->runDir);
    
    mkdir(flags->runDir,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data->dataDir,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


    /* Initialize data structures */
    data->Nchannel = 3;
    alloc_data(data, flags);

    /*
    Define and set up Orbit structure which contains spacecraft ephemerides
    */
    initialize_orbit(data, orbit, flags);
    initialize_interpolated_analytic_orbits(orbit, data->T, data->t0);

    printf("old fmin=%lg, fmax=%lg\n",data->fmin,data->fmax);
    // TODO these are not initialized because the UCB sampler set them itself later!
    // should probably have defaults... this does not feel like the source's job
    // we'll approximate command line choices
    data->lmin = (int)floor(data->fmin / data->wdm->df);
    data->lmax =  (int)ceil(data->fmax / data->wdm->df);
    //reset wavelet basis max and min ranges
    wavelet_pixel_to_index(data->wdm,0,data->lmin,&data->wdm->kmin); 
    wavelet_pixel_to_index(data->wdm,0,data->lmax,&data->wdm->kmax); 
    // this points kmin/kmax to first time pixel of first/last freq layer
    
    //recompute fmin and fmax so they align with a bin
    data->fmin = data->lmin*WAVELET_BANDWIDTH;
    data->fmax = data->lmax*WAVELET_BANDWIDTH;
    printf("new fmin=%lg, fmax=%lg\n",data->fmin,data->fmax);

    printf("total wavelet pixels: %d\n", data->N);
    if (data->wdm->NF%2==1 || data->wdm->NT%2==1) {
        fprintf(stderr, "Currently, cannot handle odd length NF or NT wdm basis!\n");
        exit(-1);
    }

    /* WDM time-axis coarse-graining factor Q (1 = full resolution). */
    int Q = flags->coarseQ;
    if (Q < 1) {
        fprintf(stderr, "--coarse-Q must be >= 1, got %d\n", Q);
        exit(-1);
    }
    if (data->wdm->NT % Q != 0) {
        fprintf(stderr, "--coarse-Q=%d must divide wdm->NT=%d\n", Q, data->wdm->NT);
        exit(-1);
    }
    if (flags->stationary && Q != data->wdm->NT) {
        printf("[stationary] overriding --coarse-Q=%d to NT=%d (stationary covariance does not vary in time; finer Q is redundant work)\n", Q, data->wdm->NT);
        Q = data->wdm->NT;
    }
    int Ncoarse = data->wdm->NT / Q;
    printf("coarse-graining: Q=%d, NT=%d -> Ncoarse=%d\n", Q, data->wdm->NT, Ncoarse);
    printf("active band: lmin=%d, lmax=%d, Nlayer=%d, total coarse pixels=%d\n",
            data->lmin, data->lmax, data->Nlayer, data->Nlayer*Ncoarse);

    /* Initialize chain structure and files */
    initialize_chain(chain, flags, &data->cseed, "w");

    /* read data */
    if(flags->strainData)
        ReadData(data,orbit,flags);
    else if(flags->simNoise) {
        // inject some noise
        struct InstrumentModel inst_inj = {0};
        initialize_instrument_model_wavelet(orbit, data, &inst_inj);
        if (flags->instDriftInjN > 0) {
            if (flags->stationaryConf) {
                fprintf(stderr,"[inst-drift-inj] requires dynamic injection; incompatible with --stationary-conf\n");
                exit(-1);
            }
            int Nlink = inst_inj.Nlink;
            if (flags->instDriftInjN != 1 && flags->instDriftInjN != 2*Nlink) {
                fprintf(stderr,"[inst-drift-inj] expects 1 value (uniform) or %d values (sacc[%d],soms[%d]), got %d\n",
                        2*Nlink, Nlink, Nlink, flags->instDriftInjN);
                exit(-1);
            }
            // linear drift injected per coefficient (or uniform when one value given).
            // 12-value layout matches the noise chain: sacc[0..5] then soms[0..5].
            inst_inj.drift_enabled = 1;
            for (int i=0; i<Nlink; i++) {
                double sacc_d = (flags->instDriftInjN==1) ? flags->instDriftInj[0] : flags->instDriftInj[i];
                double soms_d = (flags->instDriftInjN==1) ? flags->instDriftInj[0] : flags->instDriftInj[Nlink+i];
                if (fabs(sacc_d) >= 1.0 || fabs(soms_d) >= 1.0) {
                    fprintf(stderr,"[inst-drift-inj] |drift| >= 1 (sacc[%d]=%g, soms[%d]=%g) would drive the noise variance negative; use |v|<1\n",
                            i, sacc_d, i, soms_d);
                    exit(-1);
                }
                inst_inj.sacc_drift[i] = sacc_d;
                inst_inj.soms_drift[i] = soms_d;
            }
            // regenerate so the covariance "slope" reflects the injected drift
            generate_instrument_noise_model_wavelet(data->wdm, orbit, &inst_inj);
            printf("   ...injecting time-varying instrument noise (linear drift, %d value(s))\n", flags->instDriftInjN);
        }
        struct ForegroundModel conf_inj = {0};
        struct ForegroundModel *conf_ptr = NULL;
        if (flags->confNoise) {
            conf_ptr = &conf_inj;
            initialize_foreground_model_wavelet(orbit, data, conf_ptr);
        }
        struct SGWBModel sgwb_inj = {0};
        struct SGWBModel *sgwb_ptr = NULL;
        // NOTE: refactored so these pointers being NULL means excluded from cov matrix
        if (flags->sgwbTemplate>=0) {
            sgwb_ptr = &sgwb_inj;
            initialize_sgwb_model_wavelet(orbit, data, sgwb_ptr, flags->sgwbTemplate, flags->sgwbInjN ? flags->sgwbInjParams : NULL);
        }
        if (flags->stationaryConf)
            generate_full_stationary_covariance_matrix(data->wdm, &inst_inj, conf_ptr, sgwb_ptr, data->noise);
        else
            generate_full_dynamic_covariance_matrix(data->wdm, &inst_inj, conf_ptr, sgwb_ptr, data->noise);
        // alloc tdi data
        alloc_tdi(data->tdi, data->N, 3);

        //AddNoiseWavelet(data,data->tdi);
        MyAddNoiseWavelet(data,data->tdi);

        // mirror noise into dwt so the wavelet-domain state is consistent
        size_t band_bytes = (size_t)data->N * sizeof(double);
        memcpy(data->dwt->X, data->tdi->X, band_bytes);
        memcpy(data->dwt->Y, data->tdi->Y, band_bytes);
        memcpy(data->dwt->Z, data->tdi->Z, band_bytes);

        /*
        write_wdm_data(data, data->tdi->X, "./debug_wdm_noise.dat");
        write_wdm_data(data, data->noise->C[0][0], "./debug_wdm_C00.dat");
        */

    }

    /* Store DFT copy of simulated data */
    // TODO: fix below
    data->qmin = data->lmin;
    if(!flags->strainData) {
        // data->dft->X is sized for the active FFT band (2*data->NFFT doubles),
        // but wavelet_transform_inverse_freq writes a full ND+2. Likewise it
        // *reads* a full NF*NT scalogram, while data->dwt is only sized for the
        // active band. Route input and output through full-size scratch buffers
        // and copy the active slices in/out.
        int ND = data->wdm->NF * data->wdm->NT;
        double *full_fft  = malloc((ND + 2) * sizeof(double));
        double *full_wave = calloc(ND, sizeof(double));
        size_t band_bytes = (size_t)data->N * sizeof(double);
        size_t band_off   = (size_t)data->lmin * data->wdm->NT;
        size_t copy_bytes = 2 * (size_t)data->NFFT * sizeof(double);
        size_t fft_off    = 2 * (size_t)data->qmin;

        memcpy(full_wave + band_off, data->dwt->X, band_bytes);
        wavelet_transform_inverse_freq(data->wdm, full_wave, full_fft);
        memcpy(data->dft->X, full_fft + fft_off, copy_bytes);

        memcpy(full_wave + band_off, data->dwt->Y, band_bytes);
        wavelet_transform_inverse_freq(data->wdm, full_wave, full_fft);
        memcpy(data->dft->Y, full_fft + fft_off, copy_bytes);

        memcpy(full_wave + band_off, data->dwt->Z, band_bytes);
        wavelet_transform_inverse_freq(data->wdm, full_wave, full_fft);
        memcpy(data->dft->Z, full_fft + fft_off, copy_bytes);
        free(full_fft);
        free(full_wave);
    }
    
    /* print various data products for plotting */
    print_data(data, flags);

    /* Precompute sufficient statistics for the likelihood.
     * Q=1: Pij[k] = X[k]*Y[k] etc. the LL evaluates to the standard Whittle
     * Q>1: Welch/Bartlett-like average of Q fine pixels per coarse cell, Q corrects for Gamma dist
     * This is built once after data injection, shared read-only across threads.
     */
    struct CoarseStats stats;
    alloc_coarse_stats(&stats, data->Nlayer, Q, data->wdm->NT);
    coarse_grain_wavelet_data(data, &stats);
    if (Q > 1) {
        sprintf(filename,"%s/coarse_Pxx.dat", data->dataDir);
        write_coarse_wdm_data(data, (double*)stats.Pxx, filename, Q);
    }

    /* Welch-Satterthwaite effective dof, frozen at the injected covariance.
     * Reads data->noise (the full-resolution injected covariance), so it only
     * makes sense for injection runs. Without this call stats.Qeff stays at Q.
     */
    if (flags->ws_approx) {
        if (!flags->simNoise) {
            fprintf(stderr,"[ws-approx] requires an injection (--sim-noise); falling back to Qeff=Q\n");
        } else {
            precompute_coarse_Qeff(data, &stats);
            if (Q > 1) {
                sprintf(filename,"%s/coarse_Qeff.dat", data->dataDir);
                write_coarse_wdm_data(data, (double*)stats.Qeff, filename, Q);
            }
        }
    }

    /* Crop edge time pixels from the coarse likelihood. Quantized to whole coarse
     * cells: drop ceil(crop_time/Q) cells from each end so at least crop_time fine
     * pixels are removed. The P/Qeff arrays keep all Ncoarse cells; only the
     * likelihood sum (and its normalization) is restricted to [qlo,qhi). */
    if (flags->crop_time > 0) {
        int qcrop = (flags->crop_time + Q - 1) / Q;   // whole cells per end
        if (2*qcrop >= Ncoarse) {
            fprintf(stderr, "--crop-time=%d removes all coarse cells "
                    "(2*ceil(%d/%d)=%d >= Ncoarse=%d)\n",
                    flags->crop_time, flags->crop_time, Q, 2*qcrop, Ncoarse);
            exit(-1);
        }
        stats.qlo = qcrop;
        stats.qhi = Ncoarse - qcrop;
        printf("crop-time: dropping %d coarse cell(s) (= %d fine pixels) per end; "
               "likelihood uses cells [%d,%d) of %d\n",
               qcrop, qcrop*Q, stats.qlo, stats.qhi, Ncoarse);
    }


    // For now, we are not going to have full wavelet scaleograms stored everywhere!
    // we'll treat these as outer products of the spectrum and modulation
    // the modulation will be constant for everything except the confusion noise for now
    // TODO: eventually, allow for slow frequency-dependent modulations in the instrument model

    /* Initialize Instrument Noise Model */
    printf("   ...initialize instrument noise model\n");
    struct Noise **scaleogram = malloc(chain->NC*sizeof(struct Noise *)); // will save total scalograms here
    struct InstrumentModel **inst_model = malloc(chain->NC*sizeof(struct InstrumentModel *));
    struct InstrumentModel **inst_trial = malloc(chain->NC*sizeof(struct InstrumentModel *));
    for(int ic=0; ic<chain->NC; ic++)
    {
        scaleogram[ic] = malloc(sizeof(struct Noise));

        // see above note, we'll save the current scaleogram from all contributions here
        // but the individual models will only have spectrum \otimes modulation
        // Coarse layout: k = q + jrel*Ncoarse, Nlayer*Ncoarse total slots.
        // At Q=1, Ncoarse=NT and this reduces to the full-resolution indexing.
        alloc_noise(scaleogram[ic], data->Nlayer*Ncoarse, data->Nlayer, data->Nchannel);

        scaleogram[ic]->kmin = data->wdm->kmin;
        for (size_t j=0; j < (size_t)data->Nlayer; j++)
            for (size_t q=0; q < (size_t)Ncoarse; q++)
                scaleogram[ic]->f[q + j*Ncoarse] = (data->lmin + j)*data->wdm->df;

        // see above note. this is frequency-axis only
        inst_model[ic] = malloc(sizeof(struct InstrumentModel));
        inst_trial[ic] = malloc(sizeof(struct InstrumentModel));
        initialize_instrument_model_wavelet(orbit, data, inst_model[ic]);
        initialize_instrument_model_wavelet(orbit, data, inst_trial[ic]);

        // Time-varying instrument noise model: start the drift at zero (the
        // stationary state) and let the sampler explore it. Regenerate so the
        // covariance slope is consistent (zero) from the first likelihood call.
        if (flags->instDrift) {
            inst_model[ic]->drift_enabled = 1;
            inst_trial[ic]->drift_enabled = 1;
            generate_instrument_noise_model_wavelet(data->wdm, orbit, inst_model[ic]);
            generate_instrument_noise_model_wavelet(data->wdm, orbit, inst_trial[ic]);
        }
    }
    // noise model is only along freq axis for now
    sprintf(filename,"%s/instrument_noise_model.dat", data->dataDir);
    print_noise_model(inst_model[0]->psd, filename);

    /* Initialize Galactic Foreground Model */
    // note: UCB wavelet app hides this inside initialize_ucb_state if stationary noise, and in the waveform calls otherwise.
    // See GetDynamicNoiseModel and GetStationaryNoiseModel
    // remember, plan here for now is to not sample over foreground model in wavelet basis, but still sample over noise/sgwb
    if(flags->confNoise)printf("   ...initialize foreground noise model\n");
    struct ForegroundModel **conf_model = calloc(chain->NC,sizeof(struct ForegroundModel *));
    struct ForegroundModel **conf_trial = calloc(chain->NC,sizeof(struct ForegroundModel *));
    if (flags->confNoise) {
        for(int ic=0; ic<chain->NC; ic++)
        {
            conf_model[ic] = malloc(sizeof(struct ForegroundModel));
            conf_trial[ic] = malloc(sizeof(struct ForegroundModel));
            initialize_foreground_model_wavelet(orbit, data, conf_model[ic]);
            initialize_foreground_model_wavelet(orbit, data, conf_trial[ic]);
        }
        // NOTE: these do not include the modulation prefactors
        sprintf(filename,"%s/foreground_noise_model.dat",data->dataDir);
        print_noise_model(conf_model[0]->psd, filename);
    }

    /* Initialize SGWB Model */
    if(flags->sgwbTemplate>=0) printf("   ...initialize SGWB model\n");
    struct SGWBModel **sgwb_model = calloc(chain->NC,sizeof(struct SGWBModel *));
    struct SGWBModel **sgwb_trial = calloc(chain->NC,sizeof(struct SGWBModel *));
    if(flags->sgwbTemplate>=0) {
        for (int ic=0; ic<chain->NC; ic++)
        {
            sgwb_model[ic] = malloc(sizeof(struct SGWBModel));
            sgwb_trial[ic] = malloc(sizeof(struct SGWBModel));
            initialize_sgwb_model_wavelet(orbit, data, sgwb_model[ic], flags->sgwbTemplate, flags->sgwbInjN ? flags->sgwbInjParams : NULL);
            initialize_sgwb_model_wavelet(orbit, data, sgwb_trial[ic], flags->sgwbTemplate, flags->sgwbInjN ? flags->sgwbInjParams : NULL);
        }
        sprintf(filename,"%s/sgwb_noise_model.dat",data->dataDir);
        print_noise_model(sgwb_model[0]->psd, filename);
    }


    /* Combine noise components to form covariance matrix */
    for(int ic=0; ic<chain->NC; ic++)
    {
        if (flags->stationary)
            generate_full_stationary_covariance_matrix_coarse(data->wdm, Q, inst_model[ic], conf_model[ic], sgwb_model[ic], scaleogram[ic]);
        else
            generate_full_dynamic_covariance_matrix_coarse(data->wdm, Q, inst_model[ic], conf_model[ic], sgwb_model[ic], scaleogram[ic]);

        /* get initial likelihoods */
        invert_noise_covariance_matrix(scaleogram[ic]);
        double logL = my_noise_log_likelihood_wavelet_coarse(data, scaleogram[ic], &stats);
        inst_model[ic]->logL = logL;
        if (sgwb_model[ic])
            sgwb_model[ic]->logL = logL;
        if (conf_model[ic])
            conf_model[ic]->logL = logL;
    }

    sprintf(filename,"%s/full_noise_model.dat",data->dataDir);
    print_noise_model_dynamic_coarse(data, scaleogram[0], Q, filename);

    //MCMC
    printf("\n==== Noise Wavelet MCMC Sampler ====\n");
    
    sprintf(filename,"%s/noise_chain.dat",chain->chainDir);
    FILE *noiseChainFile = fopen(filename,"w");
    if (noiseChainFile == NULL) {fprintf(stderr, "Filesystem error, couldn't open %s for writing\n", filename); exit(-1);}

    FILE *foregroundChainFile = NULL;
    if(flags->confNoise)
    {
        sprintf(filename,"%s/foreground_chain.dat",chain->chainDir);
        foregroundChainFile = fopen(filename,"w");
        if (foregroundChainFile == NULL) {fprintf(stderr, "Filesystem error, couldn't open %s for writing\n", filename); exit(-1);}
    }
    FILE *sgwbChainFile = NULL;
    if(flags->sgwbTemplate>=0)
    {
        sprintf(filename,"%s/sgwb_chain.dat",chain->chainDir);
        sgwbChainFile = fopen(filename,"w");
        if (sgwbChainFile == NULL) {fprintf(stderr, "Filesystem error, couldn't open %s for writing\n", filename); exit(-1);}
    }
    
    int step = 0;
    int NC = chain->NC;

    // TODO improve paralellization
    #pragma omp parallel num_threads(flags->threads)
    {
        int threadID  = omp_get_thread_num();
        int numThreads = omp_get_num_threads();
        
        /* The MCMC loop */
        for(; step<flags->NMCMC;)
        {
            //#pragma omp barrier
            
            // (parallel) loop over chains
            for(int ic=threadID; ic<NC; ic+=numThreads)
            {
                struct Noise *psd_ptr = scaleogram[chain->index[ic]];
                struct InstrumentModel *inst_model_ptr = inst_model[chain->index[ic]];
                struct InstrumentModel *inst_trial_ptr = inst_trial[chain->index[ic]];
                struct ForegroundModel *conf_model_ptr = conf_model[chain->index[ic]];
                struct ForegroundModel *conf_trial_ptr = conf_trial[chain->index[ic]];
                struct SGWBModel       *sgwb_model_ptr = sgwb_model[chain->index[ic]];
                struct SGWBModel       *sgwb_trial_ptr = sgwb_trial[chain->index[ic]];

                for(int mc=0; mc<10; mc++)
                {
                    noise_instrument_model_mcmc_wavelet(orbit, data, inst_model_ptr, inst_trial_ptr, conf_model_ptr, sgwb_model_ptr, psd_ptr, &stats, chain, flags, ic);
                    // Note that we do not sample over the galaxy's modulation, only its spectral shape
                    if(flags->confNoise) noise_foreground_model_mcmc_wavelet(data, inst_model_ptr, conf_model_ptr, conf_trial_ptr, sgwb_model_ptr, psd_ptr, &stats, chain, flags, ic);
                    if(flags->sgwbTemplate>=0) noise_sgwb_model_mcmc_wavelet_dumb(data, inst_model_ptr, conf_model_ptr, sgwb_model_ptr, sgwb_trial_ptr, psd_ptr, &stats, chain, flags, ic);
                    // TODO: the logLs aren't tracked well at all here. We should probably refactor to have some kind of WaveletNoise struct that can handle all these components...
                }
            }// end loop over chains
            
            //Next section is single threaded. Every thread must get here before continuing
            
            #pragma omp barrier
            
            if(threadID==0)
            {
                noise_ptmcmc(inst_model, chain, flags);
                // note that this tracks the logL correctly, just swaps around chain->index
                
                if(flags->NMCMC >= 100 && step%(flags->NMCMC/100)==0)
                {
                    printf("noise_wavelet_mcmc at step %i\n",step);
                    generate_instrument_noise_model_wavelet(data->wdm, orbit, inst_model[chain->index[0]]);
                    sprintf(filename,"%s/current_instrument_noise_model.dat",data->dataDir);
                    print_noise_model(inst_model[chain->index[0]]->psd, filename);

                    if(flags->confNoise)
                    {
                        generate_galactic_foreground_model_wavelet(data->wdm, conf_model[chain->index[0]]);
                        sprintf(filename,"%s/current_foreground_noise_model.dat",data->dataDir);
                        print_noise_model(conf_model[chain->index[0]]->psd, filename);
                    }
                    if(flags->sgwbTemplate>=0) 
                    {
                        generate_sgwb_model_wavelet(data->wdm, sgwb_model[chain->index[0]]);
                        sprintf(filename,"%s/current_sgwb_noise_model.dat",data->dataDir);
                        print_noise_model(sgwb_model[chain->index[0]]->psd, filename);
                    }
                }
                
                // TODO: what should this code below do?
                /*
                if(step%data->downsample==0 && step/data->downsample < data->Nwave)
                {
                    generate_instrument_noise_model_wavelet(data->wdm, orbit, inst_model[chain->index[0]]);
                    if(flags->confNoise)
                    {
                        generate_galactic_foreground_model_wavelet(data->wdm, conf_model[chain->index[0]]);
                    }
                    if(flags->sgwbTemplate>=0) 
                    {
                        generate_sgwb_model_wavelet(data->wdm, sgwb_model[chain->index[0]]);
                    }
                    generate_full_dynamic_covariance_matrix(data->wdm, inst_model[chain->index[0]], conf_model[chain->index[0]], sgwb_model[chain->index[0]], scaleogram[chain->index[0]]);
                    
                    for(int n=0; n<data->N; n++)
                        for(int i=0; i<data->Nchannel; i++)
                            data->S_pow[n][i][step/data->downsample] = scaleogram[chain->index[0]]->C[i][i][n];
                }
                */

                fprintf(noiseChainFile,"%i %.12g ",step,inst_model[chain->index[0]]->logL);
                print_instrument_state(inst_model[chain->index[0]], noiseChainFile);
                fprintf(noiseChainFile,"\n");

                if(flags->confNoise)
                {
                    fprintf(foregroundChainFile,"%i %.12g ",step,conf_model[chain->index[0]]->logL);
                    print_foreground_state(conf_model[chain->index[0]], foregroundChainFile);
                    fprintf(foregroundChainFile,"\n");
                }
                if(flags->sgwbTemplate>=0)
                {
                    fprintf(sgwbChainFile,"%i %.12g ",step,sgwb_model[chain->index[0]]->logL);
                    print_sgwb_state(sgwb_model[chain->index[0]], sgwbChainFile);
                    fprintf(sgwbChainFile,"\n");
                }

                step++;


            }
            //Can't continue MCMC until single thread is finished
            #pragma omp barrier
            
        }// end of MCMC loop
        
    }// End of parallelization
    
    fclose(noiseChainFile);
    if(flags->confNoise)fclose(foregroundChainFile);
    if(flags->sgwbTemplate>=0)fclose(sgwbChainFile);

    generate_instrument_noise_model_wavelet(data->wdm, orbit, inst_model[chain->index[0]]);
    sprintf(filename,"%s/final_instrument_noise_model.dat", data->dataDir);
    print_noise_model(inst_model[chain->index[0]]->psd, filename);

    if(flags->confNoise)
    {
        generate_galactic_foreground_model_wavelet(data->wdm, conf_model[chain->index[0]]);
        sprintf(filename,"%s/final_foreground_noise_model.dat",data->dataDir);
        print_noise_model(conf_model[chain->index[0]]->psd, filename);
    }
    if(flags->sgwbTemplate>=0)
    {
        generate_sgwb_model_wavelet(data->wdm, sgwb_model[chain->index[0]]);
        sprintf(filename,"%s/final_sgwb_noise_model.dat",data->dataDir);
        print_noise_model(sgwb_model[chain->index[0]]->psd, filename);
    }
    if (flags->stationary)
        generate_full_stationary_covariance_matrix_coarse(data->wdm, Q, inst_model[chain->index[0]], conf_model[chain->index[0]], sgwb_model[chain->index[0]], scaleogram[chain->index[0]]);
    else
        generate_full_dynamic_covariance_matrix_coarse(data->wdm, Q, inst_model[chain->index[0]], conf_model[chain->index[0]], sgwb_model[chain->index[0]], scaleogram[chain->index[0]]);
    if(flags->confNoise || flags->sgwbTemplate>=0)
    {
        sprintf(filename,"%s/final_full_noise_model.dat",data->dataDir);
        print_noise_model_dynamic_coarse(data, scaleogram[chain->index[0]], Q, filename);
    }


    // TODO this hasn't been adapted for WDM
    //print_noise_reconstruction(data, flags);

    // print_whitened_data indexes scaleogram with the full-resolution
    // wavelet_pixel_to_index, so it is only valid at Q=1.
    if (Q == 1)
    {
        sprintf(filename,"%s/whitened_data.dat",data->dataDir);
        print_whitened_data(data, scaleogram[chain->index[0]], filename);
    }

    
    free_coarse_stats(&stats);

    //print total run time
    stop = time(NULL);

    printf(" ELAPSED TIME = %g seconds\n",(double)(stop-start));
    // TODO cleanup our mallocs?

    return 0;
}
