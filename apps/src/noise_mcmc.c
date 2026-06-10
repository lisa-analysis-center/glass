/*
 * Copyright 2024 Tyson B. Littenberg
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

/**
 @file noise_mcmc.c
 \brief Main function for stand-alone Noise parameterized model sampler 
 */

#include <glass_utils.h>
#include <glass_noise.h>

static void print_usage()
{
    print_glass_usage();
    fprintf(stdout,"EXAMPLE:\n");
    fprintf(stdout,"noise_mcmc --sim-noise --conf-noise --duration 7864320 --fmin 1e-4 --fmax 8e-3\n");
    fprintf(stdout,"\n");
    exit(0);
}

static void set_noise_defaults(struct Data *data)
{
    data->T        = 31457280; /* one "mldc years" at 15s sampling */
    data->t0       = 0.0; /* start time of data segment in seconds */
    data->sqT      = sqrt(data->T);
    data->NFFT     = (int)floor(data->T/LISA_CADENCE);   /* full data set */
    data->Nlayer   = (int)floor(1/2.0/LISA_CADENCE/WAVELET_BANDWIDTH); /* up to Nyquist */
    data->Nchannel = 3; //1=X, 2=AE, 3=XYZ
    data->qpad     = 0;
    data->fmin     = 1e-4; //Hz
    sprintf(data->basis,"fourier");
}

int main(int argc, char *argv[])
{
    fprintf(stdout, "\n================= NOISE MCMC ================\n");

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
    
    set_noise_defaults(data);
    parse_data_args(argc,argv,data,orbit,flags,chain,"fourier");
    if(flags->help) print_usage();
    
    /* Setup output directories for data and chain files */
    sprintf(data->dataDir,"%s/data",flags->runDir);
    sprintf(chain->chainDir,"%s/chains",flags->runDir);
    sprintf(chain->chkptDir,"%s/checkpoint",flags->runDir);
    
    mkdir(flags->runDir,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data->dataDir,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    /* Initialize data structures */
    alloc_data(data, flags);
    
    /* Initialize LISA orbit model */
    initialize_orbit(data, orbit, flags);
    
    /* Initialize chain structure and files */
    initialize_chain(chain, flags, &data->cseed, "a");
    
    /* read data */
    if(flags->strainData)
        ReadData(data,orbit,flags);
    else if (flags->simNoise)
        SimulateData(data, orbit, flags);
    
    /* Store DFT copy of simulated data */
    copy_tdi(data->tdi, data->dft);
    
    /* print various data products for plotting */
    print_data(data, flags);
    
    /* Initialize Instrument Noise Model */
    printf("   ...initialize instrument noise model\n");
    struct Noise **psd = malloc(chain->NC*sizeof(struct Noise *));
    struct InstrumentModel **inst_model = malloc(chain->NC*sizeof(struct InstrumentModel *));
    struct InstrumentModel **inst_trial = malloc(chain->NC*sizeof(struct InstrumentModel *));
    for(int ic=0; ic<chain->NC; ic++)
    {
        psd[ic] = malloc(sizeof(struct Noise));
        alloc_noise(psd[ic], data->NFFT, data->Nlayer, data->Nchannel);

        inst_model[ic] = malloc(sizeof(struct InstrumentModel));
        inst_trial[ic] = malloc(sizeof(struct InstrumentModel));
        initialize_instrument_model(orbit, data, inst_model[ic]);
        initialize_instrument_model(orbit, data, inst_trial[ic]);
    }
    sprintf(filename,"%s/instrument_noise_model.dat",data->dataDir);
    print_noise_model(inst_model[0]->psd, filename);


    /* Initialize Galactic Foreground Model */
    if(flags->confNoise)printf("   ...initialize foreground noise model\n");
    struct ForegroundModel **conf_model = malloc(chain->NC*sizeof(struct ForegroundModel *));
    struct ForegroundModel **conf_trial = malloc(chain->NC*sizeof(struct ForegroundModel *));
    for(int ic=0; ic<chain->NC; ic++)
    {
        conf_model[ic] = malloc(sizeof(struct ForegroundModel));
        conf_trial[ic] = malloc(sizeof(struct ForegroundModel));
        if(flags->confNoise) 
        {
           initialize_foreground_model(orbit, data, conf_model[ic]);
           initialize_foreground_model(orbit, data, conf_trial[ic]);
        }
    }
    if(flags->confNoise)
    {
        sprintf(filename,"%s/foreground_noise_model.dat",data->dataDir);
        print_noise_model(conf_model[0]->psd, filename);
    }

    /* Initialize SGWB Model */
    if(flags->sgwbTemplate>=0) printf("   ...initialize SGWB model\n");
    struct SGWBModel **sgwb_model = malloc(chain->NC*sizeof(struct SGWBModel *));
    struct SGWBModel **sgwb_trial = malloc(chain->NC*sizeof(struct SGWBModel *));
    for (int ic=0; ic<chain->NC; ic++)
    {
        sgwb_model[ic] = malloc(sizeof(struct SGWBModel));
        sgwb_trial[ic] = malloc(sizeof(struct SGWBModel));
        if(flags->sgwbTemplate>=0) 
        {
           initialize_sgwb_model(orbit, data, sgwb_model[ic], flags->sgwbTemplate, flags->sgwbInjN ? flags->sgwbInjParams : NULL);
           initialize_sgwb_model(orbit, data, sgwb_trial[ic], flags->sgwbTemplate, flags->sgwbInjN ? flags->sgwbInjParams : NULL);
        }
    }
    if(flags->sgwbTemplate>=0)
    {
        sprintf(filename,"%s/sgwb_noise_model.dat",data->dataDir);
        print_noise_model(sgwb_model[0]->psd, filename);
    }

    /* TODO: SGWB injections??? */
    

    /* Combine noise components to form covariance matrix */
    for(int ic=0; ic<chain->NC; ic++)
    {
        if(flags->confNoise) generate_full_covariance_matrix(inst_model[ic]->psd, conf_model[ic]->psd, data->Nchannel);
        if(flags->sgwbTemplate>=0) generate_full_covariance_matrix(inst_model[ic]->psd, sgwb_model[ic]->psd, data->Nchannel);
    }

    /* get initial likelihood */
    for(int ic=0; ic<chain->NC; ic++)
    {
        invert_noise_covariance_matrix(inst_model[ic]->psd);
        inst_model[ic]->logL = my_noise_log_likelihood(data, inst_model[ic]->psd);
        conf_model[ic]->logL = inst_model[ic]->logL;
        sgwb_model[ic]->logL = inst_model[ic]->logL;
    }

    sprintf(filename,"%s/full_noise_model.dat",data->dataDir);
    print_noise_model(inst_model[0]->psd, filename);

    
    
    //MCMC
    printf("\n==== Noise MCMC Sampler ====\n");
    
    sprintf(filename,"%s/noise_chain.dat",chain->chainDir);
    FILE *noiseChainFile = fopen(filename,"w");

    FILE *foregroundChainFile = NULL;
    if(flags->confNoise)
    {
        sprintf(filename,"%s/foreground_chain.dat",chain->chainDir);
        foregroundChainFile = fopen(filename,"w");
    }
    FILE *sgwbChainFile = NULL;
    if(flags->sgwbTemplate>=0)
    {
        sprintf(filename,"%s/sgwb_chain.dat",chain->chainDir);
        sgwbChainFile = fopen(filename,"w");
    }
    
    int mcmc = 0;
    int NC = chain->NC;

    /* The MCMC loop */
    for(; mcmc < flags->NMCMC;)
    {
        
        // (parallel) loop over chains
        #pragma omp parallel for num_threads(flags->threads)
        for(int ic=0; ic<NC; ic++)
        {
            struct Noise *psd_ptr = psd[chain->index[ic]];
            struct InstrumentModel *inst_model_ptr = inst_model[chain->index[ic]];
            struct InstrumentModel *inst_trial_ptr = inst_trial[chain->index[ic]];
            struct ForegroundModel *conf_model_ptr = conf_model[chain->index[ic]];
            struct ForegroundModel *conf_trial_ptr = conf_trial[chain->index[ic]];
            struct SGWBModel       *sgwb_model_ptr = sgwb_model[chain->index[ic]];
            struct SGWBModel       *sgwb_trial_ptr = sgwb_trial[chain->index[ic]];

            for(int steps=0; steps<10; steps++)
            {
                noise_instrument_model_mcmc(orbit, data, inst_model_ptr, inst_trial_ptr, conf_model_ptr, sgwb_model_ptr, psd_ptr, chain, flags, ic);
                if(flags->confNoise) noise_foreground_model_mcmc(data, inst_model_ptr, conf_model_ptr, conf_trial_ptr, sgwb_model_ptr, psd_ptr, chain, flags, ic);
                if(flags->sgwbTemplate>=0) noise_sgwb_model_mcmc(data, inst_model_ptr, conf_model_ptr, sgwb_model_ptr, sgwb_trial_ptr, psd_ptr, chain, flags, ic);
            }
        }// end (parallel) loop over chains
            
        //Next section is single threaded. Every thread must get here before continuing
            
        noise_ptmcmc(inst_model, chain, flags);
            
        if(mcmc%(flags->NMCMC/10)==0)printf("noise_mcmc at step %i\n",mcmc);
            
        // print chain files
        fprintf(noiseChainFile,"%i %.12g ",mcmc,inst_model[chain->index[0]]->logL);
        print_instrument_state(inst_model[chain->index[0]], noiseChainFile);
        fprintf(noiseChainFile,"\n");

        if(flags->confNoise) 
        {
            fprintf(foregroundChainFile,"%i %.12g ",mcmc,conf_model[chain->index[0]]->logL);
            print_foreground_state(conf_model[chain->index[0]], foregroundChainFile);
            fprintf(foregroundChainFile,"\n");
        }

        if(flags->sgwbTemplate>=0) 
        {
            fprintf(sgwbChainFile,"%i %.12g ",mcmc,sgwb_model[chain->index[0]]->logL);
            print_sgwb_state(sgwb_model[chain->index[0]], sgwbChainFile);
            fprintf(sgwbChainFile,"\n");
        }

        if(mcmc%(flags->NMCMC/10)==0)
        {
            generate_instrument_noise_model(orbit,inst_model[chain->index[0]]);
            sprintf(filename,"%s/current_instrument_noise_model.dat",data->dataDir);
            print_noise_model(inst_model[chain->index[0]]->psd, filename);

            if(flags->confNoise)
            {
                generate_galactic_foreground_model(conf_model[chain->index[0]]);
                sprintf(filename,"%s/current_foreground_noise_model.dat",data->dataDir);
                print_noise_model(conf_model[chain->index[0]]->psd, filename);
            }
            if(flags->sgwbTemplate>=0) 
            {
                generate_sgwb_model(sgwb_model[chain->index[0]]);
                sprintf(filename,"%s/current_sgwb_noise_model.dat",data->dataDir);
                print_noise_model(sgwb_model[chain->index[0]]->psd, filename);
            }
        }
            
        if(mcmc%data->downsample==0 && mcmc/data->downsample < data->Nwave)
        {
            generate_instrument_noise_model(orbit,inst_model[chain->index[0]]);
            if(flags->confNoise)
            {
                generate_galactic_foreground_model(conf_model[chain->index[0]]);
                generate_full_covariance_matrix(inst_model[chain->index[0]]->psd, conf_model[chain->index[0]]->psd, data->Nchannel);
            }
            if(flags->sgwbTemplate>=0) 
            {
                generate_sgwb_model(sgwb_model[chain->index[0]]);
                generate_full_covariance_matrix(inst_model[chain->index[0]]->psd,sgwb_model[chain->index[0]]->psd, data->Nchannel);
            }
            
            for(int n=0; n<data->NFFT; n++)
                for(int i=0; i<data->Nchannel; i++)
                    data->S_pow[n][i][mcmc/data->downsample] = inst_model[chain->index[0]]->psd->C[i][i][n];
        }

        mcmc++;
            
            
    } // end of mcmc loop
    
    
    fclose(noiseChainFile);
    if(flags->confNoise)fclose(foregroundChainFile);
    if(flags->sgwbTemplate>=0)fclose(sgwbChainFile);

    generate_instrument_noise_model(orbit,inst_model[chain->index[0]]);
    sprintf(filename,"%s/final_instrument_noise_model.dat",data->dataDir);
    print_noise_model(inst_model[chain->index[0]]->psd, filename);

    if(flags->confNoise)
    {
        generate_galactic_foreground_model(conf_model[chain->index[0]]);
        sprintf(filename,"%s/final_foreground_noise_model.dat",data->dataDir);
        print_noise_model(conf_model[chain->index[0]]->psd, filename);
        generate_full_covariance_matrix(inst_model[chain->index[0]]->psd, conf_model[chain->index[0]]->psd, data->Nchannel);
    }
    if(flags->sgwbTemplate>=0)
    {
        generate_sgwb_model(sgwb_model[chain->index[0]]);
        sprintf(filename,"%s/final_sgwb_noise_model.dat",data->dataDir);
        print_noise_model(sgwb_model[chain->index[0]]->psd, filename);
        generate_full_covariance_matrix(inst_model[chain->index[0]]->psd, sgwb_model[chain->index[0]]->psd, data->Nchannel);
    }
    if(flags->confNoise || flags->sgwbTemplate>=0)
    {
        sprintf(filename,"%s/final_full_noise_model.dat",data->dataDir);
        print_noise_model(inst_model[chain->index[0]]->psd, filename);
    }
    
    print_noise_reconstruction(data, flags);

    sprintf(filename,"%s/whitened_data.dat",data->dataDir);
    // this was here, but it was already added??
    //if(flags->confNoise)generate_full_covariance_matrix(inst_model[chain->index[0]]->psd,conf_model[chain->index[0]]->psd, data->Nchannel);
    print_whitened_data(data, inst_model[chain->index[0]]->psd, filename);

    
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g seconds\n",(double)(stop-start));
    
    
    return 0;
}

