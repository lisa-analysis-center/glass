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

#include <assert.h>
#include "glass_utils.h"
#include "glass_noise.h"
#include "gitversion.h"

#define FIXME 0

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


void printProgress (double percentage)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    fprintf(stdout, "\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}


void print_version(FILE *fptr)
{
    fprintf(fptr, "\n");
    fprintf(fptr, "=============== GLASS Version: ==============\n\n");
    //fprintf(fptr, "  Git remote origin: %s\n", GIT_URL);
    //fprintf(fptr, "  Git version: %s\n", GIT_VER);
    fprintf(fptr, "  Git commit: %s\n", GITVERSION);
    //fprintf(fptr, "  Git commit author: %s\n",GIT_AUTHOR);
    //fprintf(fptr, "  Git commit date: %s\n", GIT_DATE);
    fprintf(fptr, "\n=============================================\n\n");
}

void setup_run_directories(struct Flags *flags, struct Data *data, struct Chain *chain)
{
    
    sprintf(data->dataDir,"%s/data",flags->runDir);
    sprintf(chain->chainDir,"%s/chains",flags->runDir);
    sprintf(chain->chkptDir,"%s/checkpoint",flags->runDir);

    mkdir(flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data->dataDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

}

void initialize_orbit(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    /* Load spacecraft ephemerides */
    switch(flags->orbit)
    {
        case 0:
            initialize_analytic_orbit(orbit);
            break;
        case 1:
            initialize_numeric_orbit(orbit);
            break;
        default:
            fprintf(stderr,"unsupported orbit type\n");
            exit(1);
            break;
    }
    
    /* set approximate f/fstar for segment */
    data->sine_f_on_fstar = sin((data->fmin + (data->fmax-data->fmin)/2.)/orbit->fstar);
}

//TODO: Move file pointers to Model instead of Chain structs
void initialize_chain(struct Chain *chain, struct Flags *flags, unsigned int *seed, const char *mode)
{
    int ic;
    int NC = chain->NC;
    char filename[MAXSTRINGSIZE];

    chain->index = calloc(NC,sizeof(int));
    chain->acceptance = calloc(NC,sizeof(double));
    chain->temperature = calloc(NC,sizeof(double));
    chain->avgLogL     = calloc(NC,sizeof(double));
    chain->dimension   = calloc(NC,sizeof(int *));
    for(ic=0; ic<NC; ic++)
    {
        chain->index[ic]=ic;
        chain->acceptance[ic] = 1.0;
        chain->temperature[ic] = pow(1.2,(double)ic);
        chain->avgLogL[ic] = 0.0;
        chain->dimension[ic] = calloc(flags->DMAX,sizeof(int));
        for(int id=0; id<flags->DMAX; id++) chain->dimension[ic][id] = 0;
    }
    //set hottest chain to ~infinite temperature
    if(NC>1) chain->temperature[NC-1] = 1e12;
    chain->logLmax = 0.0;
    
    chain->r = malloc(NC*sizeof(unsigned int *));
    
    for(ic=0; ic<NC; ic++)
    {
        //set seed
        chain->r[ic] = *seed;
        
        //evolve seed
        rand_r_U_0_1(seed);
    }
    
    if(!flags->quiet)
    {
        sprintf(filename,"%s/log_likelihood_chain.dat",chain->chainDir);
        chain->likelihoodFile = fopen(filename,mode);
        
        sprintf(filename,"%s/temperature_chain.dat",chain->chainDir);
        chain->temperatureFile = fopen(filename,mode);
    }
    
    chain->chainFile = malloc(NC*sizeof(FILE *));
    sprintf(filename,"%s/model_chain.dat.0",chain->chainDir);
    chain->chainFile[0] = fopen(filename,mode);
    
    chain->parameterFile = malloc(NC*sizeof(FILE *));
    sprintf(filename,"%s/parameter_chain.dat.0",chain->chainDir);
    chain->parameterFile[0] = fopen(filename,mode);
    
    chain->dimensionFile = malloc(flags->DMAX*sizeof(FILE *));
    for(int i=0; i<flags->DMAX; i++)
    {
        /* only create these files when needed */
        chain->dimensionFile[i]=NULL;
    }
    
    chain->noiseFile = malloc(NC*sizeof(FILE *));
    sprintf(filename,"%s/noise_chain.dat.0",chain->chainDir);
    chain->noiseFile[0] = fopen(filename,mode);
    
    if(flags->confNoise)
    {
        chain->foregroundFile = malloc(NC*sizeof(FILE *));
        sprintf(filename,"%s/foreground_chain.dat.0",chain->chainDir);
        chain->foregroundFile[0] = fopen(filename,mode);
    }
    
    if(flags->calibration)
    {
        chain->calibrationFile = malloc(NC*sizeof(FILE *));
        sprintf(filename,"%s/calibration_chain.dat.0",chain->chainDir);
        chain->calibrationFile[0] = fopen(filename,mode);
    }
    
    if(flags->verbose)
    {
        for(ic=1; ic<NC; ic++)
        {
            sprintf(filename,"%s/parameter_chain.dat.%i",chain->chainDir,ic);
            chain->parameterFile[ic] = fopen(filename,mode);
            
            sprintf(filename,"%s/model_chain.dat.%i",chain->chainDir,ic);
            chain->chainFile[ic] = fopen(filename,mode);
            
            sprintf(filename,"%s/noise_chain.dat.%i",chain->chainDir,ic);
            chain->noiseFile[ic] = fopen(filename,mode);
        }
    }
}

void alloc_data(struct Data *data, struct Flags *flags)
{
    int NMCMC = flags->NMCMC;

    data->logN = log((double)(data->N*data->Nchannel));
    
    data->tdi   = malloc(sizeof(struct TDI));
    data->raw   = malloc(sizeof(struct TDI));
    data->dft   = malloc(sizeof(struct TDI));
    data->dwt   = malloc(sizeof(struct TDI));
    data->noise = malloc(sizeof(struct Noise));
            
    alloc_tdi(data->tdi, data->N, data->Nchannel);
    alloc_tdi(data->raw, data->N, data->Nchannel);
    alloc_tdi(data->dft, data->N, data->Nchannel);
    alloc_tdi(data->dwt, data->N, data->Nchannel);
    if(!strcmp(data->basis,"fourier")) alloc_noise(data->noise, data->NFFT, data->Nlayer, data->Nchannel);
    if(!strcmp(data->basis,"wavelet")) alloc_noise(data->noise, data->N, data->Nlayer, data->Nchannel);
    
    //reconstructed signal model
    int i_re,i_im;
    data->h_rec = malloc(data->N*sizeof(double **));
    data->h_res = malloc(data->N*sizeof(double **));
    if(!strcmp(data->basis,"fourier"))
    {
        data->r_pow = malloc(data->NFFT*sizeof(double **));
        data->h_pow = malloc(data->NFFT*sizeof(double **));
        data->S_pow = malloc(data->NFFT*sizeof(double **));
    }
    if(!strcmp(data->basis,"wavelet"))
    {
        data->r_pow = malloc(data->N*sizeof(double **));
        data->h_pow = malloc(data->N*sizeof(double **));
        data->S_pow = malloc(data->N*sizeof(double **));
    }    
    
    //number of waveform samples to save
    data->Nwave=100;
    
    //downsampling rate of post-burn-in samples
    data->downsample = NMCMC/data->Nwave;
    if (!data->downsample) data->downsample=1;
    
    if(!strcmp(data->basis,"fourier"))
    {
        for(int i=0; i<data->NFFT; i++)
        {
            i_re = i*2;
            i_im = i_re+1;
            
            data->S_pow[i]    = malloc(data->Nchannel*sizeof(double *));
            data->h_pow[i]    = malloc(data->Nchannel*sizeof(double *));
            data->r_pow[i]    = malloc(data->Nchannel*sizeof(double *));
            data->h_rec[i_re] = malloc(data->Nchannel*sizeof(double *));
            data->h_rec[i_im] = malloc(data->Nchannel*sizeof(double *));
            data->h_res[i_re] = malloc(data->Nchannel*sizeof(double *));
            data->h_res[i_im] = malloc(data->Nchannel*sizeof(double *));
            for(int n=0; n<data->Nchannel; n++)
            {
                data->S_pow[i][n]    = calloc(data->Nwave,sizeof(double));
                data->h_pow[i][n]    = calloc(data->Nwave,sizeof(double));
                data->r_pow[i][n]    = calloc(data->Nwave,sizeof(double));
                data->h_rec[i_re][n] = calloc(data->Nwave,sizeof(double));
                data->h_rec[i_im][n] = calloc(data->Nwave,sizeof(double));
                data->h_res[i_re][n] = calloc(data->Nwave,sizeof(double));
                data->h_res[i_im][n] = calloc(data->Nwave,sizeof(double));
            }
        }
    }

    if(!strcmp(data->basis,"wavelet"))
    {
        for(int i=0; i<data->N; i++)
        {
            
            data->S_pow[i] = malloc(data->Nchannel*sizeof(double *));
            data->h_pow[i] = malloc(data->Nchannel*sizeof(double *));
            data->r_pow[i] = malloc(data->Nchannel*sizeof(double *));
            data->h_rec[i] = malloc(data->Nchannel*sizeof(double *));
            data->h_res[i] = malloc(data->Nchannel*sizeof(double *));
            for(int n=0; n<data->Nchannel; n++)
            {
                data->S_pow[i][n] = calloc(data->Nwave,sizeof(double));
                data->h_pow[i][n] = calloc(data->Nwave,sizeof(double));
                data->r_pow[i][n] = calloc(data->Nwave,sizeof(double));
                data->h_rec[i][n] = calloc(data->Nwave,sizeof(double));
                data->h_res[i][n] = calloc(data->Nwave,sizeof(double));
            }
        }
    }
    
    //Spectrum proposal
    data->p = calloc(data->N,sizeof(double));

    // Setup wavelet basis
    if(!strcmp(data->basis,"wavelet"))
    {
        data->wdm = malloc(sizeof(struct Wavelets));
        initialize_wavelet(data->wdm, data->T);
    }
}

void alloc_noise(struct Noise *noise, int N, int Nlayer, int Nchannel)
{
    noise->N = N;
    noise->Nlayer = Nlayer;
    noise->Nchannel = Nchannel;
    
    noise->eta = calloc(Nchannel*Nlayer,sizeof(double));

    noise->f = calloc(N,sizeof(double));

    noise->C    = malloc(Nchannel*sizeof(double **));
    noise->invC = malloc(Nchannel*sizeof(double **));
    
    for(int i=0; i<Nchannel; i++)
    {
        for(int j=0; j<Nlayer; j++) noise->eta[i*Nlayer+j] = 1.0;
        noise->C[i]    = malloc(Nchannel*sizeof(double *));
        noise->invC[i] = malloc(Nchannel*sizeof(double *));
        
        for(int j=0; j<Nchannel; j++)
        {
            noise->C[i][j]    = calloc(N,sizeof(double));
            noise->invC[i][j] = calloc(N,sizeof(double));
        }
    }

    noise->logdetC  = malloc(N*sizeof(double));
    for (size_t i=0; i<N; i++) noise->logdetC[i] = -INFINITY;
    noise->transfer = calloc(N,sizeof(double));
    
    int n;
    for(n=0; n<N; n++)
    {
        for(int i=0; i<Nchannel; i++) noise->C[i][i][n] = 1.0;
        for(int i=0; i<Nchannel; i++)
        {
            for(int j=i+1; j<Nchannel; j++)
            {
                noise->C[i][j][n] = 0.0;
                noise->C[j][i][n] = 0.0;
            }
        }
        noise->transfer[n] = 1.0;
    }
}

void mpi_send_noise(struct Noise *noise, int dest)
{
    //pack up and send integers
    int int_values[4];
    int_values[0] = noise->N;
    int_values[1] = noise->kmin;
    int_values[2] = noise->Nchannel;
    int_values[3] = noise->Nlayer;
    MPI_Send(&int_values, 4, MPI_INT, dest, 0, MPI_COMM_WORLD);
    
    //pack up and send doubles
    int N = noise->N;
    int Nlayer = noise->Nlayer;
    int Nchannel = noise->Nchannel;
    int Nsend = Nchannel*Nlayer + 3*N;
    double *double_values = double_vector(Nsend);
    
    for(int n=0; n<N; n++)
    {
        double_values[n]     = noise->logdetC[n];
        double_values[n+N]   = noise->transfer[n];
        double_values[n+2*N] = noise->f[n];
    }
    for(int n=0; n<Nchannel*Nlayer; n++) double_values[n+3*N] = noise->eta[n];
    MPI_Send(double_values, Nsend, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
    
    //pack up and send covariance matrices
    int Ncij = Nchannel*Nchannel*N;
    Nsend = 2*Ncij;
    double *cij_values = double_vector(Nsend);
    int m;
    for(int i=0; i<Nchannel; i++)
    {
        for(int j=0; j<Nchannel; j++)
        {
            for(int n=0; n<N; n++)
            {
                m = i*Nchannel + j*Nchannel + n;
                cij_values[m]      = noise->C[i][j][n];
                cij_values[m+Ncij] = noise->invC[i][j][n];
            }
        }
    }
    MPI_Send(cij_values, Nsend, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);

    free_double_vector(cij_values);
    free_double_vector(double_values);
}

void mpi_receive_noise(struct Noise *noise, int source)
{
    MPI_Status status;
    
    //receive and unpack integers
    int int_values[4];
    MPI_Recv(&int_values, 4, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
    noise->N        = int_values[0];
    noise->kmin     = int_values[1];
    noise->Nchannel = int_values[2];
    noise->Nlayer   = int_values[3];
    
    //receive and unpack doubles
    int N = noise->N;
    int Nlayer = noise->Nlayer;
    int Nchannel = noise->Nchannel;
    int Nrecv = Nchannel*Nlayer + 3*N;
    double *double_values = double_vector(Nrecv);
    MPI_Recv(double_values, Nrecv, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
    for(int n=0; n<N; n++)
    {
        noise->logdetC[n]     = double_values[n];
        noise->transfer[n] = double_values[n+N];
        noise->f[n]        = double_values[n+2*N];
    }
    for(int n=0; n<Nchannel*Nlayer; n++) noise->eta[n] = double_values[n+3*N];

    //receive and unpack covariance matrices
    int Ncij = Nchannel*Nchannel*N;
    Nrecv = 2*Ncij;
    double *cij_values = double_vector(Nrecv);
    MPI_Recv(cij_values, Nrecv, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
    int m;
    for(int i=0; i<Nchannel; i++)
    {
        for(int j=0; j<Nchannel; j++)
        {
            for(int n=0; n<N; n++)
            {
                m = i*Nchannel + j*Nchannel + n;
                noise->C[i][j][n]    = cij_values[m];
                noise->invC[i][j][n] = cij_values[m+Ncij];
            }
        }
    }

    free_double_vector(double_values);
    free_double_vector(cij_values);
}

void alloc_calibration(struct Calibration *calibration)
{
    calibration->dampA = 0.0;
    calibration->dampE = 0.0;
    calibration->dampX = 0.0;
    calibration->dphiA = 0.0;
    calibration->dphiE = 0.0;
    calibration->dphiX = 0.0;
    calibration->real_dphiA = 1.0;
    calibration->real_dphiE = 1.0;
    calibration->real_dphiX = 1.0;
    calibration->imag_dphiA = 0.0;
    calibration->imag_dphiE = 0.0;
    calibration->imag_dphiX = 0.0;
}

//TODO: Expand copy_data() to include everything, and then replace where needed (NoiseWrapper.c, ...)
void copy_data(struct Data *origin, struct Data *copy)
{
    memcpy(copy->format, origin->format, sizeof(origin->format));
    memcpy(copy->basis, origin->basis, sizeof(origin->basis));
    memcpy(copy->fileName, origin->fileName, sizeof(origin->fileName));
    copy->T=origin->T;
    copy->sqT=origin->sqT;
    copy->N=origin->N;
    copy->NFFT=origin->NFFT;
    copy->Nlayer=origin->Nlayer;
    copy->Nchannel=origin->Nchannel;
    copy->qpad=origin->qpad;
    copy->cseed=origin->cseed;
    copy->nseed=origin->nseed;
    copy->iseed=origin->iseed;
    copy->t0=origin->t0;
    //TODO: need copy_wavelet
}

void mpi_send_data(struct Data *data, int dest)
{
    //pack up and send integers
    int int_values[9];
    int_values[0] = data->N;
    int_values[1] = data->NFFT;
    int_values[2] = data->Nchannel;
    int_values[3] = data->Nlayer;
    int_values[4] = data->qmin;
    int_values[5] = data->qmax;
    int_values[6] = data->qpad;
    int_values[7] = data->lmin;
    int_values[8] = data->lmax;
    MPI_Send(&int_values, 9, MPI_INT, dest, 0, MPI_COMM_WORLD);
    
    //pack up send doubles
    double double_values[9];
    double_values[0] = data->logN;
    double_values[1] = data->T;
    double_values[2] = data->sqT;
    double_values[3] = data->t0;
    double_values[4] = data->fmin;
    double_values[5] = data->fmax;
    double_values[6] = data->sine_f_on_fstar;
    double_values[7] = data->pmax;
    double_values[8] = data->SNR2;
    MPI_Send(&double_values, 9, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
    
    //send TDI structures
    mpi_send_tdi(data->tdi, dest);
    mpi_send_tdi(data->raw, dest);
    mpi_send_tdi(data->dft, dest);
    mpi_send_tdi(data->dwt, dest);
}

void mpi_receive_data(struct Data *data, int source)
{
    MPI_Status status;
    
    //receive and unpack integers
    int int_values[9];
    MPI_Recv(&int_values, 9, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
    data->N        = int_values[0];
    data->NFFT     = int_values[1];
    data->Nchannel = int_values[2];
    data->Nlayer   = int_values[3];
    data->qmin     = int_values[4];
    data->qmax     = int_values[5];
    data->qpad     = int_values[6];
    data->lmin     = int_values[7];
    data->lmax     = int_values[8];
    
    //receive and unpack doubles
    double double_values[9];
    MPI_Recv(&double_values, 9, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
    data->logN            = double_values[0];
    data->T               = double_values[1];
    data->sqT             = double_values[2];
    data->t0              = double_values[3];
    data->fmin            = double_values[4];
    data->fmax            = double_values[5];
    data->sine_f_on_fstar = double_values[6];
    data->pmax            = double_values[7];
    data->SNR2            = double_values[8];
    
    //send TDI structures
    mpi_receive_tdi(data->tdi, source);
    mpi_receive_tdi(data->raw, source);
    mpi_receive_tdi(data->dft, source);
    mpi_receive_tdi(data->dwt, source);
}


void copy_noise(struct Noise *origin, struct Noise *copy)
{
    copy->N = origin->N;
    copy->Nlayer = origin->Nlayer;
    copy->Nchannel = origin->Nchannel;
    memcpy(copy->eta,origin->eta,origin->Nchannel*origin->Nlayer*sizeof(double));

    memcpy(copy->f, origin->f, origin->N*sizeof(double));

    copy_Cij(origin->C, copy->C, origin->Nchannel, origin->N);
    copy_Cij(origin->invC, copy->invC, origin->Nchannel, origin->N);

    memcpy(copy->logdetC, origin->logdetC, origin->N*sizeof(double));
    memcpy(copy->transfer, origin->transfer, origin->N*sizeof(double));
}

void copy_Cij(double ***origin, double ***copy, int M, int N)
{
    for(int i=0; i<M; i++)
        for(int j=0; j<M; j++)
            memcpy(copy[i][j], origin[i][j], N*sizeof(double));
}

void copy_calibration(struct Calibration *origin, struct Calibration *copy)
{
    copy=origin;
    /*
    copy->dampA   = origin->dampA;
    copy->dampE   = origin->dampE;
    copy->dampX   = origin->dampX;
    copy->dphiA = origin->dphiA;
    copy->dphiE = origin->dphiE;
    copy->dphiX = origin->dphiX;
    copy->real_dphiA = origin->real_dphiA;
    copy->real_dphiE = origin->real_dphiE;
    copy->real_dphiX = origin->real_dphiX;
    copy->imag_dphiA = origin->imag_dphiA;
    copy->imag_dphiE = origin->imag_dphiE;
    copy->imag_dphiX = origin->imag_dphiX;
     */
}

void free_noise(struct Noise *noise)
{
    free(noise->eta);
    free(noise->f);
    for(int i=0; i<noise->Nchannel; i++)
    {
        for(int j=0; j<noise->Nchannel; j++)
        {
            free(noise->C[i][j]);
            free(noise->invC[i][j]);
        }
        free(noise->C[i]);
        free(noise->invC[i]);
    }
    free(noise->C);
    free(noise->invC);
    free(noise->logdetC);
    free(noise->transfer);
    free(noise);
}

void free_chain(struct Chain *chain, struct Flags *flags)
{
    free(chain->index);
    free(chain->acceptance);
    free(chain->temperature);
    free(chain->avgLogL);
    free(chain->dimension);
    free(chain->r);
    
    if(!flags->quiet)
    {
        fclose(chain->likelihoodFile);
        fclose(chain->temperatureFile);
    }

    fclose(chain->chainFile[0]);

    fclose(chain->parameterFile[0]);
    
    for(int i=0; i<flags->DMAX; i++)
    {
        /* only create these files when needed */
        if(chain->dimensionFile[i]!=NULL) fclose(chain->dimensionFile[i]);
    }
    free(chain->dimensionFile);

    fclose(chain->noiseFile[0]);
    
    if(flags->calibration)
    {
        fclose(chain->calibrationFile[0]);
        free(chain->calibrationFile);
    }

    
    if(flags->verbose)
    {
        for(int ic=1; ic<chain->NC; ic++)
        {
            fclose(chain->chainFile[ic]);
            fclose(chain->parameterFile[ic]);
            fclose(chain->noiseFile[ic]);
        }
    }
    free(chain->chainFile);
    free(chain->parameterFile);
    free(chain->noiseFile);
    
    free(chain);
}

void free_calibration(struct Calibration *calibration)
{
    free(calibration);
}

void ReadHDF5(struct Data *data, struct TDI *tdi, struct TDI *tdi_dwt, struct Flags *flags)
{
    /* LDASOFT-formatted structure for TDI data */
    struct TDI *tdi_td = malloc(sizeof(struct TDI));
        
    if(!strcmp(data->format,"frequency"))  LISA_Read_HDF5_LDC_RADLER_TDI(tdi_td, data->fileName);
    if(!strcmp(data->format,"sangria")) LISA_Read_HDF5_LDC_TDI(tdi_td, data->fileName, "/obs/tdi");
    
    
    /* Select time segment of full data set */
    double start_time = data->t0;
    double stop_time = start_time + data->T;
    double dt = tdi_td->delta;
    double Tobs = stop_time - start_time;
    int N = (int)floor(Tobs/dt);

    if(data->N>N)
    {
        printf("  Error: Requested amount of data is larger than %s\n",data->fileName);
        printf("         data->N = %i\n",data->N);
        printf("         N       = %i\n",N);
        printf("         data->Tobs = %.12g\n", data->T);
        printf("         Tobs       = %.12g\n", Tobs);
        printf("         Reduce requested --duration\n");
        exit(1);
    }
    
    /* work space for selecting and transforming time series */
    double *X = malloc(N*sizeof(double));
    double *Y = malloc(N*sizeof(double));
    double *Z = malloc(N*sizeof(double));
    
    double *Xtime = malloc(N*sizeof(double));
    double *Ytime = malloc(N*sizeof(double));
    double *Ztime = malloc(N*sizeof(double));

    /* Allocate data->tdi structure for Fourier transform output */
    alloc_tdi(tdi, N, N_TDI_CHANNELS);
    alloc_tdi(tdi_dwt, N, N_TDI_CHANNELS);
    tdi->delta = 1./Tobs;

    /* Select requested time segment */
    int n_start = (int)floor(start_time/dt); // first sample of time segment
    
    for(int n=0; n<N; n++)
    {
        int m = n_start+n;
        Xtime[n] = tdi_td->X[m];
        Ytime[n] = tdi_td->Y[m];
        Ztime[n] = tdi_td->Z[m];
    }
    
    /* lets get rid of those black holes */
    if(flags->no_mbh)
    {
        struct TDI *tdi_td_mbhb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_mbhb, data->fileName, "/sky/mbhb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            Xtime[n] -= tdi_td_mbhb->X[m];
            Ytime[n] -= tdi_td_mbhb->Y[m];
            Ztime[n] -= tdi_td_mbhb->Z[m];
        }
        free_tdi(tdi_td_mbhb);
    }
    
    /* lets get rid of the galaxy */
    if(flags->no_ucb)
    {
        struct TDI *tdi_td_dgb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_dgb, data->fileName, "/sky/dgb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            Xtime[n] -= tdi_td_dgb->X[m];
            Ytime[n] -= tdi_td_dgb->Y[m];
            Ztime[n] -= tdi_td_dgb->Z[m];
        }
        free_tdi(tdi_td_dgb);
    
        struct TDI *tdi_td_igb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_igb, data->fileName, "/sky/igb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            Xtime[n] -= tdi_td_igb->X[m];
            Ytime[n] -= tdi_td_igb->Y[m];
            Ztime[n] -= tdi_td_igb->Z[m];
        }
        free_tdi(tdi_td_igb);
    }

    /* lets get rid of the verification binaries */
    if(flags->no_vgb)
    {
        struct TDI *tdi_td_vgb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_vgb, data->fileName, "/sky/vgb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            Xtime[n] -= tdi_td_vgb->X[m];
            Ytime[n] -= tdi_td_vgb->Y[m];
            Ztime[n] -= tdi_td_vgb->Z[m];
        }
        free_tdi(tdi_td_vgb);
    }

    /* Detrend data */
    detrend(Xtime, N, (int)(FILTER_LENGTH/LISA_CADENCE));
    detrend(Ytime, N, (int)(FILTER_LENGTH/LISA_CADENCE));
    detrend(Ztime, N, (int)(FILTER_LENGTH/LISA_CADENCE));
    
    /* Tukey window time-domain TDI channels tdi_td */
    double alpha = (2.0*FILTER_LENGTH/Tobs);
    
    tukey(Xtime, alpha, N);
    tukey(Ytime, alpha, N);
    tukey(Ztime, alpha, N);
    
    /* Fourier transform time-domain TDI channels */
    for(int n=0; n<N; n++)
    {
        X[n] = Xtime[n];
        Y[n] = Ytime[n];
        Z[n] = Ztime[n];
    }

    glass_forward_real_fft(X,N);
    glass_forward_real_fft(Y,N);
    glass_forward_real_fft(Z,N);
    
    /* Normalize FD data */
    double rft_norm = 1.0;//sqrt(Tobs)/(double)N;
    
    /* Account for losses from windowing
    double tukey_norm = tukey_scale(alpha, N);
    rft_norm /= tukey_norm;
    */
    
    for(int n=0; n<N; n++)
    {
        tdi->X[n] = X[n]*rft_norm;
        tdi->Y[n] = Y[n]*rft_norm;
        tdi->Z[n] = Z[n]*rft_norm;
    }
    
    /* lets get rid of the high frequency ucbs */
    if(flags->no_ucb_hi)
    {
        double *Xgal = calloc(N,sizeof(double));
        double *Ygal = calloc(N,sizeof(double));
        double *Zgal = calloc(N,sizeof(double));

        struct TDI *tdi_td_dgb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_dgb, data->fileName, "/sky/dgb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            Xgal[n] += tdi_td_dgb->X[m];
            Ygal[n] += tdi_td_dgb->Y[m];
            Zgal[n] += tdi_td_dgb->Z[m];
        }
        free_tdi(tdi_td_dgb);
        
        struct TDI *tdi_td_igb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_igb, data->fileName, "/sky/igb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            Xgal[n] += tdi_td_igb->X[m];
            Ygal[n] += tdi_td_igb->Y[m];
            Zgal[n] += tdi_td_igb->Z[m];
        }
        free_tdi(tdi_td_igb);
        
        tukey(Xgal, alpha, N);
        tukey(Ygal, alpha, N);
        tukey(Zgal, alpha, N);

        glass_forward_real_fft(Xgal, N);
        glass_forward_real_fft(Ygal, N);
        glass_forward_real_fft(Zgal, N);

        
        /* Allocate data->tdi structure for Fourier transform output */
        struct TDI *tdi_gal = malloc(sizeof(struct TDI));
        alloc_tdi(tdi_gal, N, N_TDI_CHANNELS);
        tdi_gal->delta = 1./Tobs;

        for(int n=0; n<N; n++)
        {
            tdi_gal->X[n] = Xgal[n] * rft_norm;
            tdi_gal->Y[n] = Ygal[n] * rft_norm;
            tdi_gal->Z[n] = Zgal[n] * rft_norm;
        }

        /* remove hi-f binaries */
        for(int n=0; n<N/2; n++)
        {
            double f = (double)n/Tobs;
            if(f>0.00504)
            {
                tdi->X[2*n]   -= tdi_gal->X[2*n];
                tdi->X[2*n+1] -= tdi_gal->X[2*n+1];

                tdi->Y[2*n]   -= tdi_gal->Y[2*n];
                tdi->Y[2*n+1] -= tdi_gal->Y[2*n+1];

                tdi->Z[2*n]   -= tdi_gal->Z[2*n];
                tdi->Z[2*n+1] -= tdi_gal->Z[2*n+1];

            }
        }
        
        free_tdi(tdi_gal);
        free(Xgal);
        free(Ygal);
        free(Zgal);
    }

    /* populate AET channels because I can't let go */
    for(int n=0; n<N; n++) XYZ2AET(tdi->X[n], tdi->Y[n], tdi->Z[n], &tdi->A[n], &tdi->E[n], &tdi->T[n]);


    /* Wavelet transform time-domain TDI channels */
    if(!strcmp(data->basis,"wavelet"))
    {
        for(int n=0; n<N; n++)
        {
            X[n] = Xtime[n];
            Y[n] = Ytime[n];
            Z[n] = Ztime[n];
        }

        wavelet_transform_timefreq(data->wdm, X);
        wavelet_transform_timefreq(data->wdm, Y);
        wavelet_transform_timefreq(data->wdm, Z);
        
        for(int n=0; n<N; n++)
        {
            tdi_dwt->X[n] = X[n];
            tdi_dwt->Y[n] = Y[n];
            tdi_dwt->Z[n] = Z[n];
        }
        
        /* populate AET channels because I can't let go */
        for(int n=0; n<N; n++) XYZ2AET(tdi_dwt->X[n], tdi_dwt->Y[n], tdi_dwt->Z[n], &tdi_dwt->A[n], &tdi_dwt->E[n], &tdi_dwt->T[n]);
        
    }
    
    /* Free memory */
    free_tdi(tdi_td);
    free(X);
    free(Y);
    free(Z);
    free(Xtime);
    free(Ytime);
    free(Ztime);

}

void ReadASCII(struct Data *data, struct TDI *tdi)
{
    // for now, assume the data is saved as ascii fft, then convert to wavelet if needed

    double f;
    double junk;

    FILE *fptr = fopen(data->fileName,"r");
    
    if (data->Nchannel == 2) {
        //count number of samples
        int Nsamples = 0;
        while(!feof(fptr))
        {
            int check = fscanf(fptr,"%lg %lg %lg %lg %lg",&f,&junk,&junk,&junk,&junk);
            if(!check)
            {
                fprintf(stderr,"Error reading %s\n",data->fileName);
                exit(1);
            }
            Nsamples++;
        }
        rewind(fptr);
        Nsamples--;
        
        //load full dataset into TDI structure
        alloc_tdi(tdi, 2*Nsamples, 3);
        
        for(int n=0; n<Nsamples; n++)
        {
            int check = fscanf(fptr,"%lg %lg %lg %lg %lg",&f,&tdi->A[2*n],&tdi->A[2*n+1],&tdi->E[2*n],&tdi->E[2*n+1]);
            if(!check)
            {
                fprintf(stderr,"Error reading %s\n",data->fileName);
                exit(1);
            }
            
        }
    }
    else {
        //count number of samples
        int Nsamples = 0;
        while(!feof(fptr))
        {
            int check = fscanf(fptr,"%lg %lg %lg %lg %lg %lg %lg",&f,&junk,&junk,&junk,&junk, &junk, &junk);
            if(!check)
            {
                fprintf(stderr,"Error reading %s\n",data->fileName);
                exit(1);
            }
            Nsamples++;
        }
        rewind(fptr);
        Nsamples--;
        
        //load full dataset into TDI structure
        alloc_tdi(tdi, 2*Nsamples, 3);
        
        for(int n=0; n<Nsamples; n++)
        {
            int check = fscanf(fptr,"%lg %lg %lg %lg %lg %lg %lg",&f,&tdi->X[2*n],&tdi->X[2*n+1],&tdi->Y[2*n],&tdi->Y[2*n+1],&tdi->Z[2*n],&tdi->Z[2*n+1]);
            if(!check)
            {
                fprintf(stderr,"Error reading %s\n",data->fileName);
                exit(1);
            }
            
        }
    }
    fclose(fptr);

}

void ReadASCII_timeseries(struct Data* data, struct TDI* tdi_dft, struct TDI* tdi_dwt)
{
    double t;
    double junk;

    FILE *fptr = fopen(data->fileName,"r");
    // only can read XYZ
    data->Nchannel = 3;
    
    //count number of samples
    int Nsamples = 0;
    while(!feof(fptr))
    {
        int check = fscanf(fptr,"%lg %lg %lg %lg",&t,&junk,&junk,&junk);
        if(!check)
        {
            fprintf(stderr,"Error reading %s\n",data->fileName);
            exit(1);
        }
        Nsamples++;
    }
    rewind(fptr);
    Nsamples--;

    struct TDI tdi_td = {0};
    alloc_tdi(&tdi_td, Nsamples, 3);
    for(int n=0; n<Nsamples; n++)
    {
        int check = fscanf(fptr,"%lg %lg %lg %lg", &t, &tdi_td.X[n], &tdi_td.Y[n], &tdi_td.Z[n]);
        if(!check)
        {
            fprintf(stderr,"Error reading %s\n",data->fileName);
            exit(1);
        }
            
    }
    double Tobs = t;


    //load full dataset into TDI structure
    alloc_tdi(tdi_dft, Nsamples, 3);
    tdi_dft->delta = 1./Tobs;
    alloc_tdi(tdi_dwt, Nsamples, 3);

    // TODO: in principle, detrend, window here

    glass_forward_real_fft_outplace(tdi_td.X, tdi_dft->X, Nsamples);
    glass_forward_real_fft_outplace(tdi_td.Y, tdi_dft->Y, Nsamples);
    glass_forward_real_fft_outplace(tdi_td.Z, tdi_dft->Z, Nsamples);


    /* Wavelet transform time-domain TDI channels */
    if(!strcmp(data->basis,"wavelet"))
    {

        wavelet_transform_freq(data->wdm, tdi_dft->X, tdi_dwt->X); 
        wavelet_transform_freq(data->wdm, tdi_dft->Y, tdi_dwt->Y); 
        wavelet_transform_freq(data->wdm, tdi_dft->Z, tdi_dwt->Z); 
        
        /* populate AET channels because I can't let go */
        for(int n=0; n<Nsamples; n++) XYZ2AET(tdi_dwt->X[n], tdi_dwt->Y[n], tdi_dwt->Z[n], &tdi_dwt->A[n], &tdi_dwt->E[n], &tdi_dwt->T[n]);
        
    }
}

void ReadData(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    if(!flags->quiet)
        
    fprintf(stdout,"\n================== ReadData =================\n");

    /* load full dataset */
    struct TDI *tdi_full_dft = malloc(sizeof(struct TDI));
    struct TDI *tdi_full_dwt = malloc(sizeof(struct TDI));
    
    if(flags->hdf5Data)
        ReadHDF5(data,tdi_full_dft,tdi_full_dwt,flags);
    else {
        if (flags->timeseries) {
            ReadASCII_timeseries(data,tdi_full_dft, tdi_full_dwt);
        } else {
            ReadASCII(data,tdi_full_dft);
        }
    }

    
    
    /* select frequency segment */
    data->fmax = data->fmin + data->NFFT/data->T;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->NFFT;
    
    //store frequency segment in TDI structure
    for(int n=0; n<data->N; n++)
    {
        int m = data->qmin*2+n;
        data->dft->X[n] = tdi_full_dft->X[m];
        data->dft->Y[n] = tdi_full_dft->Y[m];
        data->dft->Z[n] = tdi_full_dft->Z[m];
        data->dft->A[n] = tdi_full_dft->A[m];
        data->dft->E[n] = tdi_full_dft->E[m];
        data->dft->T[n] = tdi_full_dft->T[m];
    }
    
    /* select wavelet layers */
    if(!strcmp(data->basis,"wavelet"))
    {
        //set lmin and lmax to represent frequency layer instead of bin
        data->lmin = (int)floor(data->fmin/WAVELET_BANDWIDTH);
        data->lmax = data->lmin + data->Nlayer;

        printf("  Minimum frequency layer=%i, maximum layer=%i\n",data->lmin,data->lmax-1);
        printf("  fmin=%lg, fmax=%lg\n",data->fmin,data->fmax);

        //reset wavelet basis max and min ranges
        wavelet_pixel_to_index(data->wdm,0,data->lmin,&data->wdm->kmin);
        wavelet_pixel_to_index(data->wdm,0,data->lmax,&data->wdm->kmax);

        //store frequency segment in TDI structure
        for(int n=0; n<data->N; n++)
        {
            int m = data->wdm->kmin+n;
            data->dwt->X[n] = tdi_full_dwt->X[m];
            data->dwt->Y[n] = tdi_full_dwt->Y[m];
            data->dwt->Z[n] = tdi_full_dwt->Z[m];
            data->dwt->A[n] = tdi_full_dwt->A[m];
            data->dwt->E[n] = tdi_full_dwt->E[m];
            data->dwt->T[n] = tdi_full_dwt->T[m];
        }
    }
    
    /* copy correct representation of data into the main tdi structure */
    for(int n=0; n<data->N; n++)
    {
        if(!strcmp(data->basis,"fourier"))
        {
            data->tdi->X[n] = data->dft->X[n];
            data->tdi->Y[n] = data->dft->Y[n];
            data->tdi->Z[n] = data->dft->Z[n];
            data->tdi->A[n] = data->dft->A[n];
            data->tdi->E[n] = data->dft->E[n];
            data->tdi->T[n] = data->dft->T[n];
        }
        else if(!strcmp(data->basis,"wavelet"))
        {
            data->tdi->X[n] = data->dwt->X[n];
            data->tdi->Y[n] = data->dwt->Y[n];
            data->tdi->Z[n] = data->dwt->Z[n];
            data->tdi->A[n] = data->dwt->A[n];
            data->tdi->E[n] = data->dwt->E[n];
            data->tdi->T[n] = data->dwt->T[n];
        }
    }
     
    //free memory
    free_tdi(tdi_full_dft);
    free_tdi(tdi_full_dwt);
    
    fprintf(stdout,"=============================================\n");
}

void GetNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    double Spm, Sop;
    
    //if you are simulating/fitting the noise
    if(!flags->psd)
    {
        for(int n=0; n<data->NFFT; n++)
        {
            double f = data->fmin + (double)(n)/data->T;
            data->noise->f[n] = f;
            data->noise->transfer[n] = noise_transfer_function(f/orbit->fstar);

            if(strcmp(data->format,"phase")==0)
            {
                data->noise->C[0][0][n] = AEnoise(orbit->L, orbit->fstar, f);
                data->noise->C[1][1][n] = AEnoise(orbit->L, orbit->fstar, f);
                data->noise->C[0][1][n] = 0.0;
                if(flags->confNoise)
                {
                    data->noise->C[0][0][n] += GBnoise(data->T,f);
                    data->noise->C[1][1][n] += GBnoise(data->T,f);
                }
            }
            else if(strcmp(data->format,"frequency")==0)
            {
                get_noise_levels("radler",f,&Spm,&Sop);
                data->noise->C[0][0][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[1][1][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[0][1][n] = 0.0;
                if(flags->confNoise)
                {
                    data->noise->C[0][0][n] += GBnoise_FF(data->T, orbit->fstar, f);
                    data->noise->C[1][1][n] += GBnoise_FF(data->T, orbit->fstar, f);
                }
            }
            else if(strcmp(data->format,"sangria")==0)
            {
                //TODO: Need a sqrt(2) to match Sangria data/noise
                get_noise_levels("sangria",f,&Spm,&Sop);
                data->noise->C[0][0][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[1][1][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[0][1][n] = 0.0;
                if(flags->confNoise)
                {
                    data->noise->C[0][0][n] += GBnoise_FF(data->T, orbit->fstar, f);
                    data->noise->C[1][1][n] += GBnoise_FF(data->T, orbit->fstar, f);
                }
            }
            else
            {
                fprintf(stderr,"Unsupported data format %s\n",data->format);
                exit(1);
            }

            //TODO: 3-channel model only has support for Sangria data conventions
            if(data->Nchannel==3)
            {
                //TODO: Need a sqrt(2) to match Sangria data/noise
                get_noise_levels("sangria",f,&Spm,&Sop);
                data->noise->C[0][0][n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[1][1][n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[2][2][n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);

                data->noise->C[0][1][n] = data->noise->C[1][0][n] = XYZcross_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[0][2][n] = data->noise->C[2][0][n] = XYZcross_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[1][2][n] = data->noise->C[2][1][n] = XYZcross_FF(orbit->L, orbit->fstar, f, Spm, Sop);

                if(flags->confNoise)
                {
                    double GBnoise=GBnoise_FF(data->T, orbit->fstar, f)/1.5; //GBnoise_FF() is hard-coded for AE channels
                    data->noise->C[0][0][n] += GBnoise;
                    data->noise->C[1][1][n] += GBnoise;
                    data->noise->C[2][2][n] += GBnoise;
                    data->noise->C[0][1][n] += -0.5*GBnoise;
                    data->noise->C[0][2][n] += -0.5*GBnoise;
                    data->noise->C[1][2][n] += -0.5*GBnoise;
                    data->noise->C[1][0][n] += -0.5*GBnoise;
                    data->noise->C[2][0][n] += -0.5*GBnoise;
                    data->noise->C[2][1][n] += -0.5*GBnoise;
                }
                
                /*normalize*/
                data->noise->C[0][0][n] /= 4.;
                data->noise->C[0][1][n] /= 4.;
                data->noise->C[0][2][n] /= 4.;
                data->noise->C[1][0][n] /= 4.;
                data->noise->C[1][1][n] /= 4.;
                data->noise->C[1][2][n] /= 4.;
                data->noise->C[2][0][n] /= 4.;
                data->noise->C[2][1][n] /= 4.;
                data->noise->C[2][2][n] /= 4.;
            }
        }
        
        invert_noise_covariance_matrix(data->noise);

    }
    //use PSD from file
    else
    {
        
        //parse input PSD file
        FILE *psdFile = fopen(flags->psdFile,"r");
        int lines=0;
        double f_temp, SnA_temp, SnE_temp;
        while(!feof(psdFile))
        {
            fscanf(psdFile,"%lg %lg %lg",&f_temp,&SnA_temp,&SnE_temp);
            lines++;
        }
        rewind(psdFile);
        lines--;
        
        double *f   = malloc(lines*sizeof(double));
        double *SnA = malloc(lines*sizeof(double));
        double *SnE = malloc(lines*sizeof(double));
        
        for(int l=0; l<lines; l++) fscanf(psdFile,"%lg %lg %lg",&f[l],&SnA[l],&SnE[l]);
        
        //interpolate input psd onto segment grid
        double *fint = malloc(data->NFFT*sizeof(double));
        for(int n=0; n<data->NFFT; n++) fint[n] = data->fmin + (double)(n)/data->T;

        CubicSplineGLASS(lines, f, SnA, data->NFFT, fint, data->noise->C[0][0]);
        CubicSplineGLASS(lines, f, SnE, data->NFFT, fint, data->noise->C[1][1]);
        
        free(f);
        free(SnA);
        free(SnE);
        fclose(psdFile);
    }
}

void MyAddNoise(struct Data *data, struct TDI *tdi)
{
    printf("   ...adding (Robbie's version of a) Gaussian noise realization\n");
    //set RNG for noise
    unsigned int r = data->nseed;

    int Nc = data->Nchannel;
    double L[Nc][Nc]; // will be Cholesky decomposition of covariance
    double C[Nc][Nc]; // non-ragged copy of covariance
    double n[Nc]; // vector of random normals
    // force 0 DC component (will probably be cropped out later anyway)
    tdi->X[0] = 0.0; tdi->X[1] = 0.0;
    tdi->Y[0] = 0.0; tdi->Y[1] = 0.0;
    tdi->Z[0] = 0.0; tdi->Z[1] = 0.0;
    for (int k=1; k<data->NFFT; k++) {
        for (int i=0; i<Nc; i++)
            for (int j=0; j<Nc; j++) {
                C[i][j] = data->noise->C[i][j][k];
            }
        my_cholesky_decomp(Nc, C, L);
        // real part
        for (int m=0; m<Nc; m++)
            n[m] = rand_r_N_0_1(&r) / M_SQRT2; // note 2 here -- variance of real or imaginary part
        for (int m=0; m<Nc; m++) {
            tdi->X[2*k] += L[0][m] * n[m];
            tdi->Y[2*k] += L[1][m] * n[m];
            tdi->Z[2*k] += L[2][m] * n[m];
        }
        // imag part
        for (int m=0; m<Nc; m++)
            n[m] = rand_r_N_0_1(&r) / M_SQRT2; // note 2 here -- variance of real or imaginary part
        for (int m=0; m<Nc; m++) {
            tdi->X[2*k+1] += L[0][m] * n[m];
            tdi->Y[2*k+1] += L[1][m] * n[m];
            tdi->Z[2*k+1] += L[2][m] * n[m];
        }
    }
}

void AddNoise(struct Data *data, struct TDI *tdi)
{
    
    printf("   ...adding Gaussian noise realization\n");
    
    //set RNG for noise
    unsigned int r = data->nseed;
    
    double n_re[data->Nchannel];
    double n_im[data->Nchannel];
    double u_re[data->Nchannel];
    double u_im[data->Nchannel];
    
    //get LU decomposition of covariance matrix
    double **L = malloc(data->Nchannel*sizeof(double*));
    double **C = malloc(data->Nchannel*sizeof(double*));
    for(int i=0; i<data->Nchannel; i++)
    {
        L[i] = malloc(data->Nchannel*sizeof(double));
        C[i] = malloc(data->Nchannel*sizeof(double));
    }
    
    
    
    for(int n=0; n<data->NFFT; n++)
    {
        for(int i=0; i<data->Nchannel; i++)
        {
            u_re[i] = rand_r_N_0_1(&r);
            u_im[i] = rand_r_N_0_1(&r);
            n_re[i] = n_im[i] = 0.0;
        }
 
        // make sure both diagonals of the covariance matrix are filled
        for(int i=0; i<data->Nchannel; i++)
            for(int j=i; j<data->Nchannel; j++)
                C[i][j] = C[j][i] = data->noise->C[i][j][n];

        cholesky_decomp(C, L, data->Nchannel);

        // n = Lu
        for(int i=0; i<data->Nchannel; i++)
        {
            for(int j=0; j<data->Nchannel; j++)
            {
                n_re[i] += L[i][j]*u_re[j]/M_SQRT2;
                n_im[i] += L[i][j]*u_im[j]/M_SQRT2;
            }
        }
        
        switch(data->Nchannel)
        {
            case 1:
                tdi->X[2*n]   += n_re[0];
                tdi->X[2*n+1] += n_im[0];
                break;
            case 2:
                tdi->A[2*n]   += n_re[0];
                tdi->A[2*n+1] += n_im[0];
                tdi->E[2*n]   += n_re[1];
                tdi->E[2*n+1] += n_im[1];
                break;
            case 3:
                tdi->X[2*n]   += n_re[0];
                tdi->X[2*n+1] += n_im[0];
                tdi->Y[2*n]   += n_re[1];
                tdi->Y[2*n+1] += n_im[1];
                tdi->Z[2*n]   += n_re[2];
                tdi->Z[2*n+1] += n_im[2];
                break;
        }
    }

    for(int i=0; i<data->Nchannel; i++)
    {
        free(L[i]);
        free(C[i]);
    }
    free(L);
    free(C);
}
void MyAddNoiseWavelet(struct Data *data, struct TDI *tdi)
{
    printf("   ...adding (Robbie's version of a) Gaussian noise realization\n");
    double L[3][3] = {0}; // will be Cholesky decomposition of covariance
    double C[3][3] = {0}; // non-ragged copy of covariance
    double n[3] = {0}; // vector of random normals
    unsigned int r = data->nseed;
    for (int i=0; i<data->wdm->NT; i++) {
        for (int j=data->lmin; j<data->lmax; j++) {
            int k;
            wavelet_pixel_to_index(data->wdm,i,j,&k);
            k -= data->wdm->kmin;
            for (int m=0; m<3; m++)
                for (int n=0; n<3; n++) {
                    C[m][n] = data->noise->C[m][n][k];
                }
            my_cholesky_decomp(3, C, L);
            for (int m=0; m<3; m++)
                n[m] = rand_r_N_0_1(&r); // real wavelet coefficient gets full variance, no sqrt(2) factor
            for (int m=0; m<3; m++) {
                tdi->X[k] += L[0][m] * n[m];
                tdi->Y[k] += L[1][m] * n[m];
                tdi->Z[k] += L[2][m] * n[m];
            }
        }
    }
}

void AddNoiseWavelet(struct Data *data, struct TDI *tdi)
{
    
    printf("   ...adding Gaussian noise realization\n");
    
    //set RNG for noise
    unsigned int r = data->nseed;
    
    double n[data->Nchannel];
    double u[data->Nchannel];
    
    //get LU decomposition of covariance matrix
    double **L = malloc(data->Nchannel*sizeof(double*));
    double **C = malloc(data->Nchannel*sizeof(double*));
    for(int i=0; i<data->Nchannel; i++)
    {
        L[i] = malloc(data->Nchannel*sizeof(double));
        C[i] = malloc(data->Nchannel*sizeof(double));
    }
    
    
    int k;
    struct Wavelets *wdm = data->wdm;
    for(int i=0; i<wdm->NT; i++)
    {
        for(int j=data->lmin; j<data->lmax; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);
            k-=data->wdm->kmin;

            for(int a=0; a<data->Nchannel; a++)
            {
                u[a] = rand_r_N_0_1(&r);
                n[a] = 0.0;
            }
    
            // make sure both diagonals of the covariance matrix are filled
            for(int a=0; a<data->Nchannel; a++)
                for(int b=a; b<data->Nchannel; b++)
                    C[a][b] = C[b][a] = data->noise->C[a][b][k];


            cholesky_decomp(C, L, data->Nchannel);

            // n = Lu
            for(int a=0; a<data->Nchannel; a++)
            {
                for(int b=0; b<data->Nchannel; b++)
                {
                    n[a] += L[a][b]*u[b];
                }
            }
            
            switch(data->Nchannel)
            {
                case 1:
                    tdi->X[k] += n[0];
                    break;
                case 2:
                    tdi->A[k] += n[0];
                    tdi->E[k] += n[1];
                    break;
                case 3:
                    tdi->X[k] += n[0];
                    tdi->Y[k] += n[1];
                    tdi->Z[k] += n[2];
                    break;
            }
        }
    }

    for(int i=0; i<data->Nchannel; i++)
    {
        free(L[i]);
        free(C[i]);
    }
    free(L);
    free(C);
}

void SimulateData(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    if(!flags->quiet) fprintf(stdout,"\n==== SimulateData ====\n");
    struct TDI *tdi = data->tdi;

    //get max and min samples
    data->fmax = data->fmin + data->NFFT/data->T;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->NFFT;

    //Get noise spectrum for data segment
    // TODO: switching here to use inst model
    //GetNoiseModel(data,orbit,flags);
    struct InstrumentModel *inst = malloc(sizeof(struct InstrumentModel));
    initialize_instrument_model(orbit, data, inst);
    copy_noise(inst->psd, data->noise);
    free_instrument_model(inst); 

    //Add Gaussian noise to injection
    if(flags->simNoise) MyAddNoise(data,tdi);

    //print various data products for plotting
    print_data(data, flags);
}

void print_data(struct Data *data, struct Flags *flags)
{
    int k;
    FILE *fptr;
    char filename[256];
    struct TDI *tdi = NULL;
    int lendir = 0;
    
    
    /* Power spectra */
    tdi = data->dft;
    lendir = strnlen(data->dataDir, MAXSTRINGSIZE);
    if (lendir==0) return;
    sprintf(filename,"%s/power_data.dat",data->dataDir);
    fptr=fopen(filename,"w");
    for(int i=0; i<data->NFFT; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        switch(data->Nchannel)
        {
            case 1:
                fprintf(fptr,"%.12g %lg\n", f, tdi->X[2*i]*tdi->X[2*i]+tdi->X[2*i+1]*tdi->X[2*i+1]);
                break;
            case 2:
                fprintf(fptr,"%.12g %lg %lg\n", f, tdi->A[2*i]*tdi->A[2*i]+tdi->A[2*i+1]*tdi->A[2*i+1], tdi->E[2*i]*tdi->E[2*i]+tdi->E[2*i+1]*tdi->E[2*i+1]);
                break;
            case 3:
                fprintf(fptr,"%.12g %lg %lg %lg\n", f, tdi->X[2*i]*tdi->X[2*i]+tdi->X[2*i+1]*tdi->X[2*i+1], tdi->Y[2*i]*tdi->Y[2*i]+tdi->Y[2*i+1]*tdi->Y[2*i+1], tdi->Z[2*i]*tdi->Z[2*i]+tdi->Z[2*i+1]*tdi->Z[2*i+1]);
                break;
        }
    }
    fclose(fptr);

    /* Scaleogram */
    if(!strcmp("wavelet",data->basis))
    {
        tdi = data->dwt;
        sprintf(filename,"%s/scaleogram_data.dat",data->dataDir);
        fptr=fopen(filename,"w");

        for(int i=0; i<data->wdm->NT; i++)
        {
            for(int j=data->lmin; j<data->lmax; j++)
            {
                wavelet_pixel_to_index(data->wdm,i,j,&k);
                k-=data->wdm->kmin;
                fprintf(fptr,"%lg %lg %.14e %.14e %.14e\n", i*data->wdm->dt, j*data->wdm->df, tdi->X[k]*tdi->X[k], tdi->Y[k]*tdi->Y[k], tdi->Z[k]*tdi->Z[k]);
            }
            fprintf(fptr,"\n");
        }
        fclose(fptr);
    }

    tdi = data->dft;
    sprintf(filename,"%s/dft_data.dat",data->dataDir);
    fptr=fopen(filename,"w");
    for(int i=0; i<data->NFFT; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        switch(data->Nchannel)
        {
            case 1:
                fprintf(fptr,"%.12g %lg %lg\n", f, tdi->X[2*i],tdi->X[2*i+1]);
                break;
            case 2:
                fprintf(fptr,"%.12g %lg %lg %lg %lg\n", f, tdi->A[2*i],tdi->A[2*i+1], tdi->E[2*i],tdi->E[2*i+1]);
                break;
            case 3:
                fprintf(fptr,"%.12g %lg %lg %lg %lg %lg %lg\n", f, tdi->X[2*i],tdi->X[2*i+1], tdi->Y[2*i],tdi->Y[2*i+1], tdi->Z[2*i],tdi->Z[2*i+1]);
                break;
        }
    }
    fclose(fptr);

    if(!strcmp("wavelet",data->basis))
    {
        tdi = data->dwt;
        sprintf(filename,"%s/dwt_data.dat",data->dataDir);
        fptr=fopen(filename,"w");
        for(int i=0; i<data->wdm->NT; i++)
        {
            for(int j=data->lmin; j<data->lmax; j++)
            {
                wavelet_pixel_to_index(data->wdm,i,j,&k);
                k-=data->wdm->kmin;
                fprintf(fptr,"%lg %lg %.14e %.14e %.14e\n", i*data->wdm->dt, j*data->wdm->df,tdi->X[k], tdi->Y[k], tdi->Z[k]);   
            }
            fprintf(fptr,"\n");
        }
        fclose(fptr);
    }
}

void print_wavelet_fourier_spectra(struct Data *data, struct TDI *tdi, char filename[])
{
    int k;
    struct Wavelets *wdm = data->wdm;
    int N = wdm->NF*wdm->NT;
    double T = N*LISA_CADENCE;

    double **freqData = double_matrix(3,N+2); //needs DF & Nyq. bins for norms
    double **waveData = double_matrix(3,N);

    //get TDI data into context
    for(int i=0; i<wdm->NT; i++)
    {
        for(int j=0; j<wdm->NF; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);
            if(k>=wdm->kmin && k<wdm->kmax)
            {
                waveData[0][k] = tdi->X[k-wdm->kmin];
                waveData[1][k] = tdi->Y[k-wdm->kmin];
                waveData[2][k] = tdi->Z[k-wdm->kmin];
            }
        }
    }

    //wavelet to frequency
    wavelet_transform_inverse_freq(wdm, waveData[0], freqData[0]);
    wavelet_transform_inverse_freq(wdm, waveData[1], freqData[1]);
    wavelet_transform_inverse_freq(wdm, waveData[2], freqData[2]);

    FILE *fptr=fopen(filename,"w");
    for(int n=0; n<N/2; n++) 
    {
        double f = (double)n/T;
        if(f>data->fmin && f < data->fmax)
        {
            fprintf(fptr,"%.14e %.14e %.14e %.14e %.14e %.14e %.14e\n", n/T, freqData[0][2*n],freqData[0][2*n+1],freqData[1][2*n],freqData[1][2*n+1], freqData[2][2*n],freqData[2][2*n+1]);
        }
    }
    fclose(fptr);

    free_double_matrix(freqData,3);
    free_double_matrix(waveData,3);
}

void wavelet_layer_to_fourier_transform(struct Data *data)
{
    int k;
    struct Wavelets *wdm = data->wdm;
    int N = wdm->NF*wdm->NT;

    double **freqData = double_matrix(3,N+2); //needs DF & Nyq. bins for norms
    double **waveData = double_matrix(3,N);

    /* store backup DWT data */
    memcpy(data->dwt->X,data->tdi->X,sizeof(double)*data->N);
    memcpy(data->dwt->Y,data->tdi->Y,sizeof(double)*data->N);
    memcpy(data->dwt->Z,data->tdi->Z,sizeof(double)*data->N);
    
    /* populate full DWT data with active layers */
    for(int i=0; i<wdm->NT; i++)
    {
        for(int j=0; j<wdm->NF; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);
            if(k>=wdm->kmin && k<wdm->kmax)
            {
                waveData[0][k] = data->dwt->X[k-wdm->kmin];
                waveData[1][k] = data->dwt->Y[k-wdm->kmin];
                waveData[2][k] = data->dwt->Z[k-wdm->kmin];
            }
        }
    }

    /* out-place DWT->DFT transform */
    wavelet_transform_inverse_freq(wdm, waveData[0], freqData[0]);
    wavelet_transform_inverse_freq(wdm, waveData[1], freqData[1]);
    wavelet_transform_inverse_freq(wdm, waveData[2], freqData[2]);

    /* copy active frequency bins into DFT struct */
    for(int n=0; n<data->NFFT; n++)
    {
        int m = n + data->qmin;
        data->dft->X[2*n]   = freqData[0][2*m];
        data->dft->X[2*n+1] = freqData[0][2*m+1];
        data->dft->Y[2*n]   = freqData[1][2*m];
        data->dft->Y[2*n+1] = freqData[1][2*m+1];
        data->dft->Z[2*n]   = freqData[2][2*m];
        data->dft->Z[2*n+1] = freqData[2][2*m+1];
    }
    
    /* Get A&E channels just in case... */
    for(int n=0; n<data->N; n++) XYZ2AE(data->dft->X[n],data->dft->Y[n],data->dft->Z[n],&data->dft->A[n],&data->dft->E[n]);

    free_double_matrix(freqData,3);
    free_double_matrix(waveData,3);
}

void print_glass_usage()
{
    fprintf(stdout,"\n");
    fprintf(stdout,"=============== GLASS Usage: ============== \n");
    fprintf(stdout,"REQUIRED:\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"OPTIONAL:\n");
    fprintf(stdout,"  -h | --help        : print help message and exit         \n");
    fprintf(stdout,"  -v | --verbose     : enable verbose output               \n");
    fprintf(stdout,"  -q | --quiet       : restrict output                     \n");
    fprintf(stdout,"  -d | --debug       : leaner settings for quick running   \n");
    fprintf(stdout,"\n");
    
    //LISA
    fprintf(stdout,"       =========== LISA =========== \n");
    fprintf(stdout,"       --orbit       : orbit ephemerides file (2.5 GM MLDC)\n");
    fprintf(stdout,"       --channels    : # of channels [1->X,2->AE,3->XYZ](3)\n");
    fprintf(stdout,"       --phase       : phase data (fractional frequency)   \n");
    fprintf(stdout,"       --sangria     : use LDC Sangria TDI conventions     \n");
    fprintf(stdout,"\n");
    
    //Data
    fprintf(stdout,"       =========== Data =========== \n");
    fprintf(stdout,"       --data        : strain data file (ASCII)            \n");
    fprintf(stdout,"       --h5-data     : strain data file (HDF5)             \n");
    fprintf(stdout,"       --h5-no-mbh   : remove mbhs from HDF5 data          \n");
    fprintf(stdout,"       --h5-no-ucb   : remove ucbs from HDF5 data          \n");
    fprintf(stdout,"       --h5-no-ucb-hi: remove high f ucbs from HDF5 data   \n");
    fprintf(stdout,"       --h5-no-vgb   : remove vgbs from HDF5 data          \n");
    fprintf(stdout,"       --h5-no-noise : remove noise from HDF5 data (TODO)  \n");
    fprintf(stdout,"       --psd         : psd data file (ASCII)               \n");
    fprintf(stdout,"       --samples     : number of DFT frequency bins (512)  \n");
    fprintf(stdout,"       --layers      : number of DWT frequency layers (1)  \n");
    fprintf(stdout,"       --padding     : number of bins padded on segment (0)\n");
    fprintf(stdout,"       --start-time  : initial time of epoch  (0)          \n");
    fprintf(stdout,"       --fmin        : minimum frequency                   \n");
    fprintf(stdout,"       --fmax        : maximum frequency (overrides --samples)\n");
    fprintf(stdout,"       --duration    : duration of epoch (31457280)        \n");
    fprintf(stdout,"       --sim-noise   : data w/out noise realization        \n");
    fprintf(stdout,"       --conf-noise  : include model for confusion noise   \n");
    fprintf(stdout,"       --stationary  : use stationary noise model in logL  \n");
    fprintf(stdout,"       --stationary-conf : inject confusion noise as stationary (no time modulation)\n");
    fprintf(stdout,"       --coarse-Q    : WDM time-axis coarse-graining factor (1)\n");
    fprintf(stdout,"       --ws-approx   : freeze Welch-Satterthwaite effective dof per coarse cell (injection only)\n");
    fprintf(stdout,"       --crop-time   : fine time pixels to crop from each end of coarse logL (0)\n");
    fprintf(stdout,"       --noiseseed   : seed for noise RNG                  \n");
    fprintf(stdout,"\n");
    
    //Chain
    fprintf(stdout,"       ========== Chains ========== \n");
    fprintf(stdout,"       --steps       : number of mcmc steps (10000)        \n");
    fprintf(stdout,"       --chainseed   : seed for MCMC RNG                   \n");
    fprintf(stdout,"       --chains      : number of parallel chains (20)      \n");
    fprintf(stdout,"       --no-burnin   : skip burn in steps                  \n");
    fprintf(stdout,"       --resume      : restart from checkpoint             \n");
    fprintf(stdout,"       --threads     : number of parallel threads (max)    \n");
    fprintf(stdout,"       --prior       : sample from prior                   \n");
    fprintf(stdout,"       --no-rj       : turn off RJMCMC                     \n");
    fprintf(stdout,"\n");
    
    //Misc.
    fprintf(stdout,"       =========== Misc =========== \n");
    fprintf(stdout,"       --rundir      : top level run directory ['./']\n");
    fprintf(stdout,"       --match-in1   : input paramaters for overlap [filename] \n");
    fprintf(stdout,"       --match-in2   : output match values [filename] \n");
    fprintf(stdout,"\n");

    /*
    fprintf(stdout,"USEFUL TOBS:\n");
    fprintf(stdout,"   1 wk: %.0f\n",62914560./2./52.);
    fprintf(stdout,"   1 mo: %.0f \n",62914560./2./12.);
    fprintf(stdout,"   2 mo: %.0f \n",62914560./2./6.);
    fprintf(stdout,"   1 yr: %.0f (default)\n",62914560./2.);
    fprintf(stdout,"   2 yr: %.0f \n",62914560.);
    fprintf(stdout,"   5 yr: %.0f \n",5.*62914560./2.);
    fprintf(stdout,"  10 yr: %.0f \n",10.*62914560./2.);
     */
    fprintf(stdout,"\n");
}

void parse_data_args(int argc, char **argv, struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, char basis[])
{
    //copy argv since getopt permutes order
    char **argv_copy=malloc((argc+1) * sizeof *argv_copy);
    copy_argv(argc,argv,argv_copy);
    opterr=0; //suppress warnings about unknown arguments
    
    //Set defaults
    flags->rj          = 1;
    flags->help        = 0;
    flags->calibration = 0;
    flags->verbose     = 0;
    flags->quiet       = 0;
    flags->simNoise    = 0;
    flags->confNoise   = 0;
    flags->sgwbTemplate=-1;
    flags->stationary  = 0;
    flags->stationaryConf = 0;
    flags->coarseQ     = 1;
    flags->ws_approx   = 0;
    flags->crop_time   = 0;
    flags->burnin      = 1;
    flags->debug       = 0;
    flags->strainData  = 0;
    flags->timeseries  = false;
    flags->hdf5Data    = 0;
    flags->psd         = 0;
    flags->orbit       = 0;
    flags->prior       = 0;
    flags->resume      = 0;
    flags->NMCMC       = 1000;
    flags->NBURN       = 1000;
    flags->threads     = 12;
    sprintf(flags->runDir,"./");
    chain->NC          = 12;//number of chains
    int set_fmax_flag  = 0; //flag watching for if fmax is set by CLI
    
    /* Simulated data building blocks */
    flags->no_mbh = 0;
    flags->no_ucb = 0;
    flags->no_vgb = 0;
    flags->no_noise = 0;
    
    /*
     default data format is 'phase'
     optional support for 'frequency' a la LDCs
     */
    sprintf(data->format,"sangria");
    sprintf(data->basis,"%s",basis);

    
    data->cseed = 150914;
    data->nseed = 151226;
    data->iseed = 151012;


    //Specifying the expected options
    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"samples",      required_argument, 0, 0},
        {"layers",       required_argument, 0, 0},
        {"padding",      required_argument, 0, 0},
        {"duration",     required_argument, 0, 0},
        {"start-time",   required_argument, 0, 0},
        {"orbit",        required_argument, 0, 0},
        {"chains",       required_argument, 0, 0},
        {"chainseed",    required_argument, 0, 0},
        {"noiseseed",    required_argument, 0, 0},
        {"data",         required_argument, 0, 0},
        {"h5-data",      required_argument, 0, 0},
        {"psd",          required_argument, 0, 0},
        {"fmin",         required_argument, 0, 0},
        {"fmax",         required_argument, 0, 0},
        {"channels",     required_argument, 0, 0},
        {"steps",        required_argument, 0, 0},
        {"threads",      required_argument, 0, 0},
        {"rundir",       required_argument, 0, 0},
        {"sgwb-template",required_argument, 0, 0},
        
        /* These options don’t set a flag.
         We distinguish them by their indices. */
        {"help",        no_argument, 0,'h'},
        {"verbose",     no_argument, 0,'v'},
        {"quiet",       no_argument, 0,'q'},
        {"debug",       no_argument, 0,'d'},
        {"resume",      no_argument, 0, 0 },
        {"sim-noise",   no_argument, 0, 0 },
        {"conf-noise",  no_argument, 0, 0 },
        {"stationary",  no_argument, 0, 0 },
        {"stationary-conf", no_argument, 0, 0 },
        {"coarse-Q",    required_argument, 0, 0 },
        {"ws-approx",   no_argument, 0, 0 },
        {"crop-time",   required_argument, 0, 0 },
        {"phase",       no_argument, 0, 0 },
        {"sangria",     no_argument, 0, 0 },
        {"prior",       no_argument, 0, 0 },
        {"no-burnin",   no_argument, 0, 0 },
        {"no-rj",       no_argument, 0, 0 },
        {"calibration", no_argument, 0, 0 },
        {"h5-no-mbh",   no_argument, 0, 0 },
        {"h5-no-ucb",   no_argument, 0, 0 },
        {"h5-no-ucb-hi",no_argument, 0, 0 },
        {"h5-no-vgb",   no_argument, 0, 0 },
        {"h5-no-noise", no_argument, 0, 0 },
        {"timeseries",  no_argument, 0, 0 },
        {0, 0, 0, 0}
    };
    
    int opt=0;
    int long_index=0;
    
    //Loop through argv string and pluck out arguments
    while ((opt = getopt_long_only(argc, argv_copy,"apl:b:", long_options, &long_index )) != -1)
    {
        switch (opt)
        {
                
            case 0:
                if(strcmp("samples",     long_options[long_index].name) == 0) data->NFFT        = atoi(optarg);
                if(strcmp("layers",      long_options[long_index].name) == 0) data->Nlayer      = atoi(optarg);
                if(strcmp("padding",     long_options[long_index].name) == 0) data->qpad        = atoi(optarg);
                if(strcmp("start-time",  long_options[long_index].name) == 0) data->t0          = (double)atof(optarg);
                if(strcmp("chains",      long_options[long_index].name) == 0) chain->NC         = atoi(optarg);
                if(strcmp("chainseed",   long_options[long_index].name) == 0) data->cseed       = (unsigned int)atoi(optarg);
                if(strcmp("noiseseed",   long_options[long_index].name) == 0) data->nseed       = (unsigned int)atoi(optarg);
                if(strcmp("injseed",     long_options[long_index].name) == 0) data->iseed       = (unsigned int)atoi(optarg);
                if(strcmp("sim-noise",   long_options[long_index].name) == 0) flags->simNoise   = 1;
                if(strcmp("conf-noise",  long_options[long_index].name) == 0) flags->confNoise  = 1;
                if(strcmp("sgwb-template",long_options[long_index].name)== 0) flags->sgwbTemplate=atoi(optarg);
                if(strcmp("stationary",  long_options[long_index].name) == 0) flags->stationary = 1;
                if(strcmp("stationary-conf", long_options[long_index].name) == 0) flags->stationaryConf = 1;
                if(strcmp("coarse-Q",    long_options[long_index].name) == 0) flags->coarseQ    = atoi(optarg);
                if(strcmp("ws-approx",   long_options[long_index].name) == 0) flags->ws_approx  = 1;
                if(strcmp("crop-time",   long_options[long_index].name) == 0) flags->crop_time  = atoi(optarg);
                if(strcmp("prior",       long_options[long_index].name) == 0) flags->prior      = 1;
                if(strcmp("no-burnin",   long_options[long_index].name) == 0) flags->burnin     = 0;
                if(strcmp("no-rj",       long_options[long_index].name) == 0) flags->rj         = 0;
                if(strcmp("calibration", long_options[long_index].name) == 0) flags->calibration= 1;
                if(strcmp("resume",      long_options[long_index].name) == 0) flags->resume     = 1;
                if(strcmp("h5-no-mbh",   long_options[long_index].name) == 0) flags->no_mbh     = 1;
                if(strcmp("h5-no-ucb",   long_options[long_index].name) == 0) flags->no_ucb     = 1;
                if(strcmp("h5-no-vgb",   long_options[long_index].name) == 0) flags->no_vgb     = 1;
                if(strcmp("h5-no-ucb-hi",long_options[long_index].name) == 0) flags->no_ucb_hi  = 1;
                if(strcmp("h5-no-noise", long_options[long_index].name) == 0) flags->no_noise   = 1;
                if(strcmp("timeseries" , long_options[long_index].name) == 0) flags->timeseries = true;
                if(strcmp("threads",     long_options[long_index].name) == 0) flags->threads    = atoi(optarg);
                if(strcmp("rundir",      long_options[long_index].name) == 0) strcpy(flags->runDir,optarg);
                if(strcmp("phase",       long_options[long_index].name) == 0) sprintf(data->format,"phase");
                if(strcmp("sangria",     long_options[long_index].name) == 0) sprintf(data->format,"sangria");
                if(strcmp("fmin",        long_options[long_index].name) == 0) sscanf(optarg, "%lg", &data->fmin);
                if(strcmp("fmax",        long_options[long_index].name) == 0)
                {
                    set_fmax_flag = 1;
                    sscanf(optarg, "%lg", &data->fmax);
                }
                if(strcmp("duration",    long_options[long_index].name) == 0)
                {   
                    data->T   = (double)atof(optarg);
                    if(!strcmp(data->basis,"wavelet")) data->T = floor(data->T/WAVELET_DURATION)*WAVELET_DURATION;
                    data->sqT = sqrt(data->T);
                }
                if(strcmp("steps",       long_options[long_index].name) == 0)
                {
                    flags->NMCMC = atoi(optarg);
                    flags->NBURN = flags->NMCMC;
                }
                if(strcmp("data", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->strainData = 1;
                    sprintf(data->fileName,"%s",optarg);
                }
                if(strcmp("h5-data", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->hdf5Data = 1;
                    flags->strainData = 1;
                    sprintf(data->fileName,"%s",optarg);
                }
                if(strcmp("psd", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->psd = 1;
                    sprintf(flags->psdFile,"%s",optarg);
                }
                if(strcmp("orbit", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->orbit = 1;
                    sprintf(orbit->OrbitFileName,"%s",optarg);
                }
                if(strcmp("channels",long_options[long_index].name) == 0)
                {
                    data->Nchannel = (int)atoi(optarg);
                    if(data->Nchannel<1 || data->Nchannel>3)
                    {
                        fprintf(stderr,"Requested umber of channels (%i) not supported\n",data->Nchannel);
                        fprintf(stderr,"Use --channels 1 for X (Michelson) data\n");
                        fprintf(stderr,"    --channels 2 for AE data\n");
                        fprintf(stderr,"    --channels 3 for XYZ data\n");
                        exit(1);
                    }
                }
                break;
            case 'd' : flags->debug = 1;
                break;
            case 'h' : flags->help = 1;
                break;
            case 'v' : flags->verbose = 1;
                break;
            case 'q' : flags->quiet = 1;
                break;
            default:
                break;
        }
    }
    if(flags->cheat || !flags->burnin) flags->NBURN = 0;

    if(flags->sgwbTemplate < -1 || flags->sgwbTemplate >= SGWB_TEMPLATE_COUNT)
    {
        fprintf(stderr,"selected SGWB template does not exist, options are:\n");
        for (int i=0; i<SGWB_TEMPLATE_COUNT; i++) {
            fprintf(stderr,"\t%d: %s\n",i,SGWB_TEMPLATE_NAMES[i]);
        }
        exit(1);
    }
    
    if(flags->verbose && flags->quiet)
    {
        fprintf(stderr,"--verbose and --quiet flags are in conflict\n");
        exit(1);
    }
    
    //Chains should be a multiple of threads for best usage of cores
    if(chain->NC % flags->threads !=0){
        chain->NC += flags->threads - (chain->NC % flags->threads);
    }
    
    omp_set_num_threads(flags->threads);

    //override size of data if fmax was requested
    if(set_fmax_flag)
    {
        if(!strcmp(data->basis,"fourier")) data->NFFT = (int)floor((data->fmax - data->fmin)*data->T); //number of frequency bins
        if(!strcmp(data->basis,"wavelet")) data->Nlayer = (int)( (ceil(data->fmax/WAVELET_BANDWIDTH) - floor(data->fmin/WAVELET_BANDWIDTH)) ); //number of frequency layers
    }

    //pad data
    if(!strcmp(data->basis,"fourier"))
    {
        data->NFFT += 2*data->qpad;
        data->fmin -= data->qpad/data->T;
    }
    if(!strcmp(data->basis,"wavelet")) 
    {
        data->Nlayer += 2; // pad with layers above and below
        data->fmin -= WAVELET_BANDWIDTH;
    }

    
    
    //map fmin to nearest bin
    if(!strcmp(data->basis,"fourier"))
    {
        data->fmin = floor(data->fmin*data->T)/data->T;
        data->fmax = data->fmin + (double)data->NFFT/data->T;
    }
    if(!strcmp(data->basis,"wavelet"))
    {
        data->fmin = floor(data->fmin/WAVELET_BANDWIDTH)*WAVELET_BANDWIDTH;
        data->fmax = data->fmin + (double)(data->Nlayer)*WAVELET_BANDWIDTH;

    }
    //printf("old fmin=%lg, fmax=%lg\n",data->fmin,data->fmax);

    //after all of that resize data
    if(!strcmp(data->basis,"fourier")) data->N = data->NFFT*2;
    if(!strcmp(data->basis,"wavelet")) 
    {
        data->N = data->Nlayer*(int)floor(data->T / WAVELET_DURATION);
        data->NFFT = data->N/2.;
    }

    //Print version control
//    sprintf(filename,"glass.log");
//    FILE *runlog = fopen(filename,"w");
//    print_version(runlog);
    
    //Report on set parameters
//    if(!flags->quiet) print_run_settings(argc, argv, data, orbit, flags, stdout);
//    print_run_settings(argc, argv, data, orbit, flags, runlog);
    
//    fclose(runlog);

    //reset opt counter
    optind = 0;

    // double check wavelet fmin
    if(!strcmp(data->basis,"wavelet")) 
    {
        if (floor(data->fmin/WAVELET_BANDWIDTH) == 0) {
            fprintf(stderr, "Error, your minimum frequency (after padding) includes the DC wavelet frequency bin!\nThis bin goes until %g Hz\nYou almost certainly don't want to do this. Recommended fix: increase fmin\n", WAVELET_BANDWIDTH);
            exit(-2);
        }
    }

    //free placeholder for argvs
    for(int i=0; i<=argc; i++)free(argv_copy[i]);
    free(argv_copy);

}

void copy_argv(int argc, char **argv, char **new_argv)
{
    for(int i = 0; i < argc; ++i)
    {
        size_t length = strlen(argv[i])+1;
        new_argv[i] = malloc(length);
        memcpy(new_argv[i], argv[i], length);
    }
    new_argv[argc] = NULL;
}

int checkfile(char filename[])
{
    FILE *fptr = fopen(filename, "r");
    if(fptr)
    {
        fclose(fptr);
        return 1;
    }
    else
    {
        fprintf(stderr,"File %s does not exist\n",filename);
        fprintf(stderr,"\n");
        exit(1);
    }
}
