/*
 *  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Kristen Lackeos, Neil J. Cornish
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

/**
 @file ucb_catalog.c
 \brief App for post-processing UCB/VGB sampler output for catalog generation
 */

#include <glass_utils.h>
#include <glass_ucb.h>

void print_usage_catalog()
{
    fprintf(stdout,"\n");
    fprintf(stdout,"================ UCB CATALOG: ============== \n");
    fprintf(stdout,"REQUIRED:\n");
    fprintf(stdout,"       --chain-file  : chain file to be sorted into catalog\n");
    fprintf(stdout,"       --sources     : maximum number of sources (10)      \n");
    fprintf(stdout,"       --fmin        : minimum frequency                   \n");
    fprintf(stdout,"\n");
    fprintf(stdout,"OPTIONAL:\n");
    fprintf(stdout,"  -h | --help        : print help message and exit         \n");
    fprintf(stdout,"       --orbit       : orbit ephemerides file (2.5 GM MLDC)\n");
    fprintf(stdout,"       --samples     : number of frequency bins (2048)     \n");
    fprintf(stdout,"       --padding     : number of bins padded on segment (0)\n");
    fprintf(stdout,"       --duration    : duration of time segment (62914560) \n");
    fprintf(stdout,"       --start-time  : initial time of epoch  (0)          \n");
    fprintf(stdout,"       --match       : match threshold for waveforms (0.8) \n");
    fprintf(stdout,"       --phase       : fractional frequency data (sangria) \n");
    fprintf(stdout,"       --sangria     : phase/tdi conventions for LDC2      \n");
    fprintf(stdout,"       --f-double-dot: include f double dot in model       \n");
    fprintf(stdout,"       --channels    : TDI channels [1->X,2->AE,3->XYZ] (2)\n");
    fprintf(stdout,"       --noise-file  : reconstructed noise model           \n");
    fprintf(stdout,"                       e.g., data/power_noise_t0_f0.dat    \n");
    fprintf(stdout,"       --catalog     : list of known sources               \n");
    fprintf(stdout,"       --Tcatalog    : observing time of previous catalog  \n");
    fprintf(stdout,"       --Nmode       : max number of GMM modes (16)        \n");
    fprintf(stdout,"       --thin        : factor for thinning chains in GMM   \n");
    fprintf(stdout,"\n");
    fprintf(stdout,"EXAMPLE:\n");
    fprintf(stdout,"glass_ucb_catalog --fmin 0.004 --samples 256 --duration 31457280 --sources 5 --chain-file chains/dimension_chain.dat.5");
    fprintf(stdout,"\n");
    fprintf(stdout,"\n");
    exit(EXIT_FAILURE);
}

void parse_catalog(int argc, char **argv, struct Data *data, struct Orbit *orbit, struct Flags *flags, int Nmax, double *Tcatalog, size_t *NMODE, size_t *NTHIN)
{
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    
    if(argc<4) print_usage_catalog();
    
    int DMAX_default = 10;
    
    //Set defaults
    flags->emPrior     = 0;
    flags->orbit       = 0;
    flags->match       = 0;
    flags->DMAX        = DMAX_default;
    flags->catalog     = 0;
    flags->verbose     = 0;
    flags->simNoise    = 1; //hijack simNoise flag for noise model
    data->pmax         = 0.5; //default match tolerance for inclusion (hijacked pmax in data structure)
    *Tcatalog          = -1.0;
    
    /*
     default data format is fractional frequency
     optional support for phase a la pre LDCs
     */
    sprintf(data->format,"frequency");
            
    data->T        = 31457280; /* two "mldc years" at 15s sampling */
    data->t0       = 0.0;
    data->NFFT     = 512;
    data->N        = data->NFFT*2;
    data->Nchannel = 2; //1=X, 2=AE
    
    data->cseed = 150914;
    data->nseed = 151226;
    data->iseed = 151012;
    
    
    
    //Specifying the expected options
    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"samples",   required_argument, 0, 0},
        {"padding",   required_argument, 0, 0},
        {"start-time",required_argument, 0, 0},
        {"duration",  required_argument, 0, 0},
        {"segments",  required_argument, 0, 0},
        {"sources",   required_argument, 0, 0},
        {"orbit",     required_argument, 0, 0},
        {"fmin",      required_argument, 0, 0},
        {"channels",  required_argument, 0, 0},
        {"chain-file",required_argument, 0, 0},
        {"noise-file",required_argument, 0, 0},
        {"match",     required_argument, 0, 0},
        {"catalog",   required_argument, 0, 0},
        {"Tcatalog",  required_argument, 0, 0},
        {"Nmode",     required_argument, 0, 0},
        {"thin",      required_argument, 0, 0},
        
        /* These options don’t set a flag.
         We distinguish them by their indices. */
        {"help",        no_argument, 0,'h'},
        {"verbose",     no_argument, 0,'v'},
        {"phase",       no_argument, 0, 0 },
        {"sangria",     no_argument, 0, 0 },
        {"f-double-dot",no_argument, 0, 0 },
        {0, 0, 0, 0}
    };
    
    int opt=0;
    int long_index=0;
    
    //Loop through argv string and pluck out arguments
    while ((opt = getopt_long_only(argc, argv,"apl:b:", long_options, &long_index )) != -1)
    {
        switch (opt)
        {
                
            case 0:
                if(strcmp("samples",     long_options[long_index].name) == 0) data->NFFT = atoi(optarg);
                if(strcmp("padding",     long_options[long_index].name) == 0) data->qpad     = atoi(optarg);
                if(strcmp("start-time",  long_options[long_index].name) == 0) data->t0       = (double)atof(optarg);
                if(strcmp("duration",    long_options[long_index].name) == 0) data->T        = (double)atof(optarg);
                if(strcmp("fmin",        long_options[long_index].name) == 0) data->fmin     = (double)atof(optarg);
                if(strcmp("Tcatalog",    long_options[long_index].name) == 0) *Tcatalog      = (double)atof(optarg);
                if(strcmp("Nmode",       long_options[long_index].name) == 0) *NMODE         = (size_t)atoi(optarg);
                if(strcmp("thin",        long_options[long_index].name) == 0) *NTHIN         = (size_t)atoi(optarg);
                if(strcmp("channels",    long_options[long_index].name) == 0) data->Nchannel = (int)atoi(optarg);
                if(strcmp("match",       long_options[long_index].name) == 0) data->pmax     = atof(optarg);
                if(strcmp("sources",     long_options[long_index].name) == 0) flags->DMAX    = atoi(optarg);
                //if(strcmp("f-double-dot",long_options[long_index].name) == 0) data->NP   = 9; TODO: get NP back into CLI

                if(strcmp("phase",   long_options[long_index].name) == 0) sprintf(data->format,"phase");
                if(strcmp("sangria", long_options[long_index].name) == 0) sprintf(data->format,"sangria");
                
                if(strcmp("orbit", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->orbit = 1;
                    sprintf(orbit->OrbitFileName,"%s",optarg);
                }
                if(strcmp("noise-file", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    sprintf(flags->noiseFile,"%s",optarg);
                    flags->simNoise = 0;
                }
                if(strcmp("chain-file", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    sprintf(data->fileName,"%s",optarg);
                }
                if(strcmp("catalog",long_options[long_index].name) == 0)
                {
                    flags->catalog = 1;
                    sprintf(flags->catalogFile,"%s",optarg);
                }
                
                break;
            case 'h' :
                print_usage_catalog();
                exit(EXIT_FAILURE);
                break;
            case 'v' : flags->verbose = 1;
                break;
            default: print_usage_catalog();
                exit(EXIT_FAILURE);
        }
    }
    
    //check if catalog is supplied so is corresponding duration
    if(flags->catalog && *Tcatalog<0)
    {
        fprintf(stderr,"Use of --catalog flag requires --Tcatalog argument\n");
        fprintf(stderr,"Supply observing time for input catalog or remove the option\n");
        exit(1);
    }
    
    //pad data
    data->NFFT += 2*data->qpad;
    data->N = data->NFFT*2;
    data->fmin -= data->qpad/data->T;
    
    //map fmin to nearest bin
    data->fmin = floor(data->fmin*data->T)/data->T;
    
    //Print command line
    char filename[128];
    sprintf(filename,"catalog_%i.sh",flags->DMAX);
    FILE *out = fopen(filename,"w");
    fprintf(out,"#!/bin/sh\n\n");
    for(opt=0; opt<argc; opt++) fprintf(out,"%s ",argv[opt]);
    fprintf(out,"\n\n");
    fclose(out);
    
    //Print version control
    FILE *runlog = fopen("gb_catalog.log","w");
    print_version(runlog);
    
    fclose(runlog);
}


static int safe_scan_source_params(struct Data *data, struct Source *source, FILE *fptr)
{

    int check = 0;
    check+=fscanf(fptr,"%lg",&source->f0);
    check+=fscanf(fptr,"%lg",&source->dfdt);
    check+=fscanf(fptr,"%lg",&source->amp);
    check+=fscanf(fptr,"%lg",&source->phi);
    check+=fscanf(fptr,"%lg",&source->costheta);
    check+=fscanf(fptr,"%lg",&source->cosi);
    check+=fscanf(fptr,"%lg",&source->psi);
    check+=fscanf(fptr,"%lg",&source->phi0);
    if(UCB_MODEL_NP>8)
        check+=fscanf(fptr,"%lg",&source->d2fdt2);
    
    if(!check)
    {
        return 1;
    }
    
    //map to parameter names (just to make code readable)
    map_params_to_array(source, source->params, data->T);
    return 0;

}

static void source_waveform_wrapper(struct Source *source, struct Data *data, struct Orbit *orbit)
{
    source->tdi = malloc(sizeof(struct TDI));
    alloc_tdi(source->tdi,data->N, data->Nchannel);
    ucb_alignment(orbit, data, source);
    ucb_waveform(orbit, data->format, data->T, data->t0, source->params, UCB_MODEL_NP, source->tdi->X, source->tdi->Y, source->tdi->Z ,source->tdi->A, source->tdi->E, source->BW, data->Nchannel);
}

int main(int argc, char *argv[])
{
    fprintf(stdout, "\n================ UCB CATALOG ================\n");

    /* ************************************************************** */
    /*             Allocate & Initialize Data Structures              */
    /* ************************************************************** */
    
    int NTEMP = 1;     //needed size of data structure
    size_t NMODE = 16; //default size of GMM
    size_t NTHIN = 1;  //thinning rate of chain
    int check;
    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Data  *data  = malloc(sizeof(struct Data));
    struct Data *data_old = malloc(sizeof(struct Data));
    
    double Tcatalog; //duration of previous analysis (for source heritage)
    
    /* Parse command line and set defaults/flags */
    parse_catalog(argc,argv,data,orbit,flags,NTEMP,&Tcatalog,&NMODE,&NTHIN);
    alloc_data(data, flags);
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin + data->NFFT;
    
    data_old->T    = Tcatalog;
    data_old->qmin = (int)(data->fmin*data_old->T);
    data_old->qmax = data_old->qmin + data->NFFT;
    
    //File containing chain samples
    FILE *chain_file = fopen(data->fileName,"r");
    
    //Orbits
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
            free_orbit(orbit);
            return(1);
            break;
    }
    
    // Get priors
    /* Load priors as used in MCMC */
    struct Model *model = malloc(sizeof(struct Model));
    alloc_model(data,model,flags->DMAX);
    set_uniform_prior(flags, model, data, 1);
    
    /* Reformat for catalog */
    //frequency & derivatives
    model->prior[0][0] /= data->T;
    model->prior[0][1] /= data->T;
    if(UCB_MODEL_NP>7)
    {
        model->prior[7][0] /= data->T*data->T;
        model->prior[7][1] /= data->T*data->T;
    }
    if(UCB_MODEL_NP>8)
    {
        model->prior[8][0] /= data->T*data->T*data->T;
        model->prior[8][1] /= data->T*data->T*data->T;
    }
    
    //alias for catalog->entry pointers used later on
    struct Entry *entry = NULL;
    
    struct Source *sample = NULL;
    sample = malloc(sizeof *sample);
    alloc_source(sample, data->N, data->Nchannel);
    
    
    //count lines in chain file
    int N=0;

    char *line;
    char buffer[16384];

    while( (line=fgets(buffer, 16384, chain_file)) != NULL) N++;

    rewind(chain_file);
    
    
    //selection criteria for catalog entries
    int matchFlag;          //track if sample matches any entries
    double Match;     //match between pairs of waveforms
    double Distance;  //distance between pairs of waveforms
    double tolerance = data->pmax; //tolerance on match to be considered associated
    double dqmax = 10;      //maximum frequency separation to try match calculation (in frequency bins)
    int downsample = 1;
    
    
    //Book-keeping
    int DMAX = flags->DMAX; //maximum number of sources per chain sample (needs to be read in from above)
    int IMAX = N/DMAX;     //maximum number of chain samples (needs to be read in from above)
    int NMAX = N; //(absurd) upper limit on total number of catalog entries
    
    struct Catalog *catalog = NULL;
    catalog = malloc(sizeof(struct Catalog));
    catalog->N = 0; //start with 0 sources in catalog
    catalog->entry = malloc(NMAX*sizeof(struct Entry*));
    
    /* ************************************************************** */
    /*             Allocate & Initialize Instrument Model             */
    /* ************************************************************** */
    
    struct Noise *noise = NULL;
    noise = malloc(sizeof(struct Noise));
    alloc_noise(noise, data->NFFT, data->Nlayer, data->Nchannel);
    
    //Noise model
    //Get noise spectrum for data segment
    if(flags->simNoise)
    {
        double Sop, Spm;
        for(int n=0; n<data->N; n++)
        {
            double f = data->fmin + (double)(n)/data->T;
            get_noise_levels("sangria", f, &Spm, &Sop);
            if(strcmp(data->format,"phase")==0)
            {
                noise->C[0][0][n] = AEnoise(orbit->L, orbit->fstar, f);
                noise->C[1][1][n] = AEnoise(orbit->L, orbit->fstar, f);
            }
            else if(strcmp(data->format,"frequency")==0)
            {
                noise->C[0][0][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                noise->C[1][1][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
            }
            else if(strcmp(data->format,"sangria")==0)
            {
                noise->C[0][0][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop)/sqrt(2.);
                noise->C[1][1][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop)/sqrt(2.);
            }
            else
            {
                fprintf(stderr,"Unsupported data format %s",data->format);
                exit(1);
            }
        }
        invert_noise_covariance_matrix(noise);
    }
    else
    {
        double junk;
        FILE *noiseFile = fopen(flags->noiseFile,"r");
        for(int n=0; n<data->N; n++)
        {
            check = 0;
            check += fscanf(noiseFile,"%lg",&junk); //);f
            check += fscanf(noiseFile,"%lg",&noise->C[0][0][n]);//A_med);
            check += fscanf(noiseFile,"%lg",&noise->C[1][1][n]);//E_med);
            if(!check)
            {
                fprintf(stderr,"Error reading %s\n",flags->noiseFile);
                exit(1);
            }
        }
        invert_noise_covariance_matrix(noise);
        fclose(noiseFile);
    }
    
    /* **************************************************************** */
    /*        First sample of the chain initializes entry list          */
    /* **************************************************************** */
    for(int d=0; d<DMAX; d++)
    {
        
        //parse source in first sample of chain file
        check = safe_scan_source_params(data, sample, chain_file);
        if(!check)
        {
            //Book-keeping of waveform in time-frequency volume
            ucb_alignment(orbit, data, sample);
            
            //calculate waveform model of sample
            ucb_waveform(orbit, data->format, data->T, data->t0, sample->params, UCB_MODEL_NP, sample->tdi->X, sample->tdi->Y, sample->tdi->Z, sample->tdi->A, sample->tdi->E, sample->BW, data->Nchannel);
            
            //check frequencies
            int q_sample = (int)floor(sample->f0 * data->T);
            if(q_sample > data->qmin+data->qpad && q_sample < data->qmax-data->qpad)
            {
                //add new source to catalog
                create_new_source(catalog, sample, noise, 0, IMAX, data->N, sample->tdi->Nchannel);
            }
        }
    }
    
    
    /* ****************************************************************/
    /*            Now loop over the rest of the chain file            */
    /* ****************************************************************/
    
    //prevent multiple templates being added to the same source (only relevent for low SNR, low match threshold)
    int *entryFlag = malloc(NMAX*sizeof(int));
    
    fprintf(stdout,"\nLooping over chain file\n");
    for(int i=1; i<IMAX; i++)
    {
        if(IMAX>100 && i%(IMAX/100)==0)printProgress((double)i/(double)IMAX);
        
        for(int n=0; n<N; n++) entryFlag[n] = 0;
        
        //check each source in chain sample
        for(int d=0; d<DMAX; d++)
        {
            //parse source parameters
            check = safe_scan_source_params(data, sample, chain_file);
            
            if(!check)
            {
                //find where the source fits in the measurement band
                ucb_alignment(orbit, data, sample);
                
                //calculate waveform model of sample
                ucb_waveform(orbit, data->format, data->T, data->t0, sample->params, UCB_MODEL_NP, sample->tdi->X, sample->tdi->Y, sample->tdi->Z, sample->tdi->A, sample->tdi->E, sample->BW, data->Nchannel);
                
                double q_sample = sample->f0 * data->T;
                
                if(q_sample < data->qmin+data->qpad || q_sample > data->qmax-data->qpad) continue;
                
                if(i%downsample!=0) continue;
                
                //calculate match of sample and all entries
                matchFlag = 0;
                for(int n=0; n<catalog->N; n++)
                {
                    entry = catalog->entry[n];
                    
                    //check frequency separation
                    double q_entry  = entry->source[0]->f0 * data->T;
                    
                    if( fabs(q_entry-q_sample) > dqmax ) Match = -1.0;
                    
                    //calculate match
                    else
                    {
                        Match = waveform_match(sample, entry->source[0], noise);
                    }
                    
                    if(Match > tolerance && !entryFlag[n])
                    {
                        matchFlag = 1;
                        entryFlag[n] = 1;
                        Distance = waveform_distance(sample, entry->source[0], noise);
                        //append sample to entry
                        entry->match[entry->Nchain] = Match;
                        entry->distance[entry->Nchain] = Distance;
                        entry->stepFlag[i] = 1;
                        append_sample_to_entry(entry, sample, IMAX, data->N, data->Nchannel);
                        
                        //stop looping over entries in catalog
                        break;
                    }
                    
                }//end loop over catalog entries
                
                
                //if the match tolerence is never met, add as new source
                if(!matchFlag)
                {
                    entryFlag[catalog->N]=1;
                    create_new_source(catalog, sample, noise, i, IMAX, data->N, data->Nchannel);
                }
                
            }
        }//end loop over sources in chain sample
        
    }//end loop over chain
    fprintf(stdout,"\n");
    fclose(chain_file);
    free(entryFlag);
    
    
    /* ****************************************************************/
    /*             Select entries that have enough weight             */
    /* ****************************************************************/
    double weight;
    double weight_threshold = 0.5; //what fraction of samples must an entry have to count?
    
    int detections = 0;            //number of entries that meet the weight threshold
    int *detection_index = malloc(NMAX*sizeof(int));//detection_index[NMAX];     //list of entry indicies for detected sources
    
    
    for(int n=0; n<catalog->N; n++)
    {
        
        weight = (double)catalog->entry[n]->Nchain/(double)(IMAX/downsample);
        if(weight > weight_threshold)
        {
            //record index of catalog entry for detected source
            detection_index[detections] = n;
            
            //increment the number of detections
            detections++;
        }
    }
    fprintf(stdout,"\nNumber of discrete sources is %d.\n",detections);
    
    
    /* *************************************************************** */
    /*           Format selected entries as L3 data products           */
    /* *************************************************************** */
    
    char outdir[MAXSTRINGSIZE];
    sprintf(outdir,"catalog_%i",DMAX);
    mkdir(outdir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    char filename[128];
    sprintf(filename,"%s/entries.dat",outdir);
    FILE *catalogFile = fopen(filename,"w");
    
    double f_med;
    double *f_vec;
    int *index;
    int i_med;
    
    fprintf(stdout,"\nPost processing events\n");
    for(int d=0; d<detections; d++)
    {
        printProgress((double)(d+1)/(double)detections);
        
        int n = detection_index[d];
        entry = catalog->entry[n];
        
        //get sample containing median frequency as identifier of source
        f_vec = calloc(entry->Nchain,sizeof(double));
        for(int i=0; i<entry->Nchain; i++) f_vec[i] = entry->source[i]->f0;
        
        index = calloc(entry->Nchain,(sizeof(int)));
        index_sort(index, f_vec, entry->Nchain);

        i_med = (int)index[entry->Nchain/2];
        free(f_vec);
        free(index);
        
        f_med = entry->source[i_med]->f0;//gsl_stats_median_from_sorted_data(f_vec, 1, entry->Nchain);
        
        entry->i = i_med;
        
        // get waveform of median sample
        source_waveform_wrapper(entry->source[i_med], data, orbit);
        
        //replace stored SNR with median sample
        entry->SNR = snr(entry->source[i_med],noise);
        
        //name source based on median frequency
        sprintf(entry->name,"LDC%010li",(long)(f_med*1e10));
        
        //print point estimate parameters at median frequency
        sprintf(filename, "%s/%s_params.dat", outdir, entry->name);
        FILE *out = fopen(filename, "w");
        print_source_params(data,entry->source[i_med],out);
        fprintf(out,"\n");
        fclose(out);
        
        sprintf(filename, "%s/%s_waveform.dat", outdir,entry->name);
        out = fopen( filename, "w");
        
        struct Source *b=entry->source[0];
        int N = b->tdi->N;
        double *b_A = malloc(N*sizeof(double));
        double *b_E = malloc(N*sizeof(double));
        for(int i=0; i<N; i++)
        {
            b_A[i] = 0.0;
            b_E[i] = 0.0;
        }
        int qmin = b->qmin - b->imin;
        
        for(int i=0; i<b->BW; i++)
        {
            int j = i+b->qmin-qmin;
            if(j>-1 && j<N/2)
            {
                int i_re = 2*i;
                int i_im = i_re+1;
                int j_re = 2*j;
                int j_im = j_re+1;
                
                b_A[j_re] = b->tdi->A[i_re];
                b_A[j_im] = b->tdi->A[i_im];
                b_E[j_re] = b->tdi->E[i_re];
                b_E[j_im] = b->tdi->E[i_im];
                double f = (double)(j+data->qmin)/data->T;
                fprintf(out,"%.12g %lg %lg %lg %lg",f,b_A[j_re],b_A[j_im],b_E[j_re],b_E[j_im]);
                fprintf(out,"\n");
            }//check that index is in range
        }//loop over waveform bins
        
        free_tdi(entry->source[i_med]->tdi);
        free(b_A);
        free(b_E);
        fclose(out);
        
        //evidence for source related to number of visits in the chain
        entry->evidence = (double)(entry->Nchain-1)/(double)(IMAX/downsample);
        
        fprintf(catalogFile,"%s %lg %lg\n",entry->name, entry->SNR, entry->evidence);
    }
    fflush(catalogFile);
    fprintf(stdout,"\n");
    
    printf("get correlation matrix\n");
    double **corr=NULL;
    corr = malloc(detections*UCB_MODEL_NP*sizeof(double *));
    for(int n=0; n<detections*UCB_MODEL_NP; n++)corr[n] = calloc(detections*UCB_MODEL_NP,sizeof(double));

    get_correlation_matrix(data, catalog, detection_index, detections, IMAX, corr);

    sprintf(filename,"%s/correlation_matrix.dat",outdir);
    FILE *corrFile = fopen(filename,"w");
    for(int n=0; n<detections*UCB_MODEL_NP; n++)
    {
        for(int m=0; m<detections*UCB_MODEL_NP; m++)
        {
            fprintf(corrFile,"%+.3f ",corr[n][m]);
        }
        fprintf(corrFile,"\n");
    }
    fclose(corrFile);

    
    /* *************************************************************** */
    /*           Check sources against previous catalog                */
    /* *************************************************************** */
    if(flags->catalog)
    {
        printf("get source ancestry\n");

        //parse file
        FILE *old_catalog_file = fopen(flags->catalogFile,"r");
        
        //count number of sources
        
        struct Source *old_catalog_entry = NULL;
        old_catalog_entry = malloc(sizeof(struct Source));
        alloc_source(old_catalog_entry, data->N, data->Nchannel);
        
        struct Source *new_catalog_entry = NULL;
        new_catalog_entry = malloc(sizeof(struct Source));
        alloc_source(new_catalog_entry, data->N, data->Nchannel);
        
        
        int Nsource = 0;
        while(!feof(old_catalog_file))
        {
            scan_source_params(data, old_catalog_entry, old_catalog_file);
            Nsource++;
        }
        Nsource--;
        rewind(old_catalog_file);
        
        sprintf(filename,"%s/history.dat",outdir);
        FILE *historyFile = fopen(filename,"w");
        for(int i=0; i<Nsource; i++)
        {
            scan_source_params(data_old, old_catalog_entry, old_catalog_file);
            
            //find where the source fits in the measurement band
            ucb_alignment(orbit, data_old, old_catalog_entry);
            
            //find central bin of catalog event for current data
            double q_old_catalog_entry = old_catalog_entry->f0 * data_old->T;
            
            //check that source lives in current data segment (allow sources in padded region)
            if(q_old_catalog_entry < data_old->qmin || q_old_catalog_entry > data_old->qmax) continue;
            
            //calculate waveform model of sample at Tcatalog
            ucb_waveform(orbit, data->format, data_old->T, data->t0, old_catalog_entry->params, UCB_MODEL_NP, old_catalog_entry->tdi->X, old_catalog_entry->tdi->Y,old_catalog_entry->tdi->Z,old_catalog_entry->tdi->A, old_catalog_entry->tdi->E, old_catalog_entry->BW, data->Nchannel);
            
            //check against new catalog
            for(int d=0; d<detections; d++)
            {
                int n = detection_index[d];
                entry = catalog->entry[n];
                
                double q_new_catalog_entry = entry->source[entry->i]->f0 * data_old->T;
                
                //check that the sources are close enough to bother looking
                if(fabs(q_new_catalog_entry - q_old_catalog_entry) > dqmax) continue;
                
                //copy_source(entry->source[entry->i], new_catalog_entry);
                
                /* parameter-only copy */
                memcpy(new_catalog_entry->params, entry->source[entry->i]->params, UCB_MODEL_NP*sizeof(double));
                
                //override q parameter
                new_catalog_entry->params[0] = entry->source[entry->i]->f0 * data_old->T;
                new_catalog_entry->params[7] = entry->source[entry->i]->dfdt * data_old->T* data_old->T;

                //re-align where the source fits in the (old) measurement band
                ucb_alignment(orbit, data_old, new_catalog_entry);
                
                //calculate waveform of entry at Tcatalog
                ucb_waveform(orbit, data->format, data_old->T, data->t0, new_catalog_entry->params, UCB_MODEL_NP, new_catalog_entry->tdi->X, new_catalog_entry->tdi->Y, new_catalog_entry->tdi->Z, new_catalog_entry->tdi->A, new_catalog_entry->tdi->E, new_catalog_entry->BW, data->Nchannel);
                
                
                Match = waveform_match(old_catalog_entry,new_catalog_entry,noise);
                
                if(Match>0.5)
                {
                    sprintf(entry->parent,"LDC%010li",(long)(old_catalog_entry->f0*1e10));
                    fprintf(historyFile,"%s %s\n",entry->parent,entry->name);
                }
                
            }
        }
        free_source(old_catalog_entry);
        free_source(new_catalog_entry);
        fclose(historyFile);
    }
    
    
    
    /* *************************************************************** */
    /*           Save source detection parameters to file              */
    /* *************************************************************** */
    
    
    omp_set_num_threads(detections);
    
    printf("get parameter chains\n");
#pragma omp parallel for default(shared) private(entry,filename)
    for(int d=0; d<detections; d++)
    {
        FILE *out;
        
        //print detection posterior samples
        int n = detection_index[d];
        entry = catalog->entry[n];
        
        sprintf(filename, "%s/%s_chain.dat", outdir,entry->name);
        out = fopen( filename, "w");
        
        //add parameters to file
        for(int k=0; k<entry->Nchain; k++)
        {
            print_source_params(data,entry->source[k],out);
            
            source_waveform_wrapper(entry->source[k], data, orbit);

            fprintf(out,"%lg %lg %lg\n",snr(entry->source[k],noise),entry->match[k],entry->distance[k]);
            
            free_tdi(entry->source[k]->tdi);
        }
        fclose(out);
    }
    
    printf("get waveform reconstructions\n");
#pragma omp parallel for default(shared) private(entry,filename)
    for(int d=0; d<detections; d++)
    {

        int n = detection_index[d];
        entry = catalog->entry[n];

        //create and print individual source waveform reconstructions
        double ***hrec = malloc(data->N * sizeof(double **));
        for(int j=0; j<data->N; j++)
        {
            hrec[j] = malloc(data->Nchannel * sizeof(double *));
            
            //allocate and zero
            for(int k=0; k<data->Nchannel; k++)
            hrec[j][k] = calloc(entry->Nchain , sizeof(double));
        }
        
        //insert waveform power
        for(int i=0; i<entry->Nchain; i++)
        {
            
            source_waveform_wrapper(entry->source[i], data, orbit);

            for(int j=0; j<entry->source[i]->BW; j++)
            {
                
                int k = j+entry->source[i]->imin;
                
                if(k>-1 && k < data->N)
                {
                    int j_re = 2*j;
                    int j_im = j_re+1;
                    
                    hrec[k][0][i] = entry->source[i]->tdi->A[j_re]*entry->source[i]->tdi->A[j_re]+entry->source[i]->tdi->A[j_im]*entry->source[i]->tdi->A[j_im];
                    hrec[k][1][i] = entry->source[i]->tdi->E[j_re]*entry->source[i]->tdi->E[j_re]+entry->source[i]->tdi->E[j_im]*entry->source[i]->tdi->E[j_im];
                }
            }
            
            free_tdi(entry->source[i]->tdi);
        }
        
        //sort reconstructed power in each frequency bin and get median, CIs
        for(int j=0; j<data->N; j++)
        {
            for(int k=0; k<data->Nchannel; k++)
            {
                double_sort(hrec[j][k],entry->Nchain);
            }
        }
        double A_med,A_lo_50,A_hi_50,A_lo_90,A_hi_90;
        double E_med,E_lo_50,E_hi_50,E_lo_90,E_hi_90;
        
        sprintf(filename, "%s/%s_power_reconstruction.dat", outdir,entry->name);
        FILE *out = fopen( filename, "w");
        
        
        for(int j=0; j<data->N; j++)
        {
            double f = (double)(j+data->qmin)/data->T;
            
            A_med   = get_quantile_from_sorted_data (hrec[j][0], entry->Nchain, 0.50);
            A_lo_50 = get_quantile_from_sorted_data (hrec[j][0], entry->Nchain, 0.25);
            A_hi_50 = get_quantile_from_sorted_data (hrec[j][0], entry->Nchain, 0.75);
            A_lo_90 = get_quantile_from_sorted_data (hrec[j][0], entry->Nchain, 0.05);
            A_hi_90 = get_quantile_from_sorted_data (hrec[j][0], entry->Nchain, 0.95);
            
            E_med   = get_quantile_from_sorted_data (hrec[j][1], entry->Nchain, 0.50);
            E_lo_50 = get_quantile_from_sorted_data (hrec[j][1], entry->Nchain, 0.25);
            E_hi_50 = get_quantile_from_sorted_data (hrec[j][1], entry->Nchain, 0.75);
            E_lo_90 = get_quantile_from_sorted_data (hrec[j][1], entry->Nchain, 0.05);
            E_hi_90 = get_quantile_from_sorted_data (hrec[j][1], entry->Nchain, 0.95);
            
            fprintf(out,"%.12g ",f);
            fprintf(out,"%lg ",A_med);
            fprintf(out,"%lg ",A_lo_50);
            fprintf(out,"%lg ",A_hi_50);
            fprintf(out,"%lg ",A_lo_90);
            fprintf(out,"%lg ",A_hi_90);
            fprintf(out,"%lg ",E_med);
            fprintf(out,"%lg ",E_lo_50);
            fprintf(out,"%lg ",E_hi_50);
            fprintf(out,"%lg ",E_lo_90);
            fprintf(out,"%lg ",E_hi_90);
            fprintf(out,"\n");
        }
        
        fclose(out);

        for(int i=0; i<data->N; i++)
        {
            for(int j=0; j<data->Nchannel; j++)
            {
                free(hrec[i][j]);
            }
            free(hrec[i]);
        }
        free(hrec);
    }
    
    printf("get gaussian mixture model\n");

#pragma omp parallel for default(shared) private(entry)
    for(int d=0; d<detections; d++)
    {
        int n = detection_index[d];
        entry = catalog->entry[n];

        /* Gaussian mixture model fit to posterior */
        //fprintf(stdout,"\nGaussian Mixture Model fit:\n");
        unsigned int r = 190521;

        int counter;
        int CMAX = 10;
        int NMODE_start = (int)NMODE;
        double BIC;
        
        counter = 0;
        while(gaussian_mixture_model_wrapper(model->prior, flags, entry, outdir, NMODE_start, NTHIN, &r, &BIC))
        {
            counter++;
            if(counter>CMAX)
            {
                fprintf(stderr,"WARNING:\n");
                fprintf(stderr,"Gaussian Mixture Model failed to converge for source %s\n",entry->name);
                NMODE_start/=2;
                counter = 0;
            }
        }
        
    }//end loop over catalog entries
    
    
    return 0;
}


