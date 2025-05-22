/*
 * Copyright 2021 Tyson B. Littenberg (MSFC-ST12)
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
 @file glass_noise_wrapper.h
 \brief Wrapper functions to noise sampler.

 Functions for handling the memory, data sharing, and sampling of the NOISE block for the global fit analysis.
 */

#ifndef NoiseWrapper_h
#define NoiseWrapper_h

struct NoiseData
{
    int mcmc_step;
    int status;
    
    int procID; //!<MPI process identifier
    int nProc;  //!<Number of processing segments
    int procID_min; //!<lowest rank MPI process for Noise block
    int procID_max; //!<highest rank MPI process for Noise block

    struct Flags *flags;
    struct Orbit *orbit;
    struct Chain *chain;
    struct Data  *data;
    struct Noise **psd;
    struct InstrumentModel **inst_model;
    struct InstrumentModel **inst_trial;
    struct ForegroundModel **conf_model;
    struct ForegroundModel **conf_trial;

    double cpu_time; //!<CPU time for single update
};

void alloc_noise_data(struct NoiseData *noise_data, struct UCBData *ucb_data, int procID, int nProc);
void setup_noise_data(struct NoiseData *noise_data, struct UCBData *ucb_data, struct VGBData *vbmcmc_data, struct MBHData *mbh_data, struct TDI *tdi_full, int procID);
void select_noise_segment(struct Noise *psd_full, struct Data *data, struct Chain *chain, struct Model **model);

void initialize_noise_sampler(struct NoiseData *noise_data);
void initialize_noise_state(struct NoiseData *noise_data);
void resume_noise_state(struct NoiseData *noise_data);

int update_noise_sampler(struct NoiseData *noise_data);

/* bad name because of naming conflict in GalacticBinaryIO.c */
void print_noise_state(struct NoiseData *noise_data, FILE *fptr, int counter);

#endif /* NoiseWrapper_h */
