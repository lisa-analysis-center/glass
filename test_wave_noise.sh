#! /usr/bin/env bash

#rm -rf ./chains
#rm -rf ./data
#OMP_NUM_THREADS=1 nice -n 15 build/apps/src/noise_wavelet_mcmc --sim-noise --conf-noise --sgwb-template 0 --steps 10000 --chains 1 --cheat --fmin 1e-4 --fmax 8e-3 --duration $((60*60*24*30*12))
#OMP_NUM_THREADS=6 build/apps/src/noise_mcmc --sim-noise --sgwb-template 0 --steps 2000 --chains 6 --cheat --fmin 1e-4 --fmax 8e-3 --duration 1966080

#OMP_NUM_THREADS=1 lldb -- build/apps/src/noise_wavelet_mcmc --sim-noise --conf-noise --sgwb-template 0 --steps 10 --chains 1 --cheat --fmin 5e-4 --fmax 8e-3 --duration $((60*60*24*30*12)) --threads 1
#OMP_NUM_THREADS=1 build/apps/src/noise_wavelet_mcmc --sim-noise --conf-noise --sgwb-template 0 --steps 10000 --chains 1 --cheat --fmin 5e-4 --fmax 8e-3 --duration $((60*60*24*30*1)) --threads 1
#OMP_NUM_THREADS=1 build/apps/src/noise_wavelet_mcmc --sim-noise --conf-noise --sgwb-template 0 --steps 500 --chains 1 --cheat --fmin 1e-3 --fmax 1.1e-3 --duration $((60*60*24*30*12)) --threads 1
#build/apps/src/noise_wavelet_mcmc \
OMP_NUM_THREADS=1 lldb -- build/apps/src/noise_wavelet_mcmc \
        --steps 1 \
        --chains 1 \
        --threads 1 \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
        --duration $((7680*338*24)) \
        --ws-approx \
        --coarse-Q 169 \
        --sgwb-template 0 \
        --conf-noise \
        --sim-noise \
        --inst-drift
#        --threads 1 \
#        --timeseries \
#        --data ./simulated_data_timeseries.dat
