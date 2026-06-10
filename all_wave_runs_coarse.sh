
# don't exit if a command errors
set +e

NSTEPS=1000
DURATION=$((7680*338*24))
NCH=6

RUNDIR=wdmruns-coarse3/stat-noise-conf-pt
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps $((10*$NSTEPS)) \
        --chains $NCH \
        --threads $NCH \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
        --duration $DURATION \
        --stationary \
        --conf-noise \
        --sgwb-template 2 \
        --coarse-Q 1 \
        --rundir $RUNDIR \
        --sim-noise

RUNDIR=wdmruns-coarse3/noise-conf-pt-169
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps $((10*$NSTEPS)) \
        --chains $NCH \
        --threads $NCH \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
        --duration $DURATION \
        --conf-noise \
        --sgwb-template 2 \
        --coarse-Q 169 \
        --rundir $RUNDIR \
        --sim-noise

#RUNDIR=wdmruns-coarse3/stat-noise-conf-pl
#mkdir -p $RUNDIR
#nice -n 15 build/apps/src/noise_wavelet_mcmc \
#        --steps $((10*$NSTEPS)) \
#        --chains $NCH \
#        --threads $NCH \
#        --cheat \
#        --fmin 5e-4 \
#        --fmax 8e-3 \
#        --duration $DURATION \
#        --stationary \
#        --conf-noise \
#        --sgwb-template 0 \
#        --coarse-Q 1 \
#        --rundir $RUNDIR \
#        --sim-noise

#RUNDIR=wdmruns-coarse3/noise-conf-pl-169
#mkdir -p $RUNDIR
#nice -n 15 build/apps/src/noise_wavelet_mcmc \
#        --steps $((10*$NSTEPS)) \
#        --chains $NCH \
#        --threads $NCH \
#        --cheat \
#        --fmin 5e-4 \
#        --fmax 8e-3 \
#        --duration $DURATION \
#        --conf-noise \
#        --sgwb-template 0 \
#        --coarse-Q 169 \
#        --rundir $RUNDIR \
#        --sim-noise

#RUNDIR=wdmruns-coarse3/stat-noise-conf-ln-detection
#mkdir -p $RUNDIR
#nice -n 15 build/apps/src/noise_wavelet_mcmc \
#        --steps $((10*$NSTEPS)) \
#        --chains $NCH \
#        --threads $NCH \
#        --cheat \
#        --fmin 5e-4 \
#        --fmax 8e-3 \
#        --duration $DURATION \
#        --stationary \
#        --conf-noise \
#        --sgwb-template 1 \
#        --coarse-Q 1 \
#        --rundir $RUNDIR \
#        --sim-noise
##        --stationary-conf \
#
#RUNDIR=wdmruns-coarse3/statws-noise-conf-ln-detection
#mkdir -p $RUNDIR
#nice -n 15 build/apps/src/noise_wavelet_mcmc \
#        --steps $((10*$NSTEPS)) \
#        --chains $NCH \
#        --threads $NCH \
#        --cheat \
#        --fmin 5e-4 \
#        --fmax 8e-3 \
#        --duration $DURATION \
#        --stationary \
#        --conf-noise \
#        --sgwb-template 1 \
#        --coarse-Q 1 \
#        --ws-approx \
#        --rundir $RUNDIR \
#        --sim-noise
##        --stationary-conf \
#
#RUNDIR=wdmruns-coarse3/noise-conf-ln-169-detection
#mkdir -p $RUNDIR
#nice -n 15 build/apps/src/noise_wavelet_mcmc \
#        --steps $((10*$NSTEPS)) \
#        --chains $NCH \
#        --threads $NCH \
#        --cheat \
#        --fmin 5e-4 \
#        --fmax 8e-3 \
#        --duration $DURATION \
#        --conf-noise \
#        --sgwb-template 1 \
#        --coarse-Q 169 \
#        --rundir $RUNDIR \
#        --sim-noise

#RUNDIR=wdmruns-coarse3/noise-conf-ln-ws169-avg-detection
#mkdir -p $RUNDIR
#nice -n 15 build/apps/src/noise_wavelet_mcmc \
#        --steps $((10*$NSTEPS)) \
#        --chains $NCH \
#        --threads $NCH \
#        --cheat \
#        --fmin 5e-4 \
#        --fmax 8e-3 \
#        --duration $DURATION \
#        --conf-noise \
#        --sgwb-template 1 \
#        --coarse-Q 169 \
#        --ws-approx \
#        --rundir $RUNDIR \
#        --sim-noise

#RUNDIR=wdmruns-coarse3/noise-conf-ln-1-detection
#mkdir -p $RUNDIR
#nice -n 15 build/apps/src/noise_wavelet_mcmc \
#        --steps $NSTEPS \
#        --chains $NCH \
#        --threads $NCH \
#        --cheat \
#        --fmin 5e-4 \
#        --fmax 8e-3 \
#        --duration $DURATION \
#        --conf-noise \
#        --sgwb-template 1 \
#        --coarse-Q 1 \
#        --rundir $RUNDIR \
#        --sim-noise
