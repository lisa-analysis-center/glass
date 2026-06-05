
# don't exit if a command errors
set +e

NSTEPS=100000
DURATION=$((7680*338*24))
NCH=4

#RUNDIR=/mnt/storage/research/wdmruns-coarse4/noise-conf-ln-169-detection
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
#
#RUNDIR=/mnt/storage/research/wdmruns-coarse4/noise-conf-ln-ws169-detection
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

RUNDIR=/mnt/storage/research/wdmruns-coarse4/stat-noise-conf-pt
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps $NSTEPS \
        --chains $NCH \
        --threads $NCH \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
	--stationary
        --duration $DURATION \
        --conf-noise \
        --sgwb-template 2 \
        --coarse-Q 1 \
        --rundir $RUNDIR \
        --sim-noise

RUNDIR=/mnt/storage/research/wdmruns-coarse4/noise-conf-pt-169
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps $NSTEPS \
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
