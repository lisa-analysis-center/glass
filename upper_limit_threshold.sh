
# don't exit if a command errors
set +e

NSTEPS=1000
DURATION=$((7680*338*24))
NCH=4

RUNDIR=/mnt/storage/research/wdmruns-coarse5/noise-conf-ln-4
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
        --sgwb-template 3 \
	--sgwb-inj -4 \
        --coarse-Q 169 \
        --rundir $RUNDIR \
        --sim-noise
RUNDIR=/mnt/storage/research/wdmruns-coarse5/noise-conf-ln-3
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
        --sgwb-template 3 \
	--sgwb-inj -3 \
        --coarse-Q 169 \
        --rundir $RUNDIR \
        --sim-noise
RUNDIR=/mnt/storage/research/wdmruns-coarse5/noise-conf-ln-25
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
        --sgwb-template 3 \
	--sgwb-inj -2.5 \
        --coarse-Q 169 \
        --rundir $RUNDIR \
        --sim-noise
