
# don't exit if a command errors
set +e

NSTEPS=2000
DURATION=$((7680*338*24))
NCH=6

# run_mcmc RUNDIR [extra args...]
run_mcmc() {
        local rundir=$1
        shift
        # skips the run if it already finished (final noise model present)
        if [ -f "$rundir/data/final_full_noise_model.dat" ]; then
                echo "Skipping $rundir (already finished)"
                return 0
        fi
        mkdir -p "$rundir"
        nice -n 15 build/apps/src/noise_wavelet_mcmc \
                --steps $((10*$NSTEPS)) \
                --chains $NCH \
                --threads $NCH \
                --cheat \
                --fmin 5e-4 \
                --fmax 8e-3 \
                --duration $DURATION \
                --conf-noise \
                --rundir "$rundir" \
                --sim-noise \
                "$@"
}

# phase transition upper limits
run_mcmc ./wdmruns-coarse6/stat-noise-conf-pt --sgwb-template 2 --coarse-Q 1 --stationary --sgwb-inj  1.0,1.0,-18.0,-2.5
run_mcmc ./wdmruns-coarse6/noise-conf-pt-169 --sgwb-template 2 --coarse-Q 169 --sgwb-inj  1.0,1.0,-18.0,-2.5

# powerlaw upper limits
run_mcmc ./wdmruns-coarse6/stat-noise-conf-pl --sgwb-template 0 --coarse-Q 1 --stationary --sgwb-inj  -20.0,0.6667
run_mcmc ./wdmruns-coarse6/noise-conf-pl-169 --sgwb-template 0 --coarse-Q 169 --sgwb-inj  -20.0,0.6667

# lognormal upper limits
run_mcmc ./wdmruns-coarse6/stat-noise-conf-ln --sgwb-template 1 --coarse-Q 1 --stationary --sgwb-inj  -4.0,-2.69,0.0
run_mcmc ./wdmruns-coarse6/noise-conf-ln-169 --sgwb-template 1 --coarse-Q 169 --sgwb-inj  -4.0,-2.69,0.0

# powerlaw detections
run_mcmc ./wdmruns-coarse6/stat-noise-conf-pl-detection --sgwb-template 0 --coarse-Q 1 --stationary --sgwb-inj  -9.5,0.6667
run_mcmc ./wdmruns-coarse6/noise-conf-pl-169-detection --sgwb-template 0 --coarse-Q 169 --sgwb-inj  -9.5,0.6667

# lognormal detections
run_mcmc ./wdmruns-coarse6/stat-noise-conf-ln-detection --sgwb-template 1 --coarse-Q 1 --stationary --sgwb-inj  -3.0,-2.69,0.0
run_mcmc ./wdmruns-coarse6/noise-conf-ln-169-detection --sgwb-template 1 --coarse-Q 169 --sgwb-inj  -3.0,-2.69,0.0
