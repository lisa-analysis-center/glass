
# don't exit if a command errors
set +e

NSTEPS=1000
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
                --sgwb-template 3 \
                --coarse-Q 169 \
                --rundir "$rundir" \
                --sim-noise \
                "$@"
}

run_mcmc ./wdmruns-coarse5/noise-conf-ln-4 --sgwb-inj -4
run_mcmc ./wdmruns-coarse5/stat-noise-conf-ln-4 --sgwb-inj -4 --stationary
run_mcmc ./wdmruns-coarse5/noise-conf-ln-3 --sgwb-inj -3
run_mcmc ./wdmruns-coarse5/stat-noise-conf-ln-3 --sgwb-inj -3 --stationary
run_mcmc ./wdmruns-coarse5/noise-conf-ln-35 --sgwb-inj -3.5
run_mcmc ./wdmruns-coarse5/stat-noise-conf-ln-35 --sgwb-inj -3.5 --stationary
