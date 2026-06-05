/*
 * Copyright 2023 Tyson B. Littenberg & Neil J. Cornish
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
 @file glass_noise_model.h
 \brief Functions and structures defining noise models.

 */
#ifndef noise_model_h
#define noise_model_h

#define MIN_SPLINE_STENCIL 5
#define MIN_SPLINE_SPACING 256.

/*!
 * \brief Data strictire for spline model used in noise model.
 *
 * The SplineModel structure stores metadata and state of spline model,
 * including size of model, current state of spline control points
 * and likelihood, and interpolated model.
 */
struct SplineModel
{
    int Nmin; //!< Minimum number of spline control points
    int Nmax; //!< Maximum number of spline control points
    int Nchannel; //!< Number of TDI channels
    double logL; //!< Log likelhood of spline model
    struct Noise *spline; //!< Spline control points
    struct Noise *psd; //!< Reconstructed noise model
};

struct InstrumentModel
{
    int Nlink;         //!< number of interferometer links (6 for 3-spacecraft constellation)
    double logL;       //!< log Likelihood of model
    double *soms;      //!< optical metrology system noise parameters
    double *sacc;      //!< acceleration noise parameters
    struct Noise *psd; //!< power and cross spectral densities

    //!< Reusable scratch model on the FFT-bin grid; populated lazily by
    //!< generate_instrument_noise_model_wavelet and reused across calls to
    //!< avoid ~NFFT-sized malloc/free churn in the MCMC inner loop. NULL
    //!< until the first wavelet-domain generate; freed by free_instrument_model.
    struct InstrumentModel *grid_cache;

    /** @name Link level noise parameters */
     ///@{
    double sacc12;
    double sacc13;
    double sacc21;
    double sacc23;
    double sacc31;
    double sacc32;
    double soms12;
    double soms13;
    double soms21;
    double soms23;
    double soms31;
    double soms32;
    ///@}
};

struct ForegroundModel
{
    int Nparams;       //!< number of foreground parameters (7)
    double Tobs;       //!< observation time (in seconds)
    double *sgal;      //!< galactic foreground parameters
    double logL;       //!< log Likelihood of model
    struct Noise *psd; //!< power and cross spectral densities

    //!< Reusable scratch model on the FFT-bin grid (see InstrumentModel::grid_cache).
    struct ForegroundModel *grid_cache;

    /** @name galactic foreground parameters from Digman & Cornish 10.3847/1538-4357/ac9139*/
     ///@{
    double Amp;   //!< overall amplitude  (~2.13e-37)
    double f1;
    double alpha; //!< spectral index of low frequency part (~1.58)
    double fk;
    double f2;    //!< some other frequency parameter? (~0.53 mHz)
     ///@}

    struct GalaxyModulation *modulation; //!< time-series of modulation
};

typedef enum {
    SGWB_TEMPLATE_POWERLAW,
    SGWB_TEMPLATE_LOGNORMAL,
    SGWB_TEMPLATE_PHASE_TRANSITION,
    SGWB_TEMPLATE_COUNT // leave this here, counts length of enum
} SGWB_t;
// note that if this array ever becomes extremely large, maybe swap to extern
static const char* SGWB_TEMPLATE_NAMES[SGWB_TEMPLATE_COUNT] = {"powerlaw", "lognormal", "phase_transition"};
static const int SGWB_TEMPLATE_NPARAMS[SGWB_TEMPLATE_COUNT] = {2, 3, 4};
_Static_assert(sizeof(SGWB_TEMPLATE_NAMES)/sizeof(SGWB_TEMPLATE_NAMES[0]) == SGWB_TEMPLATE_COUNT,
        "Did you add an SGWB template but not its name?");
_Static_assert(sizeof(SGWB_TEMPLATE_NPARAMS)/sizeof(SGWB_TEMPLATE_NPARAMS[0]) == SGWB_TEMPLATE_COUNT,
        "Did you add an SGWB template but not its parameter count?");

struct SGWBModel
{
    int Nparams;       //!< number of SGWB parameters
    double Tobs;       //!< observation time (in seconds)
    double *params;    //!< SGWB parameters
    double logL;       //!< log Likelihood of model
    struct Noise *psd; //!< power and cross spectral densities
    struct SGWBResponse *R; //!< response of the instrument to the isotropic SGWB
    SGWB_t SGWB_type;  //!< which SGWB template is used

    //!< Reusable scratch model on the FFT-bin grid (see InstrumentModel::grid_cache).
    struct SGWBModel *grid_cache;
};

/**
 \brief Converts physical noise parameters to array expected by InstrumentModel
 */
void map_noise_params_to_array(struct InstrumentModel *model);

/**
 \brief Converts array expected by InstrumentModel to
 physical noise parameters
 */
void map_array_to_noise_params(struct InstrumentModel *model);

/**
 \brief Converts phenomenological foreground parameters to array expected by ForegroundModel
 */
void map_foreground_params_to_array(struct ForegroundModel *model);

/**
 \brief Converts array expected by ForegroundModel to
 phenomenological foreground parameters
 */
void map_array_to_foreground_params(struct ForegroundModel *model);

/**
 \brief Allocates spline model structure and contents.
 */
void alloc_spline_model(struct SplineModel *model, int Ndata, int Nlayer, int Nchannel, int Nspline);

/**
 \brief Allocates instrument model structure and contents.
 */
void alloc_instrument_model(struct InstrumentModel *model, int Ndata, int Nlayer, int Nchannel);

/**
 \brief Allocates galactic foreground model structure and contents.
 */
void alloc_foreground_model(struct ForegroundModel *model, int Ndata, int Nlayer, int Nchannel);

/**
 \brief Allocates SGWB model structure and contents.
 */
void alloc_sgwb_model(struct SGWBModel *model, int Ndata, int Nlayer, int Nchannel, SGWB_t SGWB_type);

/**
 \brief Allocates SGWB response.
 */
void alloc_pop_sgwb_response(struct SGWBResponse* sgwbr, char* fname);

/**
 \brief Free allocated spline model.
 */
void free_spline_model(struct SplineModel *model);

/**
 \brief Free allocated spline model.
 */
void free_instrument_model(struct InstrumentModel *model);

/**
 \brief Free allocated foreground model.
 */
void free_foreground_model(struct ForegroundModel *model);

/**
 \brief Free allocated SGWB model.
 */
void free_sgwb_model(struct SGWBModel *model);

/**
 \brief Deep copy of SplineModel structure from `origin` into `copy`
 */
void copy_spline_model(struct SplineModel *origin, struct SplineModel *copy);

/**
 \brief Deep copy of InstrumentModel structure from `origin` into `copy`
 */
void copy_instrument_model(struct InstrumentModel *origin, struct InstrumentModel *copy);

/**
 \brief Deep copy of ForegroundModel structure from `origin` into `copy`
 */
void copy_foreground_model(struct ForegroundModel *origin, struct ForegroundModel *copy);

/**
 \brief Deep copy of SGWBModel structure from `origin` into `copy`
 */
void copy_sgwb_model(struct SGWBModel *origin, struct SGWBModel *copy);

/**
 \brief Wrapper to `CubicSplineGSL` functions for generating PSD model based on current state of `model`
 */
void generate_spline_noise_model(struct SplineModel *model);

/**
 \brief Compute instrument model contribution to noise covariance matrix based on current state of `model`
 */
void generate_instrument_noise_model(struct Orbit *orbit, struct InstrumentModel *model);
void generate_instrument_noise_model_wavelet(struct Wavelets *wdm, struct Orbit *orbit, struct InstrumentModel *model);

/**
 \brief Compute galactic foreground contribution to covariance matrix based on current state of `model`
 */
void generate_galactic_foreground_model(struct ForegroundModel *model);
void generate_galactic_foreground_model_wavelet(struct Wavelets *wdm, struct ForegroundModel *model);

/**
 \brief Compute SGWB contribution to covariance matrix based on current state of `model`
 */
void generate_sgwb_model(struct SGWBModel *model);
void generate_sgwb_model_wavelet(struct Wavelets *wdm, struct SGWBModel *model);

/**
 \brief Add components to `full` noise covariance matrix `C`
 */
void generate_full_covariance_matrix(struct Noise *full, struct Noise *component, int Nchannel);
void generate_full_dynamic_covariance_matrix(struct Wavelets *wdm, struct InstrumentModel *inst, struct ForegroundModel *conf, struct SGWBModel *sgwb, struct Noise *full);
void generate_full_stationary_covariance_matrix(struct Wavelets *wdm, struct InstrumentModel *inst, struct ForegroundModel *conf, struct SGWBModel *sgwb, struct Noise *full);

/**
 \brief Coarse-grained covariance generators for the WDM time-axis Welch/Bartlett method.

 The coarse `Noise` struct must be allocated with `Nlayer*Ncoarse` pixels, where
 `Ncoarse = wdm->NT / Q`. Pixel index is `k = q + (j-jmin)*Ncoarse`.
 Foreground modulation (dynamic only) is sampled at the coarse-cell midpoint.
 */
void generate_full_dynamic_covariance_matrix_coarse(struct Wavelets *wdm, int Q, struct InstrumentModel *inst, struct ForegroundModel *conf, struct SGWBModel *sgwb, struct Noise *coarse);
void generate_full_stationary_covariance_matrix_coarse(struct Wavelets *wdm, int Q, struct InstrumentModel *inst, struct ForegroundModel *conf, struct SGWBModel *sgwb, struct Noise *coarse);

/**
 \brief Sufficient statistics for the Welch/Bartlett-like coarse-grained likelihood.

 The six P^{(ij)}_{mq} arrays are the per-coarse-cell sample covariances
 (1/Q) sum_{i in cell q} w_{a,mi} w_{b,mi}. Each array has length
 `Nlayer * Ncoarse` (Ncoarse = NT/Q), indexed by `k = q + j*Ncoarse`.

 Lifetime contract: populated once via `coarse_grain_wavelet_data` after data
 injection / read and treated as read-only thereafter. Const-qualified pointer
 fields let multiple MCMC threads share a single instance safely. The owner
 (allocator) must keep non-const handles internally; downstream code passes
 `const struct CoarseStats *`.
 */
struct CoarseStats
{
    const double *Pxx;
    const double *Pyy;
    const double *Pzz;
    const double *Pxy;
    const double *Pxz;
    const double *Pyz;
    /**
     Per-coarse-cell effective number of degrees of freedom (length
     `Nlayer * Ncoarse`, same `k = q + j*Ncoarse` indexing as the P arrays).
     Default-filled with `Q` by `alloc_coarse_stats` so the likelihood is
     bit-for-bit unchanged unless `precompute_coarse_Qeff` overwrites it.
     When the `--ws-approx` Welch-Satterthwaite correction is active these hold
     Qeff <= Q, frozen at the injected covariance. See `precompute_coarse_Qeff`.
     */
    const double *Qeff;
    int Q;
    int Ncoarse;
    int Nlayer;
    /** Coarse-cell index range [qlo, qhi) that the likelihood actually sums over.
        Defaults to [0, Ncoarse) (all cells). Set narrower to crop edge time
        pixels (e.g. the WDM edge-effect buffer): the P/Qeff arrays still hold all
        Ncoarse cells, but cells outside [qlo, qhi) are excluded from the
        likelihood and its normalization. Cropping is quantized to whole coarse
        cells. */
    int qlo;
    int qhi;
};

void alloc_coarse_stats(struct CoarseStats *s, int Nlayer, int Q, int NT);
void free_coarse_stats(struct CoarseStats *s);

/**
 \brief Precompute and freeze a Welch-Satterthwaite effective dof (Qeff) per
 coarse cell from the injected covariance in `data->noise`.

 For each coarse cell and each diagonal channel a, with per-fine-pixel injected
 variance C_a(t_n), forms Qeff_a = (sum_n C_a)^2 / sum_n C_a^2 (the moment-matched
 effective dof of the heterogeneous-variance sum) and stores the channel-average
 (Qeff_xx + Qeff_yy + Qeff_zz)/3 into `stats->Qeff`. Requires `data->noise` to
 hold the full-resolution (Q=1) injected covariance, i.e. an injection run.
 */
void precompute_coarse_Qeff(struct Data *data, struct CoarseStats *stats);

/**
 \brief Build the six channel-pair sufficient statistics P^{(ij)}_{mq} from the
 wavelet-domain data into `stats` (which must already be allocated for this Q).
 */
void coarse_grain_wavelet_data(struct Data *data, struct CoarseStats *stats);

/**
\brief Compute spline model only where interpolant changes
 
 Interpolates spline points only in the vicinity of `new_knot`
 */
void update_spline_noise_model(struct SplineModel *model, int new_knot, int min_knot, int max_knot);

/**
 \brief Log likelihood for noise model.
 
 @param data Data structure
 @param model SplineModel structure containing current state
 @return \f$  \ln p({\rm data}|{\rm spline}) \f$
 */
double noise_log_likelihood(struct Data *data, struct Noise *noise);
double my_noise_log_likelihood(struct Data *data, struct Noise *noise);
double noise_log_likelihood_wavelet(struct Data *data, struct Noise *noise);
double my_noise_log_likelihood_wavelet(struct Data *data, struct Noise *noise);

/**
 \brief Welch/Bartlett-like coarse-grained wavelet-domain log-likelihood.

 Uses the multivariate Gaussian sufficient-statistics form:
   logL = -Q/2 sum_{m,q} [ tr(invC_{mq} * P_{mq}) + log det C_{mq} ] - (3 N_full / 2) log(2 pi)
 where N_full = Nlayer * NT is the original full-grid pixel count.

 Reduces algebraically to `my_noise_log_likelihood_wavelet` when the per-pixel
 covariance is constant in each window of Q time pixels (e.g. stationary noise).
 */
double my_noise_log_likelihood_wavelet_coarse(struct Data *data, struct Noise *coarse_noise, const struct CoarseStats *stats);

/**
 \brief Change in log likelihood for noise model.
 
 @param data Data structure
 @param model_x Current SplineModel structure containing current state
 @param model_y Current SplineModel structure containing current state
 @return \f$  \ln p({\rm data}|{\rm spline}_y)-\ln p({\rm data}|{\rm spline}_x) \f$
 */
double noise_delta_log_likelihood(struct Data *data, struct SplineModel *model_x, struct SplineModel *model_y, double fmin, double fmax, int ic);

/**
 \brief Set initial state of spline `model`
 */
void initialize_spline_model(struct Orbit *orbit, struct Data *data, struct SplineModel *model, int Nspline);

/**
 \brief Set initial state of instrument noise `model`
 */
void initialize_instrument_model(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model);
void initialize_instrument_model_wavelet(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model);
/**
 \brief Set initial state of galactic foreground `model`
 */
void initialize_foreground_model(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model);
void initialize_foreground_model_wavelet(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model);

/**
 \brief Set initial state of sgwb `model`
 */
void initialize_sgwb_model(struct Orbit *orbit, struct Data *data, struct SGWBModel *model, SGWB_t SGWB_type);
void initialize_sgwb_model_wavelet(struct Orbit *orbit, struct Data *data, struct SGWBModel *model, SGWB_t SGWB_type);


void GetDynamicNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags);
void GetStationaryNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Noise *noise);

/**
 \brief functional form for a powerlaw SGWB
 */
double sgwb_powerlaw(double f, const double* params);

/**
 \brief functional form for a (wide) lognormal scalar-induced SGWB

 Pi & Sasaki, JCAP 2020 (arXiv:2005.12306) eq. (3.29). Parameters are
 \f$\log_{10} A\f$, \f$\log_{10} f_*\f$, \f$\log_{10} \Delta\f$.
 */
double sgwb_lognormal(double f, const double* params);

/**
 \brief functional form for a first-order phase transition SGWB

 Broken power-law shape \f$M(f/f_p)\f$ with peak amplitude \f$A_p\f$ and peak
 frequency \f$f_p\f$. Parameters are \f$r_b\f$, \f$b\f$,
 \f$\log_{10} A_p\f$, \f$\log_{10} f_p\f$.
 */
double sgwb_phase_transition(double f, const double* params);

/**
 \brief default injection values for an SGWB
 */
void default_sgwb_injection(double* params, SGWB_t SGWB_type);

/**
 \brief default priors for an SGWB
 */


_Static_assert(SGWB_TEMPLATE_COUNT == 3, "Did you add an SGWB template? Add a default prior here.");
static const double default_powerlaw_prior[][2] = {
    // log Ap
    { -22.0, -7.0},
    // alpha_p
    { -2.0, 2.0},
};
static const double default_lognormal_prior[][2] = {
    // log10 A
    { -5.0, 0.0},
    // log10 fstar [Hz]
    { -5.0, 0.0},
    // log10 Delta
    { -2.0, 1.0},
};
static const double default_phase_transition_prior[][2] = {
    // rb
    { 0.1, 10.0},
    // b
    { 0.0, 8.0},
    // log10 Ap
    { -22.0, -4.0},
    // log10 fp [Hz]
    { -5.0, 0.0},
};
_Static_assert(SGWB_TEMPLATE_COUNT == 3, "Did you add an SGWB template? Edit this list of default priors, it needs to be exhaustive.");
static const double (*default_sgwb_priors[SGWB_TEMPLATE_COUNT])[2] = {
    [SGWB_TEMPLATE_POWERLAW] = default_powerlaw_prior,
    [SGWB_TEMPLATE_LOGNORMAL] = default_lognormal_prior,
    [SGWB_TEMPLATE_PHASE_TRANSITION] = default_phase_transition_prior
};


#endif /* noise_model_h */
