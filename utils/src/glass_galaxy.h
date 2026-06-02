/*
 * Copyright 2024 Tyson B. Littenberg & Neil J. Cornish
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
@file glass_galaxy.h
\brief Modulated astrophysical foreground noise from unresolved galactic binaries.
*/

///@name Galaxy Parameters
///@{
#define GALAXY_RGC 7.2 //!< distance from SSB to GC (kpc)
#define GALAXY_A  0.25 //!< bulge fraction
#define GALAXY_Rb 0.8  //!< bulge radius (kpc)
#define GALAXY_Rd 2.5  //!< disk radius (kpc)
#define GALAXY_Zd 0.4  //!< disk height (kpc)
///@}

#define NSIDE 16 //!< Healpix resolution for galaxy modulation calculations
#define LMAX 4   //!< Maximum l for spherical harmonic decomposition of galaxy modulation

/*!
 * \brief Metadata for constructing modulated foreground noise.= *
 */
struct GalaxyModulation
{
    int N;           //!<Number of time samples of orbit
    double *t;       //!<Time grid of orbit samples
    double alpha_0;
    double alphamax;
    double ***Plm;
    double ***XXR;
    double ***XXI;
    double ***YYR;
    double ***YYI;
    double ***ZZR;
    double ***ZZI;
    double ***XYR;
    double ***XYI;
    double ***YZR;
    double ***YZI;
    double ***XZR;
    double ***XZI;

    struct CubicSpline *XX_spline; //!<Modulation in PSD of X channel
    struct CubicSpline *YY_spline; //!<Modulation in PSD of Y channel
    struct CubicSpline *ZZ_spline; //!<Modulation in PSD of X channel
    struct CubicSpline *XY_spline; //!<Modulation in CSD of XY channels
    struct CubicSpline *XZ_spline; //!<Modulation in PSD of XZ channels
    struct CubicSpline *YZ_spline; //!<Modulation in PSD of YZ channels

    /* Spline evaluations cached on the wavelet time grid for a given coarsening
     * factor Q. The 6 modulation splines are stationary across an MCMC run, so
     * we evaluate them once per (wdm, Q) pair and reuse. cache_*[q] is the
     * cell-AVERAGE of the spline over the Q fine pixels in coarse cell q,
     * (1/Q) sum_{i=q*Q}^{(q+1)*Q-1} m(i*dt), which collapses to the single
     * spline value m(q*dt) when Q = 1. This matches the cell-averaged data
     * statistic (see galaxy_modulation_cache_for_Q). Q_cached = 0 means
     * uninitialized. */
    int Q_cached;
    int Nt_cached;
    int Nt_alloc;
    double *cache_XX;
    double *cache_YY;
    double *cache_ZZ;
    double *cache_XY;
    double *cache_XZ;
    double *cache_YZ;

    long Npix;        //!<Number of pixels on the healpix grid
    double *skytheta; //!<Latitude value of healpix grid pixels
    double *skyphi;   //!<Longitude value of healpix grid pixels
};

double galaxy_distribution(double *x, double bulge_to_disk, double bulge_radius, double disk_radius, double disk_height);
double galaxy_foreground(double f, double A, double f1, double alpha, double fk, double f2);
void rotate_galtoeclip(double *xg, double *xe);
void rotate_ecliptogal(double *xg, double *xe);

void initialize_galaxy_modulation(struct GalaxyModulation *gm, struct Wavelets *wdm, struct Orbit *orbit, double Tobs, double t0);
void free_galaxy_modulation(struct GalaxyModulation *gm);
void galaxy_modulation(struct GalaxyModulation *gm, double *params);

// Build (or refresh) the per-time-slice modulation cache for the given (wdm, Q).
void galaxy_modulation_cache_for_Q(struct GalaxyModulation *gm, struct Wavelets *wdm, int Q);
