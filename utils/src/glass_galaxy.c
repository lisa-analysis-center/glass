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

#include "glass_utils.h"


double galaxy_distribution(double *x, double bulge_to_disk, double bulge_radius, double disk_radius, double disk_height)
{
    
    double z   = x[2];
    double u   = sqrt(x[0]*x[0]+x[1]*x[1]);
    double rsq = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];

    //unnormalized mass density of galaxy
    double p = 0.0;

    //contribution from bulge
    p += bulge_to_disk*exp(-rsq/(bulge_radius*bulge_radius));

    //contribution from disk
    double s = 1.0/cosh(z/disk_height);
    p += (1.0-bulge_to_disk)*exp(-u/disk_radius)*s*s;
    
    return(p);
}

double galaxy_foreground(double f, double A, double f1, double alpha, double fk, double f2)
{
    double Sf = A*pow(f,-7./3.) * exp(-pow(f/f1,alpha)) * 0.5*( 1. + tanh( (fk - f)/f2 ) );
    return Sf;
}

void rotate_galtoeclip(double *xg, double *xe)
{
    xe[0] = -0.05487556043*xg[0] + 0.4941094278*xg[1] - 0.8676661492*xg[2];
    
    xe[1] = -0.99382137890*xg[0] - 0.1109907351*xg[1] - 0.00035159077*xg[2];
    
    xe[2] = -0.09647662818*xg[0] + 0.8622858751*xg[1] + 0.4971471918*xg[2];
}

void rotate_ecliptogal(double *xg, double *xe)
{
    
    xg[0] = -0.05487556043*xe[0] - 0.99382137890*xe[1] - 0.09647662818*xe[2];
    
    xg[1] = 0.4941094278*xe[0] - 0.1109907351*xe[1] + 0.8622858751*xe[2];
    
    xg[2] = -0.8676661492*xe[0] - 0.00035159077*xe[1] + 0.4971471918*xe[2];
}


/*

*/

static void trig_factors_for_modulation(double p, double *cp, double *sp)
{
    // Maple likes to express everything in terms of cosines, so only need
    // linear in sine
    cp[1] = cos(p);
    sp[1] = sin(p);
    cp[2] = cp[1]*cp[1];
    cp[3] = cp[2]*cp[1];
    cp[4] = cp[2]*cp[2];
    cp[5] = cp[2]*cp[3];
    cp[6] = cp[3]*cp[3];
    cp[7] = cp[3]*cp[4];
    cp[8] = cp[4]*cp[4];
}

// Scaled Associated Legendre functions
static void plms(double ***Plm, double *zlist)
{
    int i, l, m;
    double x;
    double sq;
    double **prefac;
    
    // scaling prefactors for Plm's in spherical harmonics
    prefac = double_matrix(LMAX+1,LMAX+1);
    // Gamma(n) = (n-1)!
    for(l=0; l <= LMAX; l++)
    {
        for(m=0; m <= l; m++)
        {
            x = (double)(2*l+1)/(4.0*M_PI)*tgamma((double)(l-m+1))/tgamma((double)(l+m+1));
            prefac[l][m] =  sqrt(x);
        }
    }
    
    for(i=1; i < 4*NSIDE; i++)
    {
        x = zlist[i];
        sq = 1.0/sqrt(1.0-x*x);
        Plm[0][0][i] = 1.0;
        Plm[1][0][i] = x;
        // recurrence to get the m=0 terms for each l
        for(l=1; l < LMAX; l++)
        {
            Plm[l+1][0][i] = 1.0/(double)(l+1)*((double)(2*l+1)*x*Plm[l][0][i]-(double)(l)*Plm[l-1][0][i]);
        }
        // recurrence to get the m>0 terms for each l
        for(l=1; l <= LMAX; l++)
        {
            for(m=0; m < l; m++)
            {
                Plm[l][m+1][i] = sq*((double)(l-m)*x*Plm[l][m][i]-(double)(l+m)*Plm[l-1][m][i]);
            }
        }
        
        // apply the prefactor and the pixel area scaling
        
        x = 4.0*M_PI/(double)(12*NSIDE*NSIDE);
        
        for(l=0; l <= LMAX; l++)
        {
            for(m=0; m <= l; m++)
            {
                Plm[l][m][i] *= (x*prefac[l][m]);
            }
        }
    }
    
    free_double_matrix(prefac,LMAX+1);
    
    
}

// get the Plms on each ring. Includes pixel area factor
static void Plmrings(double ***Plm)
{
    double N3, NN3;
    long i;
    double *zlist;
    
    NN3 = 3.0*(double)(NSIDE*NSIDE);
    N3 = 3.0*(double)(NSIDE);
    
    zlist = double_vector(4*NSIDE);
    
    // z locations for each ring
    
    // North polar cap
    for(i=1; i< NSIDE; i++)
    {
        zlist[i] = 1.0 - (double)(i*i)/NN3;
    }
    
    //North equatorial belt
    for(i=NSIDE; i<= 2*NSIDE; i++)
    {
        zlist[i] = 4.0/3.0 - (double)(2*i)/N3;
    }
    
    //Mirror to the South
    for(i=2*NSIDE+1; i< 4*NSIDE; i++)
    {
        zlist[i] = -zlist[4*NSIDE-i];
    }
    
    // get the scaled Associated Legendre functions for each ring
    plms(Plm, zlist);
    
    free_double_vector(zlist);
    
}

#define sq5 2.2360679774997897
#define sq10 3.1622776601683793
#define sq15 3.8729833462074169
#define sq30 5.4772255750516611
#define sq35 5.916079783099616
#define sq70 8.3666002653407555
#define sq105 10.2469507659595984

static void XX(double **XXR, double **XXI, double *ca, double *sa, double *cb, double *sb)
{
    
    XXR[0][0] = (12.0*RTPI)/5.0;
    
    XXR[2][0] = -3.0*sq5*RTPI/35.0;
    XXR[2][1] = (9.0*sq10*RTPI*ca[1])/35.0;
    XXI[2][1] = -(9.0*sq10*RTPI*sa[1])/35.0;
    XXR[2][2] = (9.0*sq30*RTPI*(2.0*ca[2]-1.0))/70.0;
    XXI[2][2] = -(9.0*sq30*RTPI*(2.0*ca[1]*sa[1]))/70.0;
    
    XXR[4][0] = -9.0*RTPI*((ca[2] - 0.5)*(cb[2] - 0.5)*ca[1]*cb[1]*sb[1]*sa[1] + (cb[4] - cb[2] + 0.125)*ca[4] + (-cb[4] + cb[2] - 0.125)*ca[2] + cb[4]/8.0 - cb[2]/8.0 + 11.0/630.0);
    
    
    XXR[4][1] = -1.2*sq15*RTPI*(sb[1]*(cb[2]-0.5)*(ca[4]-1.5*ca[2]+0.25)*cb[1]*sa[1]+ca[1]*((cb[4]-cb[2]+0.125)*ca[4]+(-0.25-2.0*cb[4]+2.0*cb[2])*ca[2]+0.875*cb[4]-0.875*cb[2]+19.0/168.0));
    
    XXI[4][1] = -1.2*sq15*(((-cb[4]+cb[2]-0.125)*ca[4]+0.125*cb[4]-0.125*cb[2]+1.0/84.0)*sa[1]+sb[1]*ca[1]*(cb[2]-0.5)*(ca[4]-0.5*ca[2]-0.25)*cb[1])*RTPI;
    
    XXR[4][2] = -1.2*(sb[1]*ca[1]*(cb[2]-0.5)*cb[1]*(ca[4]-ca[2]+0.75)*sa[1]+(ca[2]-0.5)*((cb[4]-cb[2]+0.125)*ca[4]+(-cb[4]+cb[2]-0.125)*ca[2]+0.625*cb[4]-0.625*cb[2]+1.0/14.0))*RTPI*sq10;
    
    XXI[4][2] = -1.2*RTPI*sq10*(-((cb[4]-cb[2]+0.125)*ca[4]+(-cb[4]+cb[2]-0.125)*ca[2]-0.375*cb[4]+0.375*cb[2]-3.0/56.0)*ca[1]*sa[1]+(cb[2]-0.5)*cb[1]*(ca[2]-0.5)*(ca[4]-ca[2]-0.5)*sb[1]);
    
    XXR[4][3] =
    -8.0/35.0*(cb[1]*(cb[2]-0.5)*(ca[6]-1.25*ca[4]+0.375*ca[2]-7.0/16.0)*sb[1]*sa[1]+(-1.0/32.0+(cb[4]-cb[2]+0.125)*ca[6]+7.0/4.0*(-cb[4]+cb[2]-0.125)*ca[4]+0.125*(0.5+7.0*cb[4]-7.0*cb[2])*ca[2]-17.0/32.0*cb[4]+17.0/32.0*cb[2])*ca[1])*sq105*RTPI;
    
    XXI[4][3] =
    -8.0/35.0*sq105*RTPI*(((-ca[6]+1.25*ca[4]-0.375*ca[2]-13.0/32.0)*cb[4]+(ca[6]-1.25*ca[4]+0.375*ca[2]+13.0/32.0)*cb[2]-0.125*ca[6]+5.0/32.0*ca[4]-1.0/16.0)*sa[1]+cb[1]*(cb[2]-0.5)*(ca[6]-7.0/4.0*ca[4]+7.0/8.0*ca[2]+5.0/16.0)*sb[1]*ca[1]);
    
    XXR[4][4] = -4.0/35.0*(cb[1]*(cb[2]-0.5)*(ca[2]-0.5)*(ca[4]-ca[2]+0.125)*sb[1]*ca[1]*sa[1]+(ca[8]-2.0*ca[6]+1.25*ca[4]-0.25*ca[2]+41.0/64.0)*cb[4]+(-ca[8]+2.0*ca[6]-1.25*ca[4]+0.25*ca[2]-41.0/64.0)*cb[2]+0.125*ca[8]-0.25*ca[6]+1.0/64.0*ca[4]+7.0/64.0*ca[2]+1.0/16.0)*RTPI*sq70;
    
    
    XXI[4][4] = -4.0/35.0*RTPI*sq70*(-(ca[2]-0.5)*((cb[4]-cb[2]+0.125)*ca[4]+(-cb[4]+cb[2]-0.125)*ca[2]+0.125*cb[4]-0.125*cb[2]-0.125)*ca[1]*sa[1]+cb[1]*sb[1]*(ca[8]-2*ca[6]+1.25*ca[4]-0.25*ca[2]-0.625)*(cb[2]-0.5));
    
    
}

static void XY(double **XYR, double **XYI, double *ca, double *sa, double *cb, double *sb)
{
    
    XYR[0][0] = -(6.0*RTPI)/5.0;
    
    XYR[2][0] = 3.0*sq5*RTPI/70.0;
    XYR[2][1] = -(9.0*sq10*RTPI*ca[1])/70.0;
    XYI[2][1] = (9.0*sq10*RTPI*sa[1])/70.0;
    XYR[2][2] = -(9.0*sq30*RTPI*(2.0*ca[2]-1.0))/140.0;
    XYI[2][2] = (9.0*sq30*RTPI*(2.0*ca[1]*sa[1]))/140.0;
    
    XYR[4][0] = 4.5*RTPI*(((cb[3]-0.5*cb[1])*sb[1]+SQ3*(cb[4]-cb[2]+0.125))*(ca[2]-0.5)*ca[1]*sa[1]-cb[1]*SQ3*(ca[4]-ca[2]+0.125)*(cb[2]-0.5)*sb[1]+(cb[4]-cb[2]+0.125)*ca[4]+(-cb[4]+cb[2]-0.125)*ca[2]+11.0/630.0-0.125*cb[2]+0.125*cb[4]);
    
    XYR[4][1] = 0.6*sq5*(3.0*(1.0/3.0*SQ3*cb[1]*(cb[2]-0.5)*sb[1]+cb[4]-cb[2]+0.125)*(ca[4]-1.5*ca[2]+0.25)*sa[1]+(-3.0*(ca[4]-2.0*ca[2]+0.875)*cb[1]*(cb[2]-0.5)*sb[1]+SQ3*((cb[4]-cb[2]+0.125)*ca[4]+(-2.0*cb[4]+2.0*cb[2]-0.25)*ca[2]+19.0/168.0+0.875*cb[4]-0.875*cb[2]))*ca[1])*RTPI;
    
    XYI[4][1] = 9.0/5.0*sq5*((cb[1]*(ca[4]-0.125)*(cb[2]-0.5)*sb[1]-1.0/3.0*SQ3*((cb[4]-cb[2]+0.125)*ca[4]-0.125*cb[4]+0.125*cb[2]-1.0/84.0))*sa[1]+ca[1]*(1.0/3.0*SQ3*cb[1]*(cb[2]-0.5)*sb[1]+cb[4]-cb[2]+0.125)*(ca[4]-0.5*ca[2]-0.25))*RTPI;
    
    XYR[4][2] = 0.6*((ca[4]-ca[2]+0.75)*ca[1]*((cb[3]-0.5*cb[1])*sb[1]+SQ3*(cb[4]-cb[2]+0.125))*sa[1]+(ca[2]-0.5)*(-cb[1]*(cb[2]-0.5)*SQ3*(ca[4]-ca[2]+0.625)*sb[1]+(cb[4]-cb[2]+0.125)*ca[4]+(-cb[4]+cb[2]-0.125)*ca[2]+0.625*cb[4]-0.625*cb[2]+1.0/14.0))*sq10*RTPI;
    
    XYI[4][2] = 0.6*(-ca[1]*(-cb[1]*(cb[2]-0.5)*(ca[4]-ca[2]-0.375)*SQ3*sb[1]+(cb[4]-cb[2]+0.125)*ca[4]+(-cb[4]+cb[2]-0.125)*ca[2]-0.375*cb[4]+0.375*cb[2]-3.0/56.0)*sa[1]+((cb[3]-0.5*cb[1])*sb[1]+SQ3*(cb[4]-cb[2]+0.125))*(ca[4]-ca[2]-0.5)*(ca[2]-0.5))*sq10*RTPI;
    
    XYR[4][3] = 4.0/35.0*sq35*RTPI*(3.0*(ca[6]-1.25*ca[4]+0.375*ca[2]-7.0/16.0)*(1.0/3.0*SQ3*cb[1]*(cb[2]-0.5)*sb[1]+cb[4]-cb[2]+0.125)*sa[1]+(-3.0*cb[1]*(cb[2]-0.5)*(ca[6]-7.0/4.0*ca[4]+0.875*ca[2]-17.0/32.0)*sb[1]+SQ3*(-1.0/32.0+(cb[4]-cb[2]+0.125)*ca[6]+7.0/4.0*(-cb[4]+cb[2]-0.125)*ca[4]+0.125*(0.5+7.0*cb[4]-7.0*cb[2])*ca[2]-17.0/32.0*cb[4]+17.0/32.0*cb[2]))*ca[1]);
    
    
    XYI[4][3] = 4.0/35.0*((((-ca[6]+1.25*ca[4]-0.375*ca[2]-13.0/32.0)*cb[4]+(ca[6]-1.25*ca[4]+0.375*ca[2]+13.0/32.0)*cb[2]-0.125*ca[6]+5.0/32.0*ca[4]-1.0/16.0)*SQ3+3.0*(cb[2]-0.5)*(ca[6]-1.25*ca[4]+0.375*ca[2]+13.0/32.0)*cb[1]*sb[1])*sa[1]+(ca[6]-7.0/4.0*ca[4]+0.875*ca[2]+5.0/16.0)*ca[1]*(0.5*cb[1]*sb[1]*(2.0*cb[2]-1.0)*SQ3+3.0*cb[4]-3.0*cb[2]+0.375))*sq35*RTPI;
    
    
    XYR[4][4] = -2.0/35.0*(-(SQ3*(cb[4]-cb[2]+0.125)+0.5*sb[1]*cb[1]*(2.0*cb[2]-1.0))*(ca[2]-0.5)*ca[1]*(ca[4]-ca[2]+0.125)*sa[1]+SQ3*(ca[8]-2.0*ca[6]+1.25*ca[4]-0.25*ca[2]+41.0/64.0)*(cb[2]-0.5)*cb[1]*sb[1]+(-ca[8]+2.0*ca[6]-1.25*ca[4]+0.25*ca[2]-41.0/64.0)*cb[4]+(ca[8]-2.0*ca[6]+1.25*ca[4]-0.25*ca[2]+41.0/64.0)*cb[2]-1.0/16.0+0.25*ca[6]-1.0/64.0*ca[4]-7.0/64.0*ca[2]-0.125*ca[8])*sq70*RTPI;
    
    XYI[4][4] = 2.0/35.0*(-(ca[2]-0.5)*(-cb[1]*SQ3*(ca[4]-ca[2]+0.125)*(cb[2]-0.5)*sb[1]+(cb[4]-cb[2]+0.125)*ca[4]+(-cb[4]+cb[2]-0.125)*ca[2]+0.125*cb[4]-0.125*cb[2]-0.125)*ca[1]*sa[1]+(ca[8]-2*ca[6]+1.25*ca[4]-0.25*ca[2]-0.625)*((cb[3]-0.5*cb[1])*sb[1]+SQ3*(cb[4]-cb[2]+0.125)))*sq70*RTPI;
    
}

void initialize_galaxy_modulation(struct GalaxyModulation *gm, struct Wavelets *wdm, struct Orbit *orbit, double Tobs, double t0)
{
    double *cbx, *sbx;
    double *cby, *sby;
    double *cbz, *sbz;
    double *ca, *sa;
    double alpha, beta;
    int i;
            
    gm->alpha_0  = (t0/YEAR)*PI2;   //initial angle of orbit
    gm->alphamax = ((t0+Tobs)/YEAR)*PI2; // angle through which orbit completes

    // mark per-Q modulation cache as empty
    gm->Q_cached = 0;
    gm->Nt_cached = 0;
    gm->Nt_alloc = 0;
    gm->cache_XX = NULL;
    gm->cache_YY = NULL;
    gm->cache_ZZ = NULL;
    gm->cache_XY = NULL;
    gm->cache_XZ = NULL;
    gm->cache_YZ = NULL;

    gm->N = (int)((Tobs/YEAR)*100.0);  // numer of samples
    
    gm->t = double_vector(gm->N);
    double dt = Tobs/gm->N;

    //pad time samples for interpolants
    dt = (Tobs + 2*dt) / (gm->N - 1);

    for(int n=0; n<gm->N; n++) gm->t[n] = t0 + dt*n - dt;

    // The Plms are always the same for a fixed values of NSIDE and LMAX
    // Precompute and store
    gm->Plm = double_tensor(LMAX+1, LMAX+1, 4*NSIDE);
    Plmrings(gm->Plm);
    
    gm->XXR = double_tensor(gm->N,LMAX+1, LMAX+1);
    gm->XXI = double_tensor(gm->N,LMAX+1, LMAX+1);
    gm->YYR = double_tensor(gm->N,LMAX+1, LMAX+1);
    gm->YYI = double_tensor(gm->N,LMAX+1, LMAX+1);
    gm->ZZR = double_tensor(gm->N,LMAX+1, LMAX+1);
    gm->ZZI = double_tensor(gm->N,LMAX+1, LMAX+1);
    gm->XYR = double_tensor(gm->N,LMAX+1, LMAX+1);
    gm->XYI = double_tensor(gm->N,LMAX+1, LMAX+1);
    gm->YZR = double_tensor(gm->N,LMAX+1, LMAX+1);
    gm->YZI = double_tensor(gm->N,LMAX+1, LMAX+1);
    gm->XZR = double_tensor(gm->N,LMAX+1, LMAX+1);
    gm->XZI = double_tensor(gm->N,LMAX+1, LMAX+1);
        
    
    cbx = double_vector(9);
    sbx = double_vector(9);
    cby = double_vector(9);
    sby = double_vector(9);
    cbz = double_vector(9);
    sbz = double_vector(9);
    ca  = double_vector(9);
    sa  = double_vector(9);
    
    beta = orbit->lambda_0;
    trig_factors_for_modulation(beta, cbx, sbx);
    beta = 2.0*M_PI/3.0 + orbit->lambda_0;
    trig_factors_for_modulation(beta, cby, sby);
    beta = 4.0*M_PI/3.0 + orbit->lambda_0;
    trig_factors_for_modulation(beta, cbz, sbz);
    
    // for the cross terms, putting beta = 2Pi/3 into XY yields YZ and beta = 4Pi/3 in XY yields XZ

    
    for(i=0; i<gm->N; i++)
    {
        alpha = gm->alpha_0 + (double)(i)*gm->alphamax/(double)(gm->N-1) + orbit->kappa_0;
        
        trig_factors_for_modulation(alpha, ca, sa);
        
        XX(gm->XXR[i], gm->XXI[i], ca, sa, cbx, sbx);
        XX(gm->YYR[i], gm->YYI[i], ca, sa, cby, sby);
        XX(gm->ZZR[i], gm->ZZI[i], ca, sa, cbz, sbz);
        
        XY(gm->XYR[i], gm->XYI[i], ca, sa, cbx, sbx);
        XY(gm->YZR[i], gm->YZI[i], ca, sa, cby, sby);
        XY(gm->XZR[i], gm->XZI[i], ca, sa, cbz, sbz);
         
    }
    
    free(cbx);
    free(sbx);
    free(cby);
    free(sby);
    free(cbz);
    free(sbz);
    free(ca);
    free(sa);
    
    //interpolants for modulation pattern
    gm->XX_spline = alloc_cubic_spline(gm->N);
    gm->YY_spline = alloc_cubic_spline(gm->N);
    gm->ZZ_spline = alloc_cubic_spline(gm->N);
    gm->XY_spline = alloc_cubic_spline(gm->N);
    gm->XZ_spline = alloc_cubic_spline(gm->N);
    gm->YZ_spline = alloc_cubic_spline(gm->N);


    double theta, phi;

    gm->Npix = 12*NSIDE*NSIDE;
           
    // set up sky grid in ecliptic coordinates
    
    gm->skytheta = (double*)malloc(sizeof(double)*(gm->Npix));
    gm->skyphi   = (double*)malloc(sizeof(double)*(gm->Npix));
    
    // galaxy in ecliptic coordinates
    
    double xe[3];
    double xg[3];
    


    for(i=0; i<gm->Npix; i++)
    {
        // get theta, phi coordinates of pixel (ecliptic)
        astropy_pix2ang_ring(NSIDE,  i, &theta, &phi);
        
        xe[0] = sin(theta)*cos(phi);
        xe[1] = sin(theta)*sin(phi);
        xe[2] = cos(theta);
        
        // get the theta and phi in galactic coords
        rotate_ecliptogal(xg, xe);
        theta = acos(xg[2]);
        phi = atan2(xg[1],xg[0]);
        if(phi<0.0) phi += 2.0*M_PI;
        
        gm->skytheta[i] = theta;
        gm->skyphi[i]   = phi;
        
    }
        
}



static double galaxy_integrand(double *params, double r, double sintheta, double costheta, double cosphi)
{    
    // A  bulge fraction
    // Rb  bulge radius (kpc)
    // Rd  disk radius (kpc)
    // Zd  disk height (kpc)
    // Rgc distance from SSB to GC (kpc)
    double A = params[0];
    double Rb = params[1];
    double Rd = params[2];
    double Zd = params[3];
    double Rgc = params[4];
    double *x = malloc(3*sizeof(double));
    
    //convert from spherical SSB galactic coordinates to cartesian galacto-centric coordinates
    double sinphi = sqrt(1.0-cosphi*cosphi);
    x[0] = r*cosphi*sintheta-Rgc;
    x[1] = r*sinphi*sintheta;
    x[2] = r*costheta;

    // unnomralized galactic mass density
    double rho = galaxy_distribution(x, A, Rb, Rd, Zd);
    
    free(x);
    return(rho);
}

static double galaxy_integration(double *params, double theta, double phi)
{
    double IG;
    double rmax;
    double r, rr, h, dh;
    double s5, s3;
    int i;
    double err, ferr, tol, min;
    
    double sintheta = sin(theta);
    double costheta = sqrt(1.0-sintheta*sintheta);
    double cosphi   = cos(phi);
    
    tol = 1.0e-6;
    min = 1.0e-10;
    
    double I5[5];
    double I3[3];
    
    rmax = 200.0;
    h = rmax/10000.0;

    r = params[5]; //integration starts at Rcut (everything closer is resolved)
    
    IG = 0.0;
    
    do
    {
        do
        {
            dh = h/4.0;
            
            for(i=0; i<5; i++)
            {
                rr = r + (double)(i)*dh;
                I5[i] = galaxy_integrand(params, rr, sintheta, costheta, cosphi);
            }
            for(i=0; i<3; i++)
            {
                I3[i] = I5[2*i];
            }
            
            // 3 point Simpson
            s3 = simpson_integration_3(I3[0],I3[1], I3[2], h);

            // 5 point Simpson
            s5 = simpson_integration_5(I5[0],I5[1],I5[2],I5[3],I5[4],h);

            // absolute  error
            err = fabs(s5-s3)/(15.0);

            // fractional  error
            ferr = err/s5;
            
            if(ferr > tol && err > min) h /= 2.0;
                        
        }while(ferr > tol && err > min);
        
        r += h;
        IG += s5;
        
        // we try a larger step for the next iteration
        // this might then have to be shrunk
        h *= 2.0;
        
    }while(r < rmax);
    return IG;
}

static void alm2map(double *skyrecon, double **almR, double **almI, double ***Plm)
{
    
    long pix, Npix;
    double theta, phi, thold, x;
    int i, l, m;
    double *cm, *sm;

    cm = double_vector(LMAX+1);
    sm = double_vector(LMAX+1);
    
    Npix = 12*NSIDE*NSIDE;
    
    x = (double)(Npix)/(4.0*M_PI);
    
    astropy_pix2ang_ring(NSIDE, 0, &thold, &phi);
    i = 1;
    
    for(pix=0; pix < Npix; pix++)
    {
        astropy_pix2ang_ring(NSIDE, pix, &theta, &phi);
        
        if(fabs(theta-thold) > 1.0e-6)
        {
            i++;  // moved to the next ring
            thold = theta;
        }
        
        // m=0 terms
        for(l=0; l <= LMAX; l++) skyrecon[pix] += almR[l][0]*Plm[l][0][i];
       
        cm[1] = cos(phi);
        sm[1] = sin(phi);
        
        // cos((m+1)phi) = cos(m*phi)cos(phi) - sin(m*phi)sin(phi)
        // sin((m+1)phi) = sin(m*phi)*cos(phi)+cos(m*pi)*sin(phi);
        for(m=1; m < LMAX; m++)
        {
            cm[m+1] = cm[m]*cm[1]-sm[m]*sm[1];
            sm[m+1] = sm[m]*cm[1]+cm[m]*sm[1];
        }
        
        for(l=0; l <= LMAX; l++)
        {
            for(m=1; m <= l; m++)
            {
                skyrecon[pix] += 2.0*(almR[l][m]*Plm[l][m][i]*cm[m]+almI[l][m]*Plm[l][m][i]*sm[m]);
            }
        }
        
        skyrecon[pix] *= x;
        
    }
    
    free_double_vector(cm);
    free_double_vector(sm);
    
}

static void map2alm(double *sky, double **almR, double **almI, double ***Plm)
{
        
        long pix, Npix;
        double theta, phi, thold;
        int i, l, m;
        double *cm, *sm;
    
        cm = double_vector(LMAX+1);
        sm = double_vector(LMAX+1);
        
        Npix = 12*NSIDE*NSIDE;
        
        // we use this to keep track of which theta ring were are on (not very efficient)
        astropy_pix2ang_ring(NSIDE, 0, &thold, &phi);
        i = 1;
        
        for(pix=0; pix < Npix; pix++)
        {
            astropy_pix2ang_ring(NSIDE, pix, &theta, &phi);
            
            if(fabs(theta-thold) > 1.0e-6)
            {
                i++;  // moved to the next ring
                thold = theta;
            }
            
            cm[0] = 1.0;
            sm[0] = 0.0;
            cm[1] = cos(phi);
            sm[1] = sin(phi);
            
            // cos((m+1)phi) = cos(m*phi)cos(phi) - sin(m*phi)sin(phi)
            // sin((m+1)phi) = sin(m*phi)*cos(phi)+cos(m*pi)*sin(phi);
            for(m=1; m < LMAX; m++)
            {
                cm[m+1] = cm[m]*cm[1]-sm[m]*sm[1];
                sm[m+1] = sm[m]*cm[1]+cm[m]*sm[1];
            }
            
            for(l=0; l <= LMAX; l++)
            {
                for(m=0; m <= l; m++)
                {
                    almR[l][m] +=  Plm[l][m][i]*cm[m]*sky[pix];
                    almI[l][m] +=  Plm[l][m][i]*sm[m]*sky[pix];
                }
            }
            
        }
    
     free_double_vector(cm);
     free_double_vector(sm);
        
}

void sphharm(double ***Plm, double **almR, double **almI, double *sky)
{
    long i;
    int l, m;
    double *skyrecon;
    long Npix, pix;
   
    // This slow version uses a simple sum over the spherical harmonics. The analysis uses two
    // iterations to reduce the errors from the uneven pixel shapes. The errors are corrected
    // by reconsructing the map from the initial alms then transforming initial minus reconstructed
    // and adding these values to the original alm estimates.
    
    Npix = 12*NSIDE*NSIDE;
    
    // reconstructed map from the initial alms
    skyrecon = double_vector((int)Npix);
    
    // initialize the alms
    for(l=0; l <= LMAX; l++)
    {
        for(m=0; m <= l; m++)
        {
            almR[l][m] =  0.0;
            almI[l][m] =  0.0;
        }
    }
    
    map2alm(sky, almR, almI, Plm);
    
    // iterate the transform to correct pixel shape errors
    for(i=0; i < 2; i++)
    {
        
        // reconstruct the sky using the alms
        for(pix=0; pix < Npix; pix++) skyrecon[pix] = 0.0;
        alm2map(skyrecon, almR, almI, Plm);
        // transform the difference to correct the errors
        for(pix=0; pix < Npix; pix++) skyrecon[pix] = sky[pix] - skyrecon[pix];
        map2alm(skyrecon, almR, almI, Plm);
        
    }
    
    // minus sign on complex due to some convention difference with Healpix
    for(l=0; l <= LMAX; l++)
    {
        for(m=0; m <= l; m++)
        {
            almI[l][m] = -almI[l][m];
        }
    }

    free(skyrecon);
    
}

void galaxy_modulation(struct GalaxyModulation *gm, double *params)
{
    double av;
    
    // galaxy in ecliptic coordinates
    double *skyeclip = double_vector((int)gm->Npix);
    for(int i=0; i<gm->Npix; i++)
        skyeclip[i] = galaxy_integration(params,gm->skytheta[i],gm->skyphi[i]);

    double **almR = double_matrix(LMAX+1, LMAX+1);
    double **almI = double_matrix(LMAX+1, LMAX+1);
    
    sphharm(gm->Plm, almR, almI, skyeclip);
    
    av = 0.0;
    
    double *xx = double_vector(gm->N);
    double *yy = double_vector(gm->N);
    double *zz = double_vector(gm->N);
    double *xy = double_vector(gm->N);
    double *yz = double_vector(gm->N);
    double *xz = double_vector(gm->N);

    for(int i=0; i<gm->N; i++)
    {
        
        xx[i] = 0.0;
        yy[i] = 0.0;
        zz[i] = 0.0;
        xy[i] = 0.0;
        yz[i] = 0.0;
        xz[i] = 0.0;
        
        //m=0 terms
        for(int l=0; l<=LMAX; l++)
        {
            xx[i] += almR[l][0]*gm->XXR[i][l][0];
            yy[i] += almR[l][0]*gm->YYR[i][l][0];
            zz[i] += almR[l][0]*gm->ZZR[i][l][0];
            xy[i] += almR[l][0]*gm->XYR[i][l][0];
            yz[i] += almR[l][0]*gm->YZR[i][l][0];
            xz[i] += almR[l][0]*gm->XZR[i][l][0];
        }
        
        //m!=0 terms
        for(int l=0; l<=LMAX; l++)
        {
            for(int m=1; m <= l; m++)
            {
                xx[i] += 2.0*(almR[l][m]*gm->XXR[i][l][m]+almI[l][m]*gm->XXI[i][l][m]);
                yy[i] += 2.0*(almR[l][m]*gm->YYR[i][l][m]+almI[l][m]*gm->YYI[i][l][m]);
                zz[i] += 2.0*(almR[l][m]*gm->ZZR[i][l][m]+almI[l][m]*gm->ZZI[i][l][m]);
                xy[i] += 2.0*(almR[l][m]*gm->XYR[i][l][m]+almI[l][m]*gm->XYI[i][l][m]);
                yz[i] += 2.0*(almR[l][m]*gm->YZR[i][l][m]+almI[l][m]*gm->YZI[i][l][m]);
                xz[i] += 2.0*(almR[l][m]*gm->XZR[i][l][m]+almI[l][m]*gm->XZI[i][l][m]);
            }
        }
                
        av += (xx[i]+yy[i]+zz[i])/3.0;
    }
    
    av /= (double)(gm->N);
    
    for(int i=0; i< gm->N; i++)
    {
        xx[i] /= av;
        yy[i] /= av;
        zz[i] /= av;
        xy[i] /= av;
        yz[i] /= av;
        xz[i] /= av;
    }
    
    initialize_cubic_spline(gm->XX_spline, gm->t, xx, SPLINE_EVEN_SAMPLED);
    initialize_cubic_spline(gm->YY_spline, gm->t, yy, SPLINE_EVEN_SAMPLED);
    initialize_cubic_spline(gm->ZZ_spline, gm->t, zz, SPLINE_EVEN_SAMPLED);
    initialize_cubic_spline(gm->XY_spline, gm->t, xy, SPLINE_EVEN_SAMPLED);
    initialize_cubic_spline(gm->XZ_spline, gm->t, xz, SPLINE_EVEN_SAMPLED);
    initialize_cubic_spline(gm->YZ_spline, gm->t, yz, SPLINE_EVEN_SAMPLED);

    // splines moved; downstream cache is stale
    gm->Q_cached = 0;

    FILE *out = fopen("modulation.dat", "w");
    for(int i=0; i<gm->N; i++)
    {
        double alpha = gm->alpha_0 + i*gm->alphamax/(double)(gm->N-1);
        fprintf(out,"%f %.10f %.10f %.10f %.10f %.10f %.10f\n", alpha, xx[i], yy[i], zz[i], xy[i], xz[i], yz[i]);
    }
    fclose(out);
   

    free_double_matrix(almR,LMAX+1);
    free_double_matrix(almI,LMAX+1);
    free_double_vector(skyeclip);
    free_double_vector(xx);
    free_double_vector(yy);
    free_double_vector(zz);
    free_double_vector(xy);
    free_double_vector(xz);
    free_double_vector(yz);
}

void galaxy_modulation_cache_for_Q(struct GalaxyModulation *gm, struct Wavelets *wdm, int Q)
{
    int Nt = wdm->NT / Q;
    if (gm->Q_cached == Q && gm->Nt_cached == Nt) return;

    if (Nt > gm->Nt_alloc) {
        free(gm->cache_XX); free(gm->cache_YY); free(gm->cache_ZZ);
        free(gm->cache_XY); free(gm->cache_XZ); free(gm->cache_YZ);
        gm->cache_XX = malloc(Nt * sizeof(double));
        gm->cache_YY = malloc(Nt * sizeof(double));
        gm->cache_ZZ = malloc(Nt * sizeof(double));
        gm->cache_XY = malloc(Nt * sizeof(double));
        gm->cache_XZ = malloc(Nt * sizeof(double));
        gm->cache_YZ = malloc(Nt * sizeof(double));
        gm->Nt_alloc = Nt;
    }

    // Cell-AVERAGE of the modulation over the Q fine pixels in each coarse cell,
    // i.e. (1/Q) sum_{i in cell} m(i*dt). This matches the data statistic
    // P_q = (1/Q) sum_i w_i^2 whose expectation is <S>_cell, the cell-average of
    // the per-pixel variance. The previous code sampled only the cell midpoint;
    // for a curved (annual) modulation that biases the modeled galactic power by
    // m(midpoint) - <m>_cell -- systematically low near the annual peaks where the
    // galaxy is brightest and the modulation is most concave. Averaging removes
    // that bias to all orders. Fine pixel i sits at t = i*dt (same mapping the
    // Q=1 path uses), so for Q=1 the average collapses to the single-pixel value
    // and this is bit-for-bit identical to the old midpoint behavior.
    double invQ = 1.0 / (double)Q;
    for (int q = 0; q < Nt; q++) {
        double axx=0, ayy=0, azz=0, axy=0, axz=0, ayz=0;
        for (int i = q*Q; i < (q+1)*Q; i++) {
            double t = i * wdm->dt;
            axx += spline_interpolation(gm->XX_spline, t);
            ayy += spline_interpolation(gm->YY_spline, t);
            azz += spline_interpolation(gm->ZZ_spline, t);
            axy += spline_interpolation(gm->XY_spline, t);
            axz += spline_interpolation(gm->XZ_spline, t);
            ayz += spline_interpolation(gm->YZ_spline, t);
        }
        gm->cache_XX[q] = axx*invQ;
        gm->cache_YY[q] = ayy*invQ;
        gm->cache_ZZ[q] = azz*invQ;
        gm->cache_XY[q] = axy*invQ;
        gm->cache_XZ[q] = axz*invQ;
        gm->cache_YZ[q] = ayz*invQ;
    }

    gm->Q_cached = Q;
    gm->Nt_cached = Nt;
}

