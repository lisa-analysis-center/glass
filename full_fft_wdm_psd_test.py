import numpy as np
from WDMWaveletTransforms.wavelet_transforms import (
    transform_wavelet_freq,
    transform_wavelet_freq_time,
    transform_wavelet_time,
    inverse_wavelet_freq,
    phitilde_vec_norm
)
import matplotlib.pyplot as plt
from WDMWaveletTransforms.transform_freq_funcs import phitilde_vec_norm

def randc(shape):
    return np.random.randn(*shape) + 1j*np.random.randn(*shape)

def testpsdfunc(x):
    fref = 1e-3
    c = 1e-30
    d = 1e-20
    return c*(x/fref)**2 + d*np.sin(x/fref)**2
def testpsdfunc_nonstat(f,t):
    fwiggle = 2.0/(3600*24*30)
    omega_wiggle = 2*np.pi*fwiggle
    fref = 1e-3 + 2e-4*np.sin(omega_wiggle*t)
    c = 1e-30
    d = 1e-20
    return c*(f/fref)**2 + d*np.sin(f/fref)**2
def wdm_psd_nonstat(psdfunc, NF, NT, dt, nx=4.0):
    """WDM-domain PSD from one-sided FFT PSD S_n(f_k) [x^2/Hz].

    Returns shape (NT, NF) array of E[w_{nm}^2].
    """
    ND = NF * NT
    phif = phitilde_vec_norm(NF, NT, nx)
    half = NT // 2
    nmax = ND // 2
    NFFT = nmax + 1
    f_dft = np.arange(NFFT) * 1.0/(NF*NT*dt)
    DeltaT = NF*dt

    # Pack into (NT, NF)
    wdmpsd = np.zeros((NT, NF))
    # Compute variance for all NF+1 transform layers (0=DC ... NF=Nyquist)
    for i in range(NT):
        ti = (i+0.5)*DeltaT
        fftpsd = psdfunc(f_dft, ti)
        centers = np.arange(NF + 1) * half
        layer_psd = np.zeros(NF + 1)
        for l in range(-(half - 1), half):
            k = centers + l
            # Fold into [0, nmax] via reflection (real signal: S(-f)=S(f))
            k = np.abs(k)
            k = k % (2 * nmax)
            k = np.where(k > nmax, 2 * nmax - k, k)
            layer_psd += phif[abs(l)] ** 2 * fftpsd[k]
        layer_psd /= dt * NT * NF
        if i == 0:
            wdmpsd[0::2, 0] = layer_psd[0]           # DC at even rows
            wdmpsd[1::2, 0] = layer_psd[NF]          # Nyquist at odd rows
        else:
            wdmpsd[i, 1:] = layer_psd[1:NF]          # regular layers

    return wdmpsd

def wdm_psd(fftpsd, NF, NT, dt, nx=4.0):
    """WDM-domain PSD from one-sided FFT PSD S_n(f_k) [x^2/Hz].

    Returns shape (NT, NF) array of E[w_{nm}^2].
    """
    ND = NF * NT
    phif = phitilde_vec_norm(NF, NT, nx)
    half = NT // 2
    nmax = ND // 2

    # Compute variance for all NF+1 transform layers (0=DC ... NF=Nyquist)
    centers = np.arange(NF + 1) * half
    layer_psd = np.zeros(NF + 1)
    for l in range(-(half - 1), half):
        k = centers + l
        # Fold into [0, nmax] via reflection (real signal: S(-f)=S(f))
        k = np.abs(k)
        k = k % (2 * nmax)
        k = np.where(k > nmax, 2 * nmax - k, k)
        layer_psd += phif[abs(l)] ** 2 * fftpsd[k]
    layer_psd /= dt * NT * NF

    # Pack into (NT, NF)
    wdmpsd = np.zeros((NT, NF))
    wdmpsd[:, 1:] = layer_psd[1:NF]          # regular layers
    wdmpsd[0::2, 0] = layer_psd[0]           # DC at even rows
    wdmpsd[1::2, 0] = layer_psd[NF]          # Nyquist at odd rows

    return wdmpsd

def nonstat_test(plots=True):
    WAVELET_DURATION = 7680
    NF = 1536
    NT = 338
    ND = NT*NF
    NFFT = ND/2 + 1
    Tobs = NT*WAVELET_DURATION
    dt = WAVELET_DURATION / NF
    DeltaF = 1/(2*dt*NF)
    DeltaT = WAVELET_DURATION
    # we need to generate this data in WDM actually
    t = np.arange(NT) * WAVELET_DURATION
    f = np.arange(NF) * DeltaF

    fcenter = f + DeltaF/2
    tcenter = t + DeltaT/2

    f_dft = np.arange(NFFT) * 1/Tobs
    if plots:
        for ti in tcenter[::50]:
            testpsd = testpsdfunc_nonstat(fcenter, ti)
            plt.loglog(fcenter, testpsd, label=f"t={ti/3600/24} days")
            plt.xlabel("f [Hz]")
            plt.ylabel("PSD [1/Hz]")
        plt.legend()
        plt.show()
    testwdmpsd_approx = np.zeros((NT,NF))
    for i, ti in enumerate(tcenter):
        testwdmpsd_approx[i,:] = testpsdfunc_nonstat(fcenter,ti)/ (2*dt)
    testwdmpsd_fullf = wdm_psd_nonstat(testpsdfunc_nonstat, NF, NT, dt) 

    wdmgen = np.random.randn(*testwdmpsd_fullf.shape) * np.sqrt(testwdmpsd_fullf)
    if plots:
        plt.imshow(np.sqrt(testwdmpsd_approx.T), aspect='auto', origin='lower',extent=(t[0],t[-1],f[0],f[-1]), norm='log')
        plt.colorbar()
        plt.title("Approximate wASD")
        plt.xlabel("t [s]")
        plt.ylabel("f [Hz]")
        plt.savefig("wpsd_approx_nonstat.png")
        plt.show()
        plt.imshow(np.sqrt(testwdmpsd_fullf.T), aspect='auto', origin='lower',extent=(t[0],t[-1],f[0],f[-1]), norm='log')
        plt.colorbar()
        plt.title("f-integrated wASD")
        plt.xlabel("t [s]")
        plt.ylabel("f [Hz]")
        plt.savefig("wpsd_nonstat.png")
        plt.show()
        plt.imshow(np.abs(wdmgen.T), aspect='auto', origin='lower',extent=(t[0],t[-1],f[0],f[-1]), norm='log')
        plt.colorbar()
        plt.title("WDM non-stationary data (|w_nm|)")
        plt.xlabel("t [s]")
        plt.ylabel("f [Hz]")
        plt.savefig("wdm_data_nonstat.png")
        plt.show()
        plt.imshow(wdmgen.T / np.sqrt(testwdmpsd_approx.T), aspect='auto', origin='lower',extent=(t[0],t[-1],f[0],f[-1]), norm='linear')
        plt.colorbar()
        plt.title("Whitened WDM data with approx wPSD")
        plt.savefig("whitened_approx_nonstat.png")
        plt.show()
        plt.imshow(wdmgen.T / np.sqrt(testwdmpsd_fullf.T), aspect='auto', origin='lower',extent=(t[0],t[-1],f[0],f[-1]), norm='linear')
        plt.colorbar()
        plt.title("Whitened WDM data with approx wPSD")
        plt.savefig("whitened_nonstat.png")
        plt.show()


def stat_test(plots=True):
    WAVELET_DURATION = 7680
    NF = 1536
    NT = 338
    ND = NT*NF
    NFFT = ND/2 + 1
    Tobs = NT*WAVELET_DURATION
    dt = WAVELET_DURATION / NF
    DeltaF = 1/(2*dt*NF)

    testf = np.arange(NFFT) * 1/Tobs

    testpsd = testpsdfunc(testf)
    testfft = randc(testpsd.shape) * np.sqrt(testpsd/2 * ND/(2*dt))

    if plots:
        plt.loglog(testf, np.abs(testfft))
        plt.loglog(testf, np.sqrt(testpsd*ND/(2*dt)))
        plt.title("Freq domain FFT and PSD")
        plt.ylabel("PSD [1/Hz]")
        plt.xlabel("f [Hz]")
        plt.savefig("fft_data.png")
        plt.show()

    testwdm = transform_wavelet_freq(testfft, NF, NT)

    t = np.arange(testwdm.shape[0]) * WAVELET_DURATION
    f = np.arange(testwdm.shape[1]) * DeltaF

    # Approximate WDM PSD
    fcenter = f + DeltaF/2
    testwdmpsd_approx = np.repeat([testpsdfunc(fcenter)], NT, axis=0) / (2*dt)

    # Exact WDM PSD
    testwdmpsd = wdm_psd(testpsd, NF, NT, dt)

    testwdmgen = np.random.randn(*testwdmpsd.shape) * np.sqrt(testwdmpsd)

    if plots:
        plt.imshow(np.abs(testwdm.T), aspect='auto', origin='lower',extent=(t[0],t[-1],f[0],f[-1]), norm='log')
        plt.colorbar()
        plt.title("WDM transform of data (|w_nm|)")
        plt.xlabel("t [s]")
        plt.ylabel("f [Hz]")
        plt.savefig("wdm_data.png")
        plt.show()
        plt.imshow(np.sqrt(testwdmpsd_approx.T), aspect='auto', origin='lower',extent=(t[0],t[-1],f[0],f[-1]), norm='log')
        plt.colorbar()
        plt.title("Approximate wASD")
        plt.xlabel("t [s]")
        plt.ylabel("f [Hz]")
        plt.savefig("wpsd_approx.png")
        plt.show()
        plt.imshow(np.sqrt(testwdmpsd.T), aspect='auto', origin='lower',extent=(t[0],t[-1],f[0],f[-1]), norm='log')
        plt.colorbar()
        plt.title("Integrated wASD")
        plt.xlabel("t [s]")
        plt.ylabel("f [Hz]")
        plt.savefig("wpsd_exact.png")
        plt.show()
        plt.imshow(np.abs(np.sqrt(testwdmpsd.T) - np.sqrt(testwdmpsd_approx.T)), aspect='auto', origin='lower',extent=(t[0],t[-1],f[0],f[-1]), norm='linear')
        plt.colorbar()
        plt.title("Difference")
        plt.savefig("wpsd_diff.png")
        plt.show()
        plt.imshow(testwdm.T / np.sqrt(testwdmpsd_approx.T), aspect='auto', origin='lower',extent=(t[0],t[-1],f[0],f[-1]), norm='linear')
        plt.colorbar()
        plt.title("Whitened WDM data with approx wPSD")
        plt.savefig("whitened_approx.png")
        plt.show()
        plt.imshow(testwdm.T / np.sqrt(testwdmpsd.T), aspect='auto', origin='lower',extent=(t[0],t[-1],f[0],f[-1]), norm='linear')
        plt.colorbar()
        plt.title("Whitened WDM data with Neil wPSD")
        plt.savefig("whitened_exact_wdm.png")
        plt.show()


    from scipy.stats import chi2

    for label, cols in [("excluding DC layer", slice(1, None)), ("including DC layer", slice(None))]:
        whitened_approx = testwdm[:,cols] / np.sqrt(testwdmpsd_approx[:,cols])
        whitened_exact = testwdm[:,cols] / np.sqrt(testwdmpsd[:,cols])
        reference = np.random.randn(*whitened_exact.shape)

        n = whitened_exact.size
        lo = np.sqrt(chi2.ppf(0.025, n-1) / (n-1))
        hi = np.sqrt(chi2.ppf(0.975, n-1) / (n-1))

        print(f"--- {label} (n={n}) ---")
        print(f"95% CI for stdev of {n} unit normals: [{lo:.6f}, {hi:.6f}]")
        print(f"Stdev (approximate PSD): {np.std(whitened_approx):.6f}")
        print(f"Stdev (exact PSD):       {np.std(whitened_exact):.6f}")
        print(f"Stdev (reference N(0,1)):{np.std(reference):.6f}")
        print()
def inner_product_test(plots=True):
    WAVELET_DURATION = 7680
    NF = 1536
    NT = 338
    ND = NT*NF
    NFFT = ND/2 + 1
    Tobs = NT*WAVELET_DURATION
    dt = WAVELET_DURATION / NF
    DeltaF = 1/(2*dt*NF)
    DeltaT = WAVELET_DURATION
    # we need to generate this data in WDM actually
    t = np.arange(NT) * WAVELET_DURATION
    f = np.arange(NF) * DeltaF

    fcenter = f + DeltaF/2
    tcenter = t + DeltaT/2

    f_dft = np.arange(NFFT) * 1/Tobs
    testwdmpsd_approx = np.zeros((NT,NF))
    for i, ti in enumerate(tcenter):
        testwdmpsd_approx[i,:] = testpsdfunc(fcenter)/ (2*dt)

    testpsd = testpsdfunc(f_dft)
    psdfft = testpsd
    testwdmpsd_fullf = wdm_psd(testpsd, NF, NT, dt) 

    datawdm = np.random.randn(*testwdmpsd_fullf.shape) * np.sqrt(testwdmpsd_fullf)

    fsin = fcenter[NF//2]
    fullt = np.arange(0,NF*NT*dt, dt)
    signal = 1e-12*np.sin(2*np.pi*fsin*fullt)
    signalwdm = transform_wavelet_freq_time(signal, NF, NT)
    signalfft = np.fft.rfft(signal)

    datawdm += signalwdm
    datafft = inverse_wavelet_freq(datawdm, NF, NT)

    # fourier inner product: DFT*dt * DFT * dt / Sn * df
    # WDM inner product: ? Michael is using DWT*DWT / Sn_layer * df
    # in general, loglike = -0.5 * <d-h|d-h> + const
    # and SNR is <d|h>/sqrt(<h|h>)
    df = 1/Tobs
    inner_prod_wdm = np.sum(datawdm*signalwdm / testwdmpsd_fullf) / np.sum(np.sqrt(signalwdm*signalwdm / testwdmpsd_fullf))
    emin = int(fsin / df) - 100
    emax = int(fsin / df) + 100
    inner_prod_fft = np.sum(signalfft[emin:emax]*datafft[emin:emax].conj() / psdfft[emin:emax]) / np.sum(signalfft[emin:emax]*signalfft[emin:emax].conj() / psdfft[emin:emax])
    usual_inner_prod_fft = 4*df*np.sum((signalfft[emin:emax]*dt)*(datafft[emin:emax]*dt).conj() / psdfft[emin:emax]).real
    print("Inner products:")
    print(f"\tUsual formula FFT:\n\t{usual_inner_prod_fft}")
    print(f"\tMy formula FFT:\n\t{inner_prod_fft}")
    print(f"\tMy formula WDM:\n\t{inner_prod_wdm}")
    if plots:
        plt.imshow(np.abs(datawdm.T), aspect='auto', origin='lower',extent=(t[0],t[-1],f[0],f[-1]), norm='log')
        plt.colorbar()
        plt.title("WDM stationary data, noise+signal (|w_nm|)")
        plt.savefig("signal_wdm_data.png")
        plt.show()
        plt.loglog(f_dft[1:], np.abs(datafft[1:])) 
        plt.loglog(f_dft[1:], np.sqrt(testpsd[1:]*ND/(2*dt)))
        plt.title("FFT stationary data, noise+signal (ASD)")
        plt.savefig("signal_fft_data.png")
        plt.show()


if __name__ == '__main__':
    stat_test()
    nonstat_test()
    inner_product_test()

