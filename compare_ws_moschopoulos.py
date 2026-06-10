#!/usr/bin/env python3
"""
Compare the Welch-Satterthwaite and Moschopoulos approximations to the
distribution of a coarse-grained WDM power statistic.

Background
----------
Coarse-graining the WDM time-frequency grid averages Q neighbouring squared
wavelet coefficients into one coarse cell (draft2.tex, "Welch and Bartlett-like
coarse-graining"):

    P = (1/Q) sum_{i=0}^{Q-1} w_i^2 ,   w_i ~ N(0, S_i)

so the *unnormalised* statistic

    Y = Q * P = sum_i S_i * X_i ,       X_i ~ chi^2_1   (i.i.d.)

is a weighted sum of chi-squares with weights S_i (the per-fine-pixel noise
variance, which drifts across the cell for non-stationary noise). When the S_i
are unequal this has no closed form. Two ways to cope:

  * Welch-Satterthwaite : a single Gamma whose first two moments match Y.
                          h_eff = (sum S_i)^2 / sum S_i^2 effective dof.
                          This is exactly noise/src/glass_noise_model.c
                          precompute_coarse_Qeff().
  * Moschopoulos (1985) : the EXACT density of a sum of independent unequal
                          Gammas, as a convergent infinite mixture of Gammas.

This script overlays both against a Monte-Carlo histogram of Y for a PSD that
drifts across the coarse cell, and sweeps the drift to show where WS breaks down.
"""

import numpy as np
from scipy import stats
from scipy.special import gammaln


# --------------------------------------------------------------------------
# The two approximations
# --------------------------------------------------------------------------
def welch_satterthwaite_gamma(S):
    """Single moment-matched Gamma for Y = sum_i S_i X_i, X_i ~ chi^2_1.

    Returns a frozen scipy Gamma with shape h_eff/2 and scale 2*g, matching
    mean(Y) = sum S_i and var(Y) = 2 sum S_i^2.

    h_eff = (sum S_i)^2 / sum S_i^2   is the Welch-Satterthwaite effective dof,
    identical to Qeff = s1*s1/s2 in precompute_coarse_Qeff().
    """
    S = np.asarray(S, float)
    s1 = S.sum()
    s2 = (S * S).sum()
    h_eff = s1 * s1 / s2              # effective degrees of freedom (<= Q)
    g = s2 / s1                       # effective scale
    # Y ~ g * chi^2_{h_eff} = Gamma(shape=h_eff/2, scale=2g)
    return stats.gamma(a=h_eff / 2.0, scale=2.0 * g), h_eff


def moschopoulos_logpdf(y, alpha, theta, Kmax=400, tol=1e-12, nterms=None):
    """Exact log-density of Y = sum_i Gamma(shape=alpha_i, scale=theta_i).

    Moschopoulos (1985): with theta_min = min theta_i,

        rho      = sum_i alpha_i
        C        = prod_i (theta_min/theta_i)^{alpha_i}
        gamma_k  = (1/k) sum_i alpha_i (1 - theta_min/theta_i)^k        (k>=1)
        delta_0  = 1, delta_{k+1} = 1/(k+1) sum_{j=1}^{k+1} j gamma_j delta_{k+1-j}
        f(y)     = C sum_k delta_k * Gamma_pdf(y; shape=rho+k, scale=theta_min)

    Done in log space for stability. Returns log f(y) for array y.

    If `nterms` is given the series is truncated at exactly that many terms
    (k = 0 .. nterms-1) with no convergence short-circuit -- used to study how
    the truncated series approaches the exact density.
    """
    alpha = np.asarray(alpha, float)
    theta = np.asarray(theta, float)
    y = np.atleast_1d(np.asarray(y, float))

    if nterms is not None:
        Kmax = nterms - 1
        tol = 0.0                     # disable early break: use exactly nterms

    theta_min = theta.min()
    rho = alpha.sum()
    logC = np.sum(alpha * np.log(theta_min / theta))

    # gamma_k for k = 1..Kmax
    one_minus = 1.0 - theta_min / theta
    gamma = np.empty(Kmax + 1)
    gamma[0] = 0.0
    for k in range(1, Kmax + 1):
        gamma[k] = np.sum(alpha * one_minus**k) / k

    # delta_k recursion
    delta = np.zeros(Kmax + 1)
    delta[0] = 1.0
    for k in range(Kmax):
        acc = 0.0
        for j in range(1, k + 2):
            acc += j * gamma[j] * delta[k + 1 - j]
        delta[k + 1] = acc / (k + 1)

    # mixture of Gamma(shape=rho+k, scale=theta_min) densities, weights C*delta_k
    logpdf = np.full_like(y, -np.inf)
    logy = np.log(y)
    for k in range(Kmax + 1):
        dk = delta[k]
        if dk <= 0:
            continue                  # delta_k >= 0 for this problem; skip tiny/neg
        shape = rho + k
        # log Gamma_pdf(y; shape, theta_min)
        log_g = (shape - 1) * logy - y / theta_min - gammaln(shape) - shape * np.log(theta_min)
        term = logC + np.log(dk) + log_g
        logpdf = np.logaddexp(logpdf, term)
        # crude convergence check on the largest contribution
        if tol > 0 and k > 5 and np.exp(np.max(term) - np.max(logpdf)) < tol:
            break
    return logpdf if logpdf.size > 1 else logpdf[0]


def stationary_gamma(S):
    """Naive 'stationary-in-window' model from draft2.tex: pretend every fine
    pixel in the cell has the mean variance Sbar, so Y = Q*P ~ Sbar * chi^2_Q.

    Returns a frozen Gamma(shape=Q/2, scale=2*Sbar). Keeps the full dof Q
    (unlike Welch-Satterthwaite), so it is too narrow once the variance drifts.
    """
    S = np.asarray(S, float)
    Q = S.size
    Sbar = S[S.shape[0]//2]
    #Sbar = S.mean()
    return stats.gamma(a=Q / 2.0, scale=2.0 * Sbar)


def terms_to_converge(S, tv_tol=1e-3, Kref=600, ngrid=2000):
    """Number of Moschopoulos series terms needed for the truncated density to
    reach the (effectively exact) Kref-term reference, measured by total-
    variation distance TV = 0.5*int|p_K - p_ref| dy < tv_tol.

    Also returns the geometric ratio r = 1 - theta_min/theta_max that sets the
    asymptotic term decay, and rho*theta_max/theta_min, the rough term count the
    series must climb to before its base Gamma (anchored at theta_min) reaches
    the bulk of the mass.
    """
    S = np.asarray(S, float)
    Q = S.size
    alpha = np.full(Q, 0.5)
    theta = 2.0 * S
    rho = alpha.sum()
    r = 1.0 - theta.min() / theta.max()
    climb = rho * theta.max() / theta.min()   # ~ where base Gamma must reach

    # grid out to the far tail of the reference
    ref_dist = stats.gamma(a=rho, scale=theta.max())   # heaviest single component
    y = np.linspace(1e-6, ref_dist.ppf(0.9995) * 1.3, ngrid)
    dy = y[1] - y[0]
    p_ref = np.exp(moschopoulos_logpdf(y, alpha, theta, nterms=Kref))

    K = 1
    while K <= Kref:
        p_K = np.exp(moschopoulos_logpdf(y, alpha, theta, nterms=K))
        tv = 0.5 * np.sum(np.abs(p_K - p_ref)) * dy
        if tv < tv_tol:
            break
        K = K + 1 if K < 8 else int(np.ceil(K * 1.25))   # coarsen the search
    return K, r, climb


# --------------------------------------------------------------------------
# Build a drifting-PSD coarse cell and Monte-Carlo the true Y
# --------------------------------------------------------------------------
def make_cell(Q, drift, S0=1.0):
    """Per-fine-pixel variances S_i across one coarse cell of Q pixels.

    drift = (S_max - S_min)/S_mean is the fractional variance spread across the
    cell.  drift=0 -> stationary (S_i all equal); larger drift -> WS degrades.
    """
    x = np.linspace(-0.5, 0.5, Q)
    S = S0 * (1.0 + drift * x)        # linear drift in variance across the cell
    return S


def monte_carlo_Y(S, n=400000, rng=None):
    rng = rng or np.random.default_rng(0)
    # Y = sum_i S_i * X_i, X_i ~ chi^2_1  (i.e. w_i = sqrt(S_i) * N(0,1), Y=sum w_i^2)
    w = rng.standard_normal((n, S.size)) * np.sqrt(S)
    return (w * w).sum(axis=1)


def plot_likelihood_shape(Q=16, drift=1.5, n_real=400000, n_rug=40,
                          conv_terms=(1, 2, 3, 6, 25, 50, 100, 200, 300), rng=None,
                          out="ws_moschopoulos_likelihood_shape.png"):
    """Single figure, two panels, at one in-cell drift:

      left : the measured likelihood shape (histogram over many noise
             realizations of Y) vs the stationary, Welch-Satterthwaite, and
             (converged) Moschopoulos model densities;
      right: the Moschopoulos series truncated at increasing numbers of terms,
             on a log-y axis so tail convergence is visible.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rng = rng or np.random.default_rng(7)
    S = make_cell(Q, drift)
    alpha = np.full(Q, 0.5)
    theta = 2.0 * S

    # "several noise realizations" -> sampling distribution of the statistic Y.
    y = monte_carlo_Y(S, n=n_real, rng=rng)
    y_rug = rng.choice(y, size=n_rug, replace=False)   # a few individual draws

    edges = np.linspace(0.0, np.percentile(y, 99.7), 140)
    centers = 0.5 * (edges[:-1] + edges[1:])
    p_true, _ = np.histogram(y, bins=edges, density=True)

    stat = stationary_gamma(S)
    ws, h_eff = welch_satterthwaite_gamma(S)
    p_stat = stat.pdf(centers)
    p_ws = ws.pdf(centers)
    p_mosch = np.exp(moschopoulos_logpdf(centers, alpha, theta))  # auto-converged

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

    # ---- left: method comparison ----
    axL.fill_between(centers, p_true, alpha=0.3, color="0.5", step="mid",
                     label=f"measured ({n_real:,} realizations)")
    axL.plot(centers, p_stat, lw=2, color="tab:green",
             label=r"stationary  $\Gamma(Q/2,\,2\bar S)$, dof$=Q$=%d" % Q)
    axL.plot(centers, p_ws, lw=2, color="tab:orange",
             label=r"Welch-Satterthwaite ($h_{\rm eff}$=%.1f)" % h_eff)
    axL.plot(centers, p_mosch, lw=2, ls="--", color="tab:blue",
             label="Moschopoulos (converged)")
    axL.plot(y_rug, np.full_like(y_rug, -0.002 * p_true.max()), "|",
             color="0.3", ms=10, label="individual realizations")
    axL.set_title(f"likelihood shape, Q={Q}, in-cell drift={drift:.0%}")
    axL.set_xlabel(r"$Y = Q\,P = \sum_i S_i\,\chi^2_{1,i}$")
    axL.set_ylabel("density  (=likelihood vs. data)")
    axL.legend(fontsize=8)

    # ---- right: Moschopoulos convergence in number of terms ----
    axR.fill_between(centers, p_true, alpha=0.3, color="0.5", step="mid",
                     label="measured")
    cmap = plt.cm.viridis(np.linspace(0.15, 0.9, len(conv_terms)))
    for c, K in zip(cmap, conv_terms):
        p_K = np.exp(moschopoulos_logpdf(centers, alpha, theta, nterms=K))
        axR.plot(centers, p_K, lw=1.8, color=c, label=f"{K} term" + ("s" if K > 1 else ""))
    axR.set_yscale("log")
    axR.set_ylim(p_true[p_true > 0].min() * 0.3, p_true.max() * 2)
    axR.set_title("Moschopoulos truncation convergence")
    axR.set_xlabel(r"$Y$")
    axR.set_ylabel("density (log)")
    axR.legend(fontsize=8, title="series terms")

    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"saved likelihood-shape plot -> {out}")


def analyze_convergence(Q=16, out="moschopoulos_convergence.png"):
    """Why does Moschopoulos need many terms at large drift, and how does it
    behave at moderate (10-20%) drift?  Table + plot of terms-to-converge vs the
    in-cell variance spread, against the analytic scaling.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    drifts = [0.05, 0.10, 0.15, 0.20, 0.30, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8]
    print(f"\nMoschopoulos convergence vs in-cell drift (Q={Q}):")
    print(f"{'drift':>6} {'Smin/Smax':>10} {'r=1-tmin/tmax':>14} "
          f"{'h_eff':>7} {'climb~rho*Smax/Smin':>20} {'terms@TV<1e-3':>14}")

    rows = []
    for drift in drifts:
        S = make_cell(Q, drift)
        _, h_eff = welch_satterthwaite_gamma(S)
        K, r, climb = terms_to_converge(S)
        smin_smax = S.min() / S.max()
        rows.append((drift, smin_smax, r, h_eff, climb, K))
        print(f"{drift:6.2f} {smin_smax:10.3f} {r:14.3f} {h_eff:7.1f} "
              f"{climb:20.1f} {K:14d}")

    drifts_a = np.array([r[0] for r in rows])
    climb_a = np.array([r[4] for r in rows])
    K_a = np.array([r[5] for r in rows])

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(drifts_a, K_a, "o-", color="tab:blue", label="measured (TV < 1e-3)")
    ax.plot(drifts_a, climb_a, "s--", color="tab:gray",
            label=r"analytic climb $\rho\,S_{\max}/S_{\min}$")
    ax.axvspan(0.10, 0.20, color="tab:green", alpha=0.12,
               label="moderate drift (10-20%)")
    ax.set_xlabel("in-cell variance drift")
    ax.set_ylabel("Moschopoulos terms to converge")
    ax.set_yscale("log")
    ax.set_title(f"Series term count is set by $S_{{\\min}}$ (Q={Q})")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"saved convergence plot -> {out}")


def check_midpoint_vs_average(modfile="modulation.dat", Q=169, NT=8112,
                              Nlayer_gal=44, n_real=200, rng=None,
                              out="midpoint_vs_average_check.png"):
    """Offline check of the C fix (galaxy_modulation_cache_for_Q: midpoint ->
    cell-average) using the real galaxy modulation in `modfile`.

    The coarse data statistic P_q = (1/Q) sum_i w_i^2 has mean <S>_cell, the
    cell-average of the per-pixel variance. The OLD model sampled the modulation
    at the cell midpoint m(mid); the NEW model uses <m>_cell. Treating the band as
    pure galaxy (amplitude 1) and fitting a single galactic amplitude A across all
    galaxy-dominated cells, the per-cell MAP gives

        A_hat = mean_{layer,cell} P_{lq} / m_model_q ,
        E[A_hat] = (1/Ncoarse) sum_q <m>_q / m_model_q ,

    which is exactly 1 for the cell-average model and biased by the curvature for
    the midpoint model. This quantifies the residual Q=169-vs-Q=1 amplitude shift
    the fix removes, against the (tiny) statistical error of the fit.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        d = np.loadtxt(modfile)
    except Exception as e:
        print(f"(midpoint-vs-average check skipped: cannot read {modfile}: {e})")
        return

    rng = rng or np.random.default_rng(11)
    Ncoarse = NT // Q
    # col0 is the orbital angle alpha; the real run spans ~2 revolutions (~2 yr).
    # NOTE: running noise_wavelet_mcmc rewrites modulation.dat in the CWD for
    # *its* duration, so a short test run will clobber the real file -- warn if so.
    revs = (d[-1, 0] - d[0, 0]) / (2 * np.pi)
    print(f"\n(modulation.dat spans {revs:.2f} orbital revolutions; "
          f"the real 2-yr run is ~1.97)")
    if revs < 1.0:
        print("  WARNING: this looks like a short/clobbered run, not the 2-yr "
              "grid; regenerate with the real --duration before trusting numbers.")
    tfrac = np.linspace(0, 1, len(d))
    tfine = np.linspace(0, 1, NT)
    m = np.interp(tfine, tfrac, d[:, 1])           # XX modulation, per fine pixel (=S, amp 1)

    cells = m[:Ncoarse * Q].reshape(Ncoarse, Q)
    m_avg = cells.mean(axis=1)                       # <m>_cell  -> the NEW (fixed) model
    mid = Q // 2
    m_mid = 0.5 * (cells[:, mid - 1] + cells[:, mid])  # midpoint  -> the OLD model
    cell_bias = (m_mid - m_avg) / m_avg              # signed per-cell model error (old)

    # analytic expected recovered amplitude (infinite data per cell)
    Ahat_mid_exp = np.mean(m_avg / m_mid)
    Ahat_avg_exp = np.mean(m_avg / m_avg)            # == 1 exactly

    # Monte-Carlo scatter of A_hat over the galaxy-dominated layers (same m(t) in
    # every layer; only the realization differs). Pure-galaxy, single channel.
    Ah_mid, Ah_avg = [], []
    for _ in range(n_real):
        w = rng.standard_normal((Nlayer_gal, Ncoarse * Q)) * np.sqrt(m[:Ncoarse * Q])
        P = (w * w).reshape(Nlayer_gal, Ncoarse, Q).mean(axis=2)   # (Nlayer, Ncoarse)
        Ah_mid.append((P / m_mid[None, :]).mean())
        Ah_avg.append((P / m_avg[None, :]).mean())
    Ah_mid = np.array(Ah_mid); Ah_avg = np.array(Ah_avg)
    sig = Ah_avg.std()

    print(f"\nMidpoint-vs-average offline check (Q={Q}, Ncoarse={Ncoarse}, "
          f"{Nlayer_gal} galaxy-dominated layers, XX):")
    print(f"  per-cell model error (m_mid-<m>)/<m>: median={np.median(np.abs(cell_bias)):.2%}"
          f"  worst={np.max(np.abs(cell_bias)):.2%}")
    print(f"  recovered galactic amplitude A_hat (true=1):")
    print(f"    OLD midpoint model : E[A]={Ahat_mid_exp:.5f}  (bias {Ahat_mid_exp-1:+.2%}), "
          f"MC {Ah_mid.mean():.5f} +/- {Ah_mid.std():.5f}")
    print(f"    NEW average  model : E[A]={Ahat_avg_exp:.5f}  (bias {Ahat_avg_exp-1:+.2%}), "
          f"MC {Ah_avg.mean():.5f} +/- {Ah_avg.std():.5f}")
    print(f"  statistical sigma(A)~{sig:.5f}; midpoint bias is "
          f"{abs(Ahat_mid_exp-1)/sig:.1f} sigma -> coherent, does not average down")

    # ----- plot -----
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))
    tdays = np.arange(Ncoarse) * Q * 7680 / 86400.0

    axL.axhline(0, color="0.7", lw=1)
    axL.bar(tdays, 100 * cell_bias, width=Q * 7680 / 86400.0 * 0.9,
            color=np.where(cell_bias < 0, "tab:red", "tab:blue"), alpha=0.7,
            label="old model error per cell")
    axL.set_xlabel("cell time [days]")
    axL.set_ylabel(r"$(m_{\rm mid}-\langle m\rangle)/\langle m\rangle$  [%]")
    axt = axL.twinx()
    axt.plot(tdays, m_avg, "k-", lw=1.5, label=r"$\langle m\rangle_{\rm cell}$")
    axt.set_ylabel("galaxy modulation (mean 1)")
    axL.set_title("midpoint sampling error vs annual modulation")
    axL.legend(loc="upper left", fontsize=8); axt.legend(loc="upper right", fontsize=8)

    axR.hist(Ah_mid, bins=40, alpha=0.6, color="tab:orange",
             label=f"OLD midpoint  ({Ah_mid.mean():.4f})")
    axR.hist(Ah_avg, bins=40, alpha=0.6, color="tab:blue",
             label=f"NEW average   ({Ah_avg.mean():.4f})")
    axR.axvline(1.0, color="k", ls="--", lw=1.5, label="true A = 1")
    axR.set_xlabel(r"recovered galactic amplitude $\hat A$")
    axR.set_ylabel(f"MC realizations ({n_real})")
    axR.set_title("amplitude recovery: midpoint biased, average unbiased")
    axR.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"saved midpoint-vs-average check -> {out}")


# --------------------------------------------------------------------------
# Main: numbers + plot
# --------------------------------------------------------------------------
def main():
    Q = 16
    rng = np.random.default_rng(1)

    print(f"Coarse cell with Q = {Q} fine pixels; alpha_i = 1/2, theta_i = 2 S_i\n")
    print(f"{'drift':>6} {'h_eff (WS dof)':>16} {'mean true/WS/Mosch':>30} "
          f"{'KL(true||WS)':>14} {'KL(true||Mosch)':>16}")

    # drift is the fractional variance spread (S_max-S_min)/S_mean across the
    # cell; must stay < 2 so every S_i remains positive.
    drifts = [0.0, 0.6, 1.2, 1.8]
    results = {}
    for drift in drifts:
        S = make_cell(Q, drift)
        alpha = np.full(Q, 0.5)
        theta = 2.0 * S

        y = monte_carlo_Y(S, n=400000, rng=rng)
        ws, h_eff = welch_satterthwaite_gamma(S)

        # histogram-based KL of the MC truth against each model density
        lo, hi = 0.0, np.percentile(y, 99.9)
        edges = np.linspace(lo, hi, 120)
        centers = 0.5 * (edges[:-1] + edges[1:])
        p_true, _ = np.histogram(y, bins=edges, density=True)

        p_ws = ws.pdf(centers)
        p_mo = np.exp(moschopoulos_logpdf(centers, alpha, theta))

        def kl(p, q):
            m = (p > 0) & (q > 0)
            dx = centers[1] - centers[0]
            return np.sum(p[m] * np.log(p[m] / q[m])) * dx

        kl_ws = kl(p_true, p_ws)
        kl_mo = kl(p_true, p_mo)

        print(f"{drift:6.1f} {h_eff:16.2f} "
              f"{y.mean():9.3f}/{ws.mean():.3f}/{np.trapezoid(centers*p_mo, centers):.3f}  "
              f"{kl_ws:14.4e} {kl_mo:16.4e}")

        results[drift] = (centers, p_true, p_ws, p_mo, h_eff)

    # ----- plot -----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=False)
        for ax, drift in zip(axes.ravel(), drifts):
            centers, p_true, p_ws, p_mo, h_eff = results[drift]
            ax.fill_between(centers, p_true, step=None, alpha=0.25,
                            color="0.5", label="Monte-Carlo (truth)")
            ax.plot(centers, p_ws, lw=2, color="tab:orange",
                    label=f"Welch-Satterthwaite ($h_{{eff}}$={h_eff:.1f})")
            ax.plot(centers, p_mo, lw=2, ls="--", color="tab:blue",
                    label="Moschopoulos (exact)")
            ax.set_title(f"in-cell variance drift = {drift:.0%}")
            ax.set_xlabel(r"$Y = Q\,P = \sum_i S_i\chi^2_{1,i}$")
            ax.set_ylabel("density")
            ax.legend(fontsize=8)
        fig.suptitle(f"Coarse-grained WDM power statistic, Q={Q} fine pixels")
        fig.tight_layout()
        out = "ws_vs_moschopoulos.png"
        fig.savefig(out, dpi=130)
        print(f"\nsaved plot -> {out}")
    except Exception as e:
        print(f"\n(plot skipped: {e})")

    # requested figure: measured likelihood shape vs stationary / WS /
    # Moschopoulos, with Moschopoulos truncation-convergence panel
    try:
        plot_likelihood_shape(Q=Q, drift=1.5, rng=rng)
        # moderate-drift regime (the one the paper actually targets): WS is
        # essentially exact and Moschopoulos converges in a handful of terms.
        plot_likelihood_shape(Q=Q, drift=0.15, rng=rng,
                              conv_terms=(1, 2, 3, 5, 8, 12),
                              out="ws_moschopoulos_likelihood_shape_moderate.png")
        analyze_convergence(Q=Q)
        # offline check of the C midpoint->cell-average fix on the real modulation
        check_midpoint_vs_average(rng=rng)
    except Exception as e:
        print(f"(likelihood-shape plot skipped: {e})")


if __name__ == "__main__":
    main()
