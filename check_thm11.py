"""
Numerical verification for Theorem (guarantee neighborhood) -- large-band set.

This script verifies the scalar inequalities and neighborhood bounds for the
predictor--corrector analysis of the paper with the widened parameter choice

    eta    = 0.016
    beta_l = 0.9
    beta_u = 0.95
    omega  = 0.008
    alpha_c = 0.96
    sigma  = 0.9058

The band [beta_l, beta_u] = [0.9, 0.95] is ten times wider than the original
[0.9, 0.905].  All scalar inequalities are checked at the worst case nu = 1
(monotonicity in nu is established in the paper's appendix by direct
differentiation) and, independently, on a fine grid nu in [1, 1e6].
"""

import numpy as np

# ----------------------------------------------------------------------
#  parameters
# ----------------------------------------------------------------------
eta     = 0.016
beta_l  = 0.9
beta_u  = 0.95
omega   = 0.008
alpha_c = 0.96
sigma   = 0.9058

# uniform lower bound for beta_l^+(nu) over all nu >= 1 (reported in the paper)
beta_hat_l = 0.8998

# ----------------------------------------------------------------------
#  derived constants
# ----------------------------------------------------------------------
delta = eta**2 + 2*eta
q = omega**2 * (1 - delta)

def gamma(nu):
    return np.sqrt((eta + np.sqrt(nu))**2 / (1 - delta) + beta_u / 2.0)

def alpha_p(nu):
    return omega * np.sqrt(1 - delta) / gamma(nu)

def beta_l_p(nu):
    a = alpha_p(nu)
    return ((1 - a) * beta_l - q) / (1 - a + q)

def beta_u_p(nu):
    a = alpha_p(nu)
    return ((1 - a + a**2 / 4.0) * beta_u) / (1 - a - q)

eta_p = (omega / (1 - omega)**2
         + (eta + (1 - delta) * (omega + omega**2) + np.sqrt(1 - delta**2) * omega)
           / ((1 - omega - (1 - delta) * omega**2) * (1 - omega)))
delta_p = eta_p**2 + 2 * eta_p

def h(v):
    return (v - sigma)**2 / (4.0 * v)

def theta(nu):
    return max(h(beta_l_p(nu)), h(beta_u_p(nu)))

def gamma_p(nu):
    return np.sqrt(eta_p**2 / (1 - delta_p) + 2 * theta(nu))

def omega_p(nu):
    return alpha_c * gamma_p(nu) / np.sqrt(1 - delta_p)

def eta_tilde(nu):
    return np.sqrt(1 + delta_p) * gamma_p(nu) / (1 - omega_p(nu))

def eta_pp(nu):
    c = (alpha_c * gamma_p(nu))**2
    w = omega_p(nu)
    return (1 / (1 - c)) * (
        (w / (1 - w))**2
        + (1 - alpha_c + delta_p) * w / (alpha_c * (1 - w))
        + (1 - alpha_c) * eta_tilde(nu)
        + c)

def beta_l_pp(nu):
    c = (alpha_c * gamma_p(nu))**2
    return ((1 - alpha_c) * beta_l_p(nu) + sigma * alpha_c - c) / (1 + c / nu)

def beta_u_pp(nu):
    c = (alpha_c * gamma_p(nu))**2
    return ((1 - alpha_c) * beta_u_p(nu) + sigma * alpha_c + alpha_c**2 * theta(nu)) \
           / (1 - c / nu)

# ----------------------------------------------------------------------
#  worst-case values at nu = 1
# ----------------------------------------------------------------------
a1    = alpha_p(1.0)
bl1   = beta_l_p(1.0)
bu1   = beta_u_p(1.0)
t1    = theta(1.0)
gp1   = gamma_p(1.0)
wp1   = omega_p(1.0)
c1    = (alpha_c * gp1)**2
et1   = eta_tilde(1.0)
epp1  = eta_pp(1.0)
blpp1 = beta_l_pp(1.0)
bupp1 = beta_u_pp(1.0)

bl_inf = (beta_l - q) / (1 + q)
bu_inf = beta_u / (1 - q)

alpha_c_bound = (np.sqrt((sigma - bl1)**2 + 4 * bl1 * gp1**2) + sigma - bl1) \
                / (2 * gp1**2)

c_p  = omega * (1 - delta) / np.sqrt((1 + eta)**2 + (1 - delta) * beta_u / 2.0)
c_mu = alpha_c * (sigma - bl1)
xi   = c_p - c_mu

# paper's conservative complexity constants (Theorem 4.14 proof):
#   c_mu := alpha_c*(sigma - beta_hat_l^+) = 0.00576,  xi := c_p - c_mu > 0.00057
c_mu_paper = alpha_c * (sigma - beta_hat_l)
xi_paper   = c_p - c_mu_paper

# ----------------------------------------------------------------------
#  checks
# ----------------------------------------------------------------------
checks = {
    "0 < eta < sqrt(2)-1          ": 0 < eta < np.sqrt(2) - 1,
    "0 < beta_l <= beta_u         ": 0 < beta_l <= beta_u,
    "0 < omega < 1                ": 0 < omega < 1,
    "0 < alpha_c <= 1             ": 0 < alpha_c <= 1,
    "1-omega-(1-delta)omega^2 > 0 ": 1 - omega - (1 - delta) * omega**2 > 0,
    "pred. step cond. (nu=1)      ": (1 - a1) * beta_l - q > 0,
    "eta^+ < sqrt(2)-1, delta^+<1 ": (eta_p < np.sqrt(2) - 1) and (delta_p < 1),
    "sigma > lim beta_l^+         ": sigma > bl_inf,
    "sigma < lim beta_u^+         ": sigma < bu_inf,
    "omega^+(1) < 1               ": wp1 < 1,
    "corrector step cond. (nu=1)  ": (1 - alpha_c) * bl1 + sigma * alpha_c - c1 > 0,
    "beta_l^++(1) > beta_l        ": blpp1 > beta_l,
    "beta_u^++(1) < beta_u        ": bupp1 < beta_u,
    "eta^++(1) < eta              ": epp1 < eta,
    "xi = c_p - c_mu > 0          ": xi > 0,
    # complexity constants reported in the proof of Theorem 4.14
    "beta_hat_l^+ uniform lower bound": bl1 > beta_hat_l,
    "c_p > 0.00633                ": c_p > 0.00633,
    "c_mu = alpha_c(sigma-beta_hat) = 0.00576": abs(c_mu_paper - 0.00576) < 1e-12,
    "xi = c_p - c_mu > 0.00057    ": xi_paper > 0.00057,
}


print("=" * 66)
print("Fixed constants (large-band set)")
print(f"  eta={eta}  beta_l={beta_l}  beta_u={beta_u}  omega={omega}")
print(f"  alpha_c={alpha_c}  sigma={sigma}   (band width {beta_u-beta_l:.2f})")
print("=" * 66)
for k, v in checks.items():
    print(f"  [{'OK' if v else 'FAIL'}] {k}")
print()

print(f"delta is : {delta:.12f}")
print(f"q is : {q:.12e}")
print(f"1-omega-(1-delta)*omega^2 is : {1 - omega - (1 - delta) * omega**2:.12f}")
print(f"gamma(1) is : {gamma(1.0):.12f}")
print(f"alpha_p(1) is : {a1:.12f}")
print(f"alpha_p upper bound is (>alpha_p(1)): {a1:.12f}")
print(f"(1-alpha_p(1))*beta_l-q is : {(1 - a1) * beta_l - q:.12f}")
print(f"eta^+ is : {eta_p:.12f}")
print(f"d beta_l^+ is : {-(1 + beta_l) * q / (1 - a1 + q)**2:.12e}")
print(f"d beta_u^+ is : {beta_u * (2 - a1) * (a1 + 2 * q) / (4 * (1 - a1 - q)**2):.12e}")
print(f"beta_l^+(1) is : {bl1:.12f}")
print(f"beta_u^+(1) is : {bu1:.12f}")
print(f"lim beta_l^+ is : {bl_inf:.12f}")
print(f"lim beta_u^+ is : {bu_inf:.12f}")
print(f"c_p is : {c_p:.12f}")
print(f"c_mu is : {c_mu:.12f}")
print(f"xi := c_p - c_mu is : {xi:.12f}")
print(f"beta_hat_l^+ is : {beta_hat_l}")
print(f"c_mu (paper, using beta_hat_l^+) is : {c_mu_paper:.12f}")
print(f"xi (paper) is : {xi_paper:.12f}")
print()
print(f"delta^+ is : {delta_p:.12f}")
print(f"theta(1) is : {t1:.12e}")
print(f"gamma^+(1) is : {gp1:.12f}")
print(f"omega^+(1) is : {wp1:.12f}")
print(f"alpha_c upper bound is (>alpha_c): {alpha_c_bound:.12f}")
print(f"eta_tilde(1) is : {et1:.12f}")
print(f"eta^++(1) is : {epp1:.12f}")
print(f"corrector step LHS is : {(1 - alpha_c) * bl1 + sigma * alpha_c - c1:.12f}")
print(f"beta_l^++(1) is : {blpp1:.12f}")
print(f"beta_u^++(1) is : {bupp1:.12f}")
print()


# check: z^+ in the reported intermediate neighborhood
if eta_p < 0.0405 and bl1 > 0.8998 and bu1 < 0.9501:
    print("z^+  lies in N(0.0405, 0.8998, 0.9501)")
else:
    print("z^+  does not lie in N(0.0405, 0.8998, 0.9501)")
# check: z^++ in the reported intermediate neighborhood
if epp1 < 0.0154 and blpp1 > 0.9006 and bupp1 < 0.9105:
    print("z^++ lies in N(0.0154, 0.9006, 0.9105)")
else:
    print("z^++ does not lie in N(0.0154, 0.9006, 0.9105)")
# check: z := z^++ lies in the original neighborhood
if epp1 < eta and blpp1 > beta_l and bupp1 < beta_u:
    print(f"z    lies in N({eta:.3f}, {beta_l:.1f}, {beta_u:.2f})")
else:
    print(f"z    does not lie in N({eta:.3f}, {beta_l:.1f}, {beta_u:.2f})")
print()

# ----------------------------------------------------------------------
#  independent grid check over nu in [1, 1e6]
# ----------------------------------------------------------------------
nus = [1.0] + list(np.logspace(0, 6, 500))
grid_ok = True
worst = {}
for nu in nus:
    b = {
        "pred step": (1 - alpha_p(nu)) * beta_l - q,
        "sigma_l": sigma - beta_l_p(nu),
        "sigma_u": beta_u_p(nu) - sigma,
        "omega_p": 1 - omega_p(nu),
        "corr step": (1 - alpha_c) * beta_l_p(nu) + sigma * alpha_c
                     - (alpha_c * gamma_p(nu))**2,
        "ret_l": beta_l_pp(nu) - beta_l,
        "ret_u": beta_u - beta_u_pp(nu),
        "eta_ret": eta - eta_pp(nu),
    }
    grid_ok = grid_ok and all(v > 0 for v in b.values())
    for k, v in b.items():
        worst[k] = min(worst.get(k, np.inf), v)

print("Grid check over nu in [1, 1e6]: " + ("ALL OK" if grid_ok else "FAILED"))
for k, v in sorted(worst.items()):
    print(f"    {k:10s} min margin = {v:.6e}")
print()

print("All checks passed:", all(checks.values()) and grid_ok)
