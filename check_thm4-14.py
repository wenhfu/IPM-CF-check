import numpy as np

eta = 0.02
beta_l = 0.90
beta_u = 0.905
omega = 0.005
alpha_c = 0.85
sigma = 0.9025
delta = eta**2 + 2*eta

eta_p_fixed = 0.03518
beta_l_p_fixed = 0.89995
beta_u_p_fied = 0.90503
eta_pp_fixed = 0.01665
beta_l_pp_fixed = 0.90028
beta_u_pp_fied = 0.90376

def fun_gamma(nu):
    return np.sqrt((eta+np.sqrt(nu))**2/(1-delta) + beta_u/2)
def fun_alpha_p(nu):
    return omega*np.sqrt(1-delta)/fun_gamma(nu)

verbose = 1


nu = 1.0
alpha_p = fun_alpha_p(nu)
gamma = fun_gamma(nu)
omega_bound = 1 - omega - (1 - delta) * omega**2
alpha_p_bound = ( np.sqrt(beta_l**2+4*beta_l*gamma**2) - (beta_l) ) / (2*gamma**2)
eta_p = omega / (omega - 1)**2 + (eta + (1 - delta) * (omega + omega**2) + np.sqrt(1 - delta**2) * omega) / ((1 - omega - (1 - delta) * omega**2) * (1 - omega))
beta_l_p = ( (1-alpha_p) * beta_l - alpha_p**2*gamma**2 ) / (1 - alpha_p + gamma**2 * alpha_p**2)
beta_u_p = ( (1 - alpha_p + (alpha_p**2)/4) * beta_u ) / (1 - alpha_p - gamma**2 * alpha_p**2)
delta_p = eta_p**2 + 2 * eta_p

# \[
#     c_p:=\frac{\omega(1-\delta)}{\sqrt{(1+\eta)^2+(1-\delta)\beta_u/2}},
#     \quad
#     c_\mu:=\alpha_c(\sigma-\widehat\beta_l^+),
# \]
c_p = omega * np.sqrt(1 - delta) / fun_gamma(1.0)
c_mu = alpha_c * (sigma - beta_l_p)

delta_p = eta_p**2 + 2 * eta_p
theta = max( ( (beta_l_p - sigma)**2) / (4 * beta_l_p), ( (beta_u_p - sigma)**2) / (4 * beta_u_p) )
gamma_p = np.sqrt( (eta_p**2)/(1 - delta_p) + 2 * theta )
omega_p = ( alpha_c * gamma_p ) / ( np.sqrt(1 - delta_p) )
alpha_c_bound = ( np.sqrt( (sigma - beta_l_p)**2 + 4 * beta_l_p * gamma_p**2) +sigma - beta_l_p) / (2 * gamma_p**2)
eta_tilde = ( np.sqrt(1 + delta_p) * gamma_p ) / ( 1 - omega_p )
eta_pp = (1 - (alpha_c * gamma_p)**2)**(-1) * ( ( (omega_p / (1 - omega_p) )**2 ) + ( (1 - alpha_c + delta_p) * omega_p ) / ( alpha_c * (1 - omega_p) ) + (1 - alpha_c) * eta_tilde + (alpha_c * gamma_p)**2 )
beta_l_pp = ( (1 - alpha_c) * beta_l_p - (alpha_c*gamma_p)**2 + sigma * alpha_c ) * ( 1 + (alpha_c * gamma_p)**2 / nu )**(-1)
beta_u_pp = ( (1 - alpha_c) * beta_u_p + (alpha_c)**2 * theta + sigma * alpha_c ) * ( 1 - (alpha_c * gamma_p)**2 / nu )**(-1)


if verbose == 1:
    print(f"nu is : {nu:.5f}")
    print(f"delta is : {delta:.12f}")
    print(f"omega upper bound is (<1): {omega_bound:.12f}")
    print(f"gamma(1) is : {gamma:.12f}")
    print(f"alpha_p(1) is : {alpha_p:.12f}")
    print(f"alpha_p upper bound is (>alpha_p): {alpha_p_bound:.12f}")
    print(f"eta^+ is : {eta_p:.12f}")
    print(f"beta_l^+ is : {beta_l_p:.12f}")
    print(f"beta_u^+ is : {beta_u_p:.12f}")
    print(f"c_p is : {c_p:.12f}")
    print(f"c_mu is : {c_mu:.12f}")
    print(f"xi:=c_p-c_mu is : {c_p - c_mu:.12f}")
    print()
    print(f"delta^+ is : {delta_p:.12f}")
    print(f"theta is : {theta:.12e}")
    print(f"gamma^+ is : {gamma_p:.12f}")
    print(f"omega^+ is : {omega_p:.12f}")
    print(f"alpha_c upper bound is (>1): {alpha_c_bound:.12f}")
    print(f"eta_tilde is : {eta_tilde:.12f}")
    print(f"eta^++ is : {eta_pp:.12f}")
    print(f"beta_l^++ is : {beta_l_pp:.12f}")
    print(f"beta_u^++ is : {beta_u_pp:.12f}")
    print()
elif verbose == 2:
    print(f"nu is : {nu:.5f}")
    # check: z_+ in N
    if eta_p < eta_p_fixed and beta_l_p > beta_l_p_fixed and beta_u_p < beta_u_p_fied:
        print(f"z^+  lies in N({eta_p_fixed:.5f},{beta_l_p_fixed:.5f},{beta_u_p_fied:.5f})")
    else:
        print(f"z^+  does not lie in N({eta_p_fixed:.5f},{beta_l_p_fixed:.5f},{beta_u_p_fied:.5f})")

    # check: z_++ in N
    if eta_pp < eta_pp_fixed and beta_l_pp > beta_l_pp_fixed and beta_u_pp < beta_u_pp_fied:
        print(f"z^++ lies in N({eta_pp_fixed:.5f},{beta_l_pp_fixed:.5f},{beta_u_pp_fied:.5f})")
    else:
        print(f"z^++ does not lie in N({eta_pp_fixed:.5f},{beta_l_pp_fixed:.5f},{beta_u_pp_fied:.5f})")

    # check: z:=z_++ in the original N
    if eta_pp < eta and beta_l_pp > beta_l and beta_u_pp < beta_u:
        print(f"z    lies in N({eta:.5f},{beta_l:.5f},{beta_u:.5f})")
    else:
        print(f"z    does not lie in N({eta:.5f},{beta_l:.5f},{beta_u:.5f})")
    
    print()


