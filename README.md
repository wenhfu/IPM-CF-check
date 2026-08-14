# Numerical Verification: IPM-CF-check

This repository contains the supplementary numerical verification script for the paper:

> **"A Primal--Dual Interior-Point Method for Nonsymmetric Conic Optimization with Conjugate-Free Scaling"**
> *by Rui-Jin Zhang, Wenhao Fu, and Yu-Hong Dai.*

The main script (`check_thm4_14_large.py`) rigorously verifies the scalar inequalities
and neighborhood bounds for both the predictor and corrector steps presented in the
paper's theoretical analysis, with the widened parameter choice

$$
\eta=0.016,\quad \beta_l=0.9,\quad \beta_u=0.95,\quad \omega=0.008,\quad \alpha_c=0.96,\quad \sigma=0.9058.
$$

The band $[\beta_l,\beta_u]=[0.9,0.95]$ is ten times wider than the original
$[0.9,0.905]$. As established via direct differentiation in the paper (Appendix B),
all relevant quantities attain their worst-case values at the endpoint $\nu=1$ or at
the limit $\nu\to\infty$. The script therefore

1. checks every scalar inequality at the worst-case endpoint $\nu=1$, together with
   the $\sigma$-range conditions at the limits $\lim_{\nu\to\infty}\beta_l^+(\nu)$ and
   $\lim_{\nu\to\infty}\beta_u^+(\nu)$; and
2. independently re-verifies all inequalities on a fine grid $\nu\in[1,10^6]$.

## Verification Overview

### 1. Predictor Step

The script computes the worst-case scalar quantities for the predictor point $z^+$ at
$\nu=1$ (e.g., step length $\alpha_p(1)$, $\eta^+$, $\beta_l^+(1)$, $\beta_u^+(1)$, and
$\delta^+$) and verifies

$$
z^+ \in \mathcal{N}(0.0405,\,0.8998,\,0.9501).
$$

### 2. Corrector Step

Using the fixed corrector step length $\alpha_c=0.96$, the script evaluates the
necessary conditions and bounds for the corrected point $z^{++}$ at $\nu=1$,
confirming that the iterate is pulled back into the original predefined neighborhood:

$$
z^{++} \in \mathcal{N}(0.0154,\,0.9006,\,0.9105)
      \subset \mathcal{N}(0.016,\,0.9,\,0.95).
$$

### 3. Complexity constants

The script also verifies the constants used in the complexity analysis:
$c_p>0.00633$, $c_\mu=0.00576$, and $\xi=c_p-c_\mu>0.00057$.

## Usage

The repository consists solely of the verification script. You can execute it directly
using Python (requires `numpy`):

```bash
python check_thm4_14_large.py
```

The legacy script `check_thm4-14.py` verifies the original narrow parameter choice
$(\eta,\beta_l,\beta_u,\omega,\alpha_c,\sigma)=(0.02,0.9,0.905,0.005,0.85,0.9025)$
and is kept for reference.

## Expected Output

<<<<<<< HEAD
Running the script will print the computed scalar quantities for both the predictor
and corrector steps, boolean flags verifying that all strict inequalities hold, the
reported neighborhood membership statements, and the result of the independent grid
check over $\nu\in[1,10^6]$.
=======
Running the script will print the computed scalar quantities for both the predictor and corrector steps, along with boolean flags verifying that all strict inequalities hold true. These numerical results are consistent with the calculations in the proof of Theorem 4.14 in Appendix B of the paper.
>>>>>>> 1c1b0f8dcedc59f96efc6e60bb5463171a894b74
