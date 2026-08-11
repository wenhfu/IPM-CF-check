
# Numerical Verification: IPM-CF-check

This repository contains the supplementary numerical verification script for the paper:

> **"A Primal--Dual Interior-Point Method for Nonsymmetric Conic Optimization with Conjugate-Free Scaling"**
> *by Rui-Jin Zhang, Wenhao Fu, and Yu-Hong Dai.*

The single Python script provided here (`check_thm4-14.py`) rigorously verifies the scalar inequalities and neighborhood bounds for both the predictor and corrector steps presented in the paper's theoretical analysis.

## Verification Overview

Our theoretical analysis defines a specific neighborhood $\mathcal{N}(\eta,\beta_l,\beta_u)$ controlled by the following fixed constants:

$$
\eta=0.02,\quad \beta_l=0.9,\quad \beta_u=0.905,\quad \omega=0.005,\quad \alpha_c=0.85,\quad \sigma=0.9025
$$

As established via direct differentiation in the paper, all relevant quantities attain their worst-case values at the endpoint $\nu = 1$. The script focuses on evaluating the exact results at this critical endpoint.

### 1. Predictor Step

The script computes the worst-case scalar quantities for the predictor point $z^+$ at $\nu=1$ (e.g., step length $\alpha_p(1)$, $\eta^+$, $\beta_l^+(1)$, $\beta_u^+(1)$, and $\delta^+$) and verifies:

$$
z^+ \in \mathcal{N}(0.03518,0.89995,0.90503)
$$

### 2. Corrector Step

Using the fixed corrector step length $\alpha_c = 0.85$, the script evaluates the necessary conditions and bounds for the corrected point $z^{++}$ at $\nu=1$, confirming that the iterate is pulled back into the original predefined neighborhood:

$$
z^{++} \in \mathcal{N}(0.01665,0.90028,0.90376) \subset \mathcal{N}(0.02,0.9,0.905)
$$

## Usage

The repository consists solely of the verification script. You can execute it directly using Python (requires `numpy`):

```bash
python check_thm4-14.py
```

## Expected Output

Running the script will print the computed scalar quantities for both the predictor and corrector steps, along with boolean flags verifying that all strict inequalities hold true. These numerical results are consistent with the calculations in the proof of Theorem 4.14 in Appendix B of the paper.
