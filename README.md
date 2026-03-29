# Parameter Estimation for a Stochastic Delay-Differential Equation (SDDE)

Simulation-based parameter estimation for stochastic delay-differential equations (SDDEs), recovering drift and volatility from noisy sample paths.

---

## Overview
This project demonstrates simulation-based recovery of drift and volatility parameters for a stochastic delay-differential equation (SDDE) from noisy sample paths.

The objective is to infer latent dynamics from stochastic time-series data and validate the inferred model via forward simulation.

We consider the SDDE
\[
dX(t) = \bigl(\mu_0 X(t) + \mu_1 X(t-\tau)\bigr)\,\mathrm{d}t + \sigma\,\mathrm{d}W(t),
\]
with known delay \( \tau \) and unknown parameters \( \mu_0, \mu_1, \sigma \).

---

## Method

### Data generation
Multiple independent sample paths are simulated using an Euler–Maruyama discretisation with a fixed delay.

### Volatility estimation
The volatility parameter \( \sigma \) is estimated via quadratic variation of the process increments.

### Drift estimation
The drift coefficients \( (\mu_0, \mu_1) \) are estimated using least-squares regression on the discretised dynamics.

### Validation
The learned parameters are validated by simulating the corresponding deterministic delayed system and comparing:
- simulated stochastic trajectories,
- deterministic dynamics using true parameters,
- deterministic dynamics using learned parameters,
- uncertainty bands induced by the estimated volatility.

---

## Results
The inferred parameters closely match the ground truth across multiple stochastic trials.  
Forward simulations using the learned model reproduce the qualitative behaviour of the original system, including delayed feedback effects and stochastic spread.

The figure below shows parameter recovery and model validation.

---

## Parameter recovery and model validation

> **Note:** GitHub does not render EPS or PDF files inline.  
> The figure below should be viewed by downloading the PDF.

[Download parameter-estimation figure](SDDE_Parameter_Estimation_100-eps-converted-to.pdf)

---

## Why This Matters
This experiment mirrors a common quantitative research workflow:
- infer latent model structure from noisy data,
- separate signal from stochastic variability,
- stress-test learned dynamics via simulation.

The same principles apply broadly to problems involving noisy time-series data, delayed response, and model validation under uncertainty.

---

## Implementation Notes
- **Language:** MATLAB  
- **Numerical method:** Euler–Maruyama  
- **Focus:** robustness and interpretability rather than optimisation or machine learning

---

## Files
- `SDDE_Parameter_Estimation.m` — main simulation, inference, and validation script  
- `SDDE_Parameter_Estimation_100-eps-converted-to.pdf` — parameter recovery and validation figure
