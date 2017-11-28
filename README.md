# Thesis_Codes
This repo contains (some) codes for the back-calculation models discussed in the thesis "Estimating HIV incidence from multiple sources of data" by Francesco Brizzi (University of Cambridge).

Age-independent and age-dependent CD4-based multi-state Bayesian back-calculation models are written in the Stan language. This is interfaced with R using the rstan package.

# Preliminaries
Three RData files contain three template data-sets.
- TempDataAI.RData (age-independent quarterly data)
- TempDataAD.RData (age-dependent yearly data)
- TempDataADQt.RData (age-dependent quarterly data)

# R scripts
Two R scripts provide working examples:
- RunAI.R (runs the age-independent back-calculation model).
- RunAD.R (runs the age-dependent back-calculation model).

It is possible to fit different variants of the back-calculation models (e.g. using a yearly and a quarterly time scale, using different models for the infection process) by calling different stan scripts. For further details see the stan scripts section or the instructions in the R file.

# Stan scripts
The stan scripts are named using the following conventions:

"AI" denotes age-independent back-calculation models. These are of four types: 
- RW1AI.stan (random walk from an intermediate point of the epidemic)
- RW1978AI.stan (random walk from the beginning of the epidemic)
- GPAI.stan (Gaussian Process)
- splAI.stan (splines).

The other models refer to age-dependent back-calculation. The two simplest models are:
- "ptens.stan" denotes the age-dependent back-calculation model, using a tensor product spline to model incidence. 
- "tps.stan" denotes the age-dependent back-calculation model, using a thin plate spline to model incidence.

Both models consider age-independent diagnosis probabilities and a yearly time scale. However more complex models can be considered. These are described by the following conventions:
- "quar" the model uses a quarterly (rather than yearly) time scale
- "age_diag" and "age_diag1" use age-dependent (rather than age-independent) diagnosis probabilities. "age-diag" and "age-diag1" respectively refer to age and age-and-state specific intercept for the diagnosis process (see Section 8.4.2)

Note that a pdf of the thesis will be uploaded when corrections are approved by the University of Cambridge.
