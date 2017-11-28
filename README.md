# Thesis_Codes
This repo contains (some) codes for the back-calculation models discussed in the thesis "Estimating HIV incidence from multiple sources of data" by Francesco Brizzi (University of Cambridge).

The stan scripts are named using the following conventions:

"AI" denotes age-independent back-calculation models. These are of four types: 
- RW1AI.stan (random walk from an intermediate point of the epidemic)
- RW1978AI.stan (random walk from the beginning of the epidemic)
- GPAI.stan (Gaussian Process)
- splAI.stan (splines).

"ptens" denotes the age-dependent back-calculation model, using a tensor product spline to model incidence.
"tps" denotes the age-dependent back-calculation model, using a thin plate spline to model incidence.

The scripts with the following suffix
"quar"
"age_diag" and "age_diag1" denotes the m
