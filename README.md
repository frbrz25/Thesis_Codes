# Thesis_Codes
This repo contains (some) codes for the back-calculation models discussed in the thesis "Estimating HIV incidence from multiple sources of data" by Francesco Brizzi (University of Cambridge).

The stan scripts are named using the following conventions:

"AI" denotes age-independent back-calculation models. These are of four types: 
- RW1AI.stan (random walk from an intermediate point of the epidemic)
- RW1978AI.stan (random walk from the beginning of the epidemic)
- GPAI.stan (Gaussian Process)
- splAI.stan (splines).

The other models refer to the age-dependent back-calculation models. The two simplest models are:
- "ptens.stan" denotes the age-dependent back-calculation model, using a tensor product spline to model incidence. 
- "tps.stan" denotes the age-dependent back-calculation model, using a thin plate spline to model incidence.
Both models consider age-independent diagnosis probabilities and a yearly time scale.

More complex models are described by the following conventions
- "quar" the model uses a quarterly (rather than yearly) time scale
- "age_diag" and "age_diag1" use age-dependent (rather than age-independent) diagnosis probabilities. "age-diag" and "age-diag1" respectively refer to age and age-and-state specific intercept for the diagnosis process (see Section 8.4.2)

Note that a pdf of the thesis will be uploaded when corrections are approved by the university of Cambridge.
