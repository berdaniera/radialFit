# radialFit
This code can be used to fit new observations of radial sap flux to generate new profiles, as a supplement for Berdanier et al. (2016) Tree Physiology. 

Currently, it is set up to accommodate trees in different groups. It generates parameter predictions for each group for the best fitting model from our analysis, the gamma model (see Berdanier et al. 2016). Fitting is done with a Bayesian regression of the natural log of the sap flux observations on the absolute depth into the xylem. Sampling from posterior distributions is done with gibbs sampling. Estimates of profile parameters are done with a Metropolis step and the uncertainty and individual random effects are sampled directly from conjugate distributions. Right now the priors are weak or uninformative. For extensions, you could use our fitted parameters as priors, which would require modifying the code in the gmaUpdate() function. 

Contact me (Aaron Berdanier, aaron.berdanier@gmail.com) if you have any questions. I will add additional documentation if requested.
