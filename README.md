# Random-generation-of-biological-networks
Some scripts I have developed for the random generation of enzymatic networks (scale-free and random topologies) to explore changes in their responsiveness under different parameter perturbations.
Note that random networks can be generated with scale-free topology or with a random topology. 
Furthermore network responsiveness as judged by peaks and troughs can be quantified in terms of the absolute values reached by
the peaks/troughs or their difference from basal levels.
The uploaded scripts reflect this (although there is no fold-change analysis of networks with random topologies - but the codes can be easily adapted for this).
When generated models are analysed through the 'RandNetGen_analysis_.m' file, the definition of network response (ie. absolute or relative values of peaks/troughs)
must be the same as that defined in the 'RandomNetworkGeneration_.m' file that was run to generate the models.
