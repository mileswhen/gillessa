# Gillessa
Python implementation of **Gille**spie **s**tochastic-**s**imulation **a**lgorithm, for simulating [FKBP-F36V](https://www.rcsb.org/structure/1BL4) homodimerization in presence of [AP20187](https://www.sigmaaldrich.com/BE/en/product/sigma/sml2838). Contains several experimental notebooks for conducting basic reaction-diffusion simulations on a 2D surface, which corresponds to a typical [TIRF microscopy setup](https://www.leica-microsystems.com/science-lab/total-internal-reflection-fluorescence-tirf-microscopy/).

## Structure
* **SRDS.ipynb**: simple stochastic simulations of reactions. Contains tests of objects and classes for incorporating diffusion, as well as some prelim results.
* **SSAmodel**: complete code for generating 2D reaction-diffusion sims with plotting functions. Note these are compartment-based.
* **gillespie.py**: incomplete script for generating 2D reaction-diffusion sims
* **particle_sims.ipynb**:  particle-based simulations of homodimerization. Further checks to validate trajectory [correlation algorithm](https://www.researchgate.net/publication/256072559_Correlation_of_Dual_Colour_Single_Particle_Trajectories_for_Improved_Detection_and_Analysis_of_Interactions_in_Living_Cells) in single-particle tracking experiments.

## References
1. [Exact stochastic simulation of coupled chemical reactions](https://pubs.acs.org/doi/abs/10.1021/j100540a008)
2. [A practical guide to stochastic simulations of reaction-diffusion processes](https://arxiv.org/abs/0704.1908)
3. [GillespieSSA: Implementing the Stochastic
Simulation Algorithm in R](https://www.deenaschmidt.com/Teaching/Fa17/Gillespie-paper.pdf)
