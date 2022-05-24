# Gene Regulatory Network with Developmental Memory Capacity
This repository includes scripts that implement experiments from "The Evolution of Phenotypic Correlations and Developmental Memory" by Watson et al. (https://onlinelibrary.wiley.com/doi/full/10.1111/evo.12337). One can replicate experiments 1 and 2 depicted in figures 1 and 2 of the original paper.

## Reproduction

Each script replicates a specific aspect of experiments 1 and 2 while helper functions generate the respective plots. Namely:

- `single_env.py`: serves as a placeholder for the gene regulatory network (GRN) and implements experiment 1 with a single phenotypic target.

    <img src="readme_source\single_env.png" width="500" height="300">

    **Figure 1**: Single Environment Selection Pressure

- `varying_env.py`: implements experiment 2 with a varying phenotypic target.

    <img src="readme_source\two_target_env.png" width="800" height="300">

    **Figure 2**: Two-target Environment Selection Pressure


- `hebb_net_exp2.py`: serves as the implementation of evolution via the Hebb rule: $\Delta w = x_{in}\cdot y $
- `generate_plots.py`: includes helper functions that generate the plots shown in figure 1 & 2 of the original paper.
- `varying_env_maskB.py`: serves as an extension of the experiments presented in the original paper and explores the capacity of GRNs when one gradually prevents genes from interacting.

    <img src="readme_source\extension.png" width="900" height="150">

    **Figure 3**: Masking the GRN



GPL-3.0-only
