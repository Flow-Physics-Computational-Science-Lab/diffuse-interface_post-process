# Volume correction for ACDI

In the paper, they address the challenge of accurately calculating droplet/bubble
properties (e.g., volume, number) in diffuse-interface two-phase flow simulations.
Currently, flood-fill algorithms can truncate a significant portion of the volume
of droplets/bubbles contained within the diffuse interface region or artificially
merge multiple droplets/bubbles. This error is also dependent on the volume fraction
cutoff value, which is typically chosen to be 0.5 arbitrarily, in the flood-fill
algorithms. They propose a simple volume-correction approach that incorporates an
analytical approximation of the truncated volume to correct for the missing
droplet/bubble volumes. The proposed method results in accurately recovering the
dispersed phase volumes with minimal volume error over a wide range of volume
fraction cutoff values, and hence, also accurately recovering the number of
droplets/bubbles. It can be a valuable tool for accurate calculation of drop/bubble
size distributions for analysis and for Eulerian-to-Lagrangian conversion of the
dispersed phase in multi-scale modeling approaches.

# TODO:

1. Generalize input;

# Cite

```bib
@article{nathan:2025,
    title = {Accurate calculation of bubble and droplet properties in diffuse-interface two-phase simulations},
    journal = {Journal of Computational Physics},
    volume = {538},
    pages = {114190},
    year = {2025},
    issn = {0021-9991},
    doi = {https://doi.org/10.1016/j.jcp.2025.114190},
    author = {Pranav J. Nathan and Suhas S. Jain},}
```
