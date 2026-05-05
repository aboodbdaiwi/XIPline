# decay-modeled-cs-paper
![Figure1](Fig1.png)

A repository that contains non-Cartesian reconstruction software for hyperpolarized xenon MRI. Written by Joseph Plummer (`joseph.plummer@nih.gov`, GitHub username: `joeyplum`). 

Please post all issues, questions, or collaboration requests on the Issues tab. Alternatively, feel free to contact me directly via email. 

## Installation:

In a Linux, WSL2, or Mac terminal, run the following commands in sequence to install (it is recommended to set up `ssh` keys (see GitHub help or ask ChatGPT)).

1. `git clone git@github.com:cchmc-cpir/decay-modeled-cs-paper.git`
2. `cd decay-modeled-cs-paper`
3. `conda update -n base -c defaults conda`
4. `make conda`
5. `conda activate decay-modeled-cs-paper`
6. `make pip`

**Troubleshooting**:

1. This repository was tested on an NVIDIA GPU. If running on a system without
   an NVIDIA GPU, please remove (comment out with `#` key) the following packages from `environment.yaml`:
   - `cudnn`
   - `nccl`
   - `cupy`
2. Additionally, if not using an NVIDIA GPU, please set `devnum = -1` for each
   reconstruction script. Otherwise, you will encounter a bug with `cupy`.

## Running the scripts: 

It is recommended to run all scripts using the `Run Current File in Interactive Window' tool in VScode so that a reconstructions can be monitored and figures can be easily viewed. However, the scripts also work in command line. 
1. `simulation_recon_2d.py`
2. `ventilation_recon_2d.py`
3. `ventilation_recon_3d.py`
4. `gasex_recon_3d.py`

## Optimizations:

Iterative reconstructions have many knobs that can be turned to speed up and improve image quality. There is a chance that the settings used in this code may not apply as efficiently to your own site's data. Therefore, I recommend changing the following settings:

1. `lamda` --> the regularization parameter for each regularization function.
2. `num_normal` --> the number of $A^H A$ operations (i.e. the number of iterations to solve the inner $||Ax-y||^2_2$ loop).
3. `num_iters` --> number of outer iterations in ADMM algorithm.
4. `rho` --> ADMM step size (sometimes, this can be increased to speed up convergence, although this directly multiplies `lamda`)
5. `scan_resolution`/`recon_resolution` --> tune for your own application. `resolution` = matrix size (Philips terminology--please don't hate!)

If you want to talk about the algorithms, ADMM, other optimization theory, or simply want help with getting the code running for your application, please post on the Issues tab.

## Uninstall:

To uninstall, run the following commands:

1. `conda activate`
2. `make clean`


## DOI:
https://doi.org/10.1002/mrm.30188
