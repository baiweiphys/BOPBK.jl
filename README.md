# BOPBK.jl

**Version 1.0** (Released on 2026-05-28)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20434445.svg)](https://doi.org/10.5281/zenodo.20434445)

[GitHub Repository](https://github.com/baiweiphys/BOPBK.jl/) | [Releases](https://github.com/baiweiphys/BOPBK/releases)


**`BOPBK.jl`** (BO–Product-Bi-Kappa) is a Julia numerical code developed by Wei Bai, under the supervision of Dr. Huasheng Xie, who also proposed the BO framework. It is designed for analyzing plasma waves and instabilities in both space and laboratory plasmas, specializing in obliquely propagating waves in magnetized, multi-species hot plasmas.


The code supports a wide range of multi-component velocity distribution functions, such as:
- Anisotropic drift loss-cone product-bi-Kappa (PBK)
- Anisotropic drift loss-cone Kappa-Maxwellian (KM)
- Anisotropic drift loss-cone bi-Maxwellian (BM)
- Hybrid combinations of these distributions


## Developer & Supervisor
**Wei Bai** (Developer)  
College of Electrical and Power Engineering  
Taiyuan University of Technology  
Taiyuan 030024, China  
Email: baiweiphys@gmail.com  

**Huasheng Xie** (Supervisor – proposed the BO framework)  
ENN Science and Technology Development Co., Ltd., Langfang 065001, China  
& VeloAlpha Technology Co., Ltd., Beijing 100080, China  
Email: huashengxie@gmail.com  


## Program Description

**BO-PBK** is a comprehensive code for simultaneously solving all significant roots of the dispersion relation for obliquely propagating waves in magnetized multi-species plasmas.
The code supports various velocity distribution functions, including:
- Anisotropic drift loss-cone product-bi-Kappa (Type I PBK)
- Anisotropic drift loss-cone kappa-Maxwellian (Type I KM)
- Anisotropic drift loss-cone bi-Maxwellian (BM)
- Hybrid combinations of these distributions

Benchmark tests confirm that the solver accurately reproduces standard kinetic results and efficiently resolves waves and instabilities. BO-PBK provides a robust and computationally efficient tool for wave and stability analysis in complex plasma systems relevant to both space and laboratory applications.

## Input Parameters
The background magnetic field is denoted as $B_0$ (in Tesla), the oblique propagation angle as $\theta$ (in degrees), $J$ as a non-negative integer for $J$-pole expansion, and `nk` as the number of grid points.

The `bopbk.in` file requires the following parameters for each species $s$:
- $q_s$: charge (in elementary charge units)
- $m_s$: mass (relative to electron mass)
- $n_s$: number density (in m⁻³)
- $T_{\parallel s}$ and $T_{\perp s}$: parallel and perpendicular temperatures (in eV)
- $u_{s0}/c$: drift velocity (normalized by speed of light)
- $\sigma_s$: loss-cone factor
- $\kappa_{\parallel s}$ and $\kappa_{\perp s}$: parallel and perpendicular kappa indices
- $\kappa_{\parallel s,\text{th}}$: threshold value for $\kappa_{\parallel s}$
- $N_s$: Bessel function truncation order for species $s$


## Benchmark Examples
Five benchmark cases are provided as representative examples in `runall.ipynb`:
1. R-, L-, and P-mode waves
2. Whistler instability
3. Firehose instability
4. EMIC waves
5. Beam-mode

## When referencing `BO-PBK` in a publication, please cite:
```bibtex
@article{BAI2026110281,
title = {BO-PBK: A new solver for dispersion relations of obliquely propagating waves in multi-species plasmas with anisotropic loss-cone drift product-bi-kappa distributions},
journal = {Computer Physics Communications},
volume = {327},
pages = {110281},
year = {2026},
issn = {0010-4655},
doi = {https://doi.org/10.1016/j.cpc.2026.110281},
url = {https://www.sciencedirect.com/science/article/pii/S0010465526002638},
author = {Wei Bai and Huasheng Xie},
}
```

[Computer Physics Communications 327 (2026) 110281]: https://doi.org/10.1016/j.cpc.2026.110281
[CPC Library link to program files]: https://doi.org/10.17632/zy3nhcr7xf.1


## License
This project is licensed under the BSD 3-Clause License — see the [LICENSE](LICENSE) file for details.

---

Created by Wei Bai on 2026-05-25 
