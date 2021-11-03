# aHCAMB

## Background
This is a modified version of CAMB for the axi-Higgs model and the standard axion cosmology, written based on CAMB (version June 2021) and axionCAMB.

The axi-Higgs model is proposed recently in [arxiv: 2102.11257](https://arxiv.org/abs/2102.11257), in which the Higgs is coupled with an ultra-light axion. With a perturbative analysis, it has been shown in [arxiv: 2105.01631](https://arxiv.org/abs/2105.01631) that it can address many tensions in cosmology, including:

* Hubble tension
* S8 tension
* Primordial lithium problem
* Cosmic birefringence

Complementary to the perturbative analysis, I present a modification of CAMB to calculate the CMB spectrum and the matter power spectrum in the axi-higgs model, which shows a reasonable agreement with the previous perturbative analysis, meanwhile uncovering more cosmological information from the full spectra. 


## Compile and Run
A ``` make ``` file is provided, which would automatically handle the linkage of CAMB and Recfast++ code and compile them all at once. 

Compile the code with make: 
```bash
cd ./fortran && make clean && make
```
Run	the code with an input file: 
```bash
cd ./fortran && ./camb ahCAMB_input.ini
```
Example of inputting cosmological parameters, with explanations can be found in "./fortran/aHCAMB_input.ini".
    
### Dependencies
The code is based on CAMB and axionCAMB, which requires a fortran compiler to compile.
In addition, a modified Recfast++ code is used to implement the recombination physics. To compile, a C++ compiler is required. The GSL library for scientific computation in C should be installed properly; otherwise, an error message: `` /usr/bin/ld: cannot find -lgsl, -lgslcblas `` may be raised.

The required codebases CAMB and Recfast++ have already been included in this repository, there is no need to install them separately. 

## Remark
Python wrapper has not been modified in this version. 
Note that the default recombination mode is Recfast++ to account for the variation of the electron mass.
Axion physics are included in "AxionStandard.f90" and "AxionHiggs.f90".
The changes in source files are marked with the keyword "Nhan".

## Citation
When using this code, please consider citing:
* The original CAMB and axionCAMB papers: 
	* [arXiv:astro-ph/9911177](https://arxiv.org/abs/astro-ph/9911177)
	* [arXiv:1201.3654](https://arxiv.org/abs/1201.3654)
	* [arXiv:1410.2896](https://arxiv.org/abs/1410.2896')
* The original Recfast and Recfast++ papers: 
	* [arXiv:astro-ph/9909275](https://arxiv.org/abs/astro-ph/9909275)
	* [arXiv:1003.4928](https://arxiv.org/abs/1003.4928)
	* [arXiv:1010.3631](https://arxiv.org/abs/1010.3631)
* The original aHCAMB paper (if you use the axi-Higgs mode):
	* [arXiv:2111.01347](https://arxiv.org/abs/2111.01347)
