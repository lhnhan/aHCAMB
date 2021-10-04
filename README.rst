===================
aHCAMB
===================

The modified version of CAMB for the axi-Higgs model and the standard axion cosmology, written based on CAMB (version June 2021) and axionCAMB.
Python wrapper has not been modified in this version.
Compile the code with fortran: "cd ./fortran && make clean && make".
Run the code with an input file: "cd ./fortran && ./camb input.ini".
Example with explanations can be found in "./fortran/aHCAMB_input.ini".
Note that the default recombination mode is Recfast++ to account for the variation of the electron mass.

When using this code, please consider citing:
- The original CAMB and axionCAMB paper at: arXiv:astro-ph/9911177, arXiv:1201.3654, arXiv:1410.2896
- The original Recfast and Recfast++ paper at: arXiv:astro-ph/9909275, arXiv:1003.4928, arXiv:1010.3631
