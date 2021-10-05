===================
aHCAMB
===================
<div>
	<h2>
		Background
	</h2>
	This is a modified version of CAMB for the axi-Higgs model and the standard axion cosmology, written based on CAMB (version June 2021) and axionCAMB. <br>

	The axi-Higgs model is proposed recently in <a href='https://arxiv.org/abs/2102.11257'> arxiv: 2102.11257</a>, in which the Higgs is coupled with an ultra-light axion. With the linear perturbative analysis, it has been shown in <a href='https://arxiv.org/abs/2105.01631'> arxiv: 2105.01631</a> that it can address many tensions in cosmology, including:
	<ul>
		<li>Hubble tension</li>
		<li>S8 tension</li>
		<li>Primordial lithium problem</li>
		<li>Cosmic birefringence</li>
	</ul>
	Complementary to the perturbative analysis, I present a modification of CAMB to calculate the CMB spectrum and the matter power spectrum in the axi-higgs model, which shows a reasonable agreement with the previous analysis while uncovering more cosmological information from the full spectra.
</div>

<div>
	<h2> Compile and Run </h2>
    A <code> make </code> file is provided, which would automatically handle the linkage of CAMB and Recfast++ code and compile them all at once. <br>

	Compile the code with make: <br>
	<code>
		cd ./fortran && make clean && make
	</code> <br>
	Run	the code with an input file: <br>
	<code>
		cd ./fortran && ./camb ahCAMB_input.ini
	</code><br>
	Example of inputting cosmological parameters, with explanations can be found in "./fortran/aHCAMB_input.ini".
    <br><br>
	<h4> Dependencies </h4>
	The code is based on CAMB and axionCAMB, which requires a fortran compiler to compile. <br>
	In addition, a modified Recfast++ code is used to implement the recombination physics. To compile, a C++ compiler is required. The GSL library for scientific computation in C should be installed properly; otherwise, an error message: <code> /usr/bin/ld: cannot find -lgsl, -lgslcblas </code> may be raised.
    <br>
    The required codebases CAMB and Recfast++ have already been included in this repository , there is no need to install them separately. 
</div>

<div>
	<h2> Remark </h2>
	Python wrapper has not been modified in this version. 
	Note that the default recombination mode is Recfast++ to account for the variation of the electron mass.
	Axion physics are included in "AxionStandard.f90" and "AxionHiggs.f90". The changes in source files are marked with the keyword "Nhan".
</div>

<div> 
	<h2> Citation </h2>
	When using this code, please consider citing: ahCAMB (Luu 2021, in preparation), and
	<ul>
		<li> The original CAMB and axionCAMB paper at: 
			<a href='https://arxiv.org/abs/astro-ph/9911177'> arXiv:astro-ph/9911177 </a>
			, <a href='https://arxiv.org/abs/1201.3654'> arXiv:1201.3654 </a>
			, <a href='https://arxiv.org/abs/1410.2896'> arXiv:1410.2896  </a>
		</li>
		<li> The original Recfast and Recfast++ paper at: 
			<a href='https://arxiv.org/abs/astro-ph/9909275'> arXiv:astro-ph/9909275 </a>
			, <a href='https://arxiv.org/abs/1003.4928'> arXiv:1003.4928 </a>
			, <a href='https://arxiv.org/abs/1010.3631'> arXiv:1010.3631 </a>
		</li>
	</ul>
</div>
