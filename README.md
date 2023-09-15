# This package contains the data and the computer code to reproduce the results in Zhang, Pan, Zhu, Fang and Wang's paper titled "A Latent Space Model for Weighted Keyword Co-occurrence Networks with Applications in Knowledge Discovery in Statistics".
Step 1: Load the following R-packages with the specified versions (The first three packages are used for parallel computing. If parallelism is not required, they do not need to be loaded.): 
	parallel version: 4.1.2
	foreach version: 1.5.2
	doParallel version: 1.0.16
	irlba version: 2.3.5
	tidyverse version: 1.3.1
	Rspectra version: 0.16.1

Step 2: Simulation Example 1 (Figure 3 (a) in paper): Use the ‘Simulation’ folder:
	Run the ‘simulation-N.R’ code

Step 3: Simulation Example 2 (Figure 3 (b) in paper): Use the ‘Simulation’ folder:
	Run the ‘simulation-d.R’ code

Step 4: Simulation Example 3 (Figure 4 (a) in paper): Use the ‘Simulation’ folder:
	Run the ‘simulation-N.R’ code

Step 5: Simulation Example 4 (Figure 4 (b) in paper): Use the ‘Simulation’ folder:
	Run the ‘simulation-T.R’ code

Step 6: Case Study Part: Use the ‘CaseStudy’ folder:
	Run the ‘real_data-model.R’ code
