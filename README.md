This repository contains the code (written in R) used to create Figure one in the manuscript "Strategies for determining heteroaggregation attachment efficiencies of engineered nanoparticles in aquatic environments", by Praetorius et al. submitted to ES Nano in 2019. The code was written by Arnaud Clavier and Antonia Praetorius.

The code consists only of 1 file: r_script_to_create_Fig1
To recreate the 4 plots in Figure 1, run the code as it is in R

If using this code to create material for publication, please reference the paper listed above. 

This code is used to evaluate the relative importance of the 3 concurrent mechanisms (perikinetic, orthokinetic and differential sedimentation) contributing to the overall collision rate between 2 aggregating particles of different sizes and densities. This process is very relevant for heteroaggregation of particles in the nano- or micrometer size range.  
The size and density of particle i are fixed while size and density of particle j vary within a specified range. The fixed values and the ranges can be varied to explore different scenarios. The default values provide in the code represent the parameters used for creating Figure 1 as represented int he paper. They are chosen to represent realistic cases a metal oxide type engineered nanoparticle (i) heteroaggregating with a natural suspended matter (SPM) paticle (j).
The collision rate constant for heteroaggregation is calculated according to Equation 2 in the paper following classical colloid theory. This represents the most simple case of rectilinear collision mechanisms and Stokes Law is used to calculate settling
In theory, this code can be extended to include curvilinear (or other) correction terms and/or to use different settling equations, 
but these options are not implemented in the current code. For examples of such correction terms, see references in relevant sectiosn in the paper. 
Part 1 of the code calculates the collision individual rate constants and saves them in the file: kcollspecific.dat
Part 2 plots the indivudal rate constants as well as their sum in contour plots. 
Part 3 (commented out) is optional for saving the plots and requires installation of orca: https://github.com/plotly/orca. To use this part, install orca and uncomment code. Alternatively, other saving methods can be used to create different figures files. 

