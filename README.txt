%% Code   : OneBandScattering-nonPerturbative
%% Authors: Indrajit Wadgaonkar, Marco Battiato
%% Date   : 17 July 2021
%%
%% Matlab >>non-perturbative<< implementation of single band Boltzmann scattering
%%
%% This code is a Matlab test implementation of the second order variation [1] of the algorithm 
%% introduced in [2] and extended in [3]
%%
%% [1] I. Wadgaonkar, M. Wais, and M. Battiato, 
%% Numerical Solver for the out-of-equilibrium time dependent Boltzmann Collision operator: Application to 2D materials
%% under review
%% [2] M. Wais, K. Held, M. Battiato, 
%% Numerical solver for the time-dependent far-from-equilibrium Boltzmann equation, 
%% Comput. Phys. Commun. 264, 107877  (2021) 
%% [3] I. Wadgaonkar, R. Jain, M. Battiato, 
%% Numerical scheme for the far-out-of-equilibrium time-dependent Boltzmann collision operator: 1D second-degree momentum discretisation and adaptive time stepping, 
%% Comput. Phys. Commun 263, 107863 (2021).
%%
%% If used, please cite the work above



The following code is a Matlab implementation of the numerical algorithm proposed in the article.
It does not provide full functionalities and it is specialised for a single band case (even though the algorithm is not limited to that).


The Code is run in two main steps:
A) Generating the Scattering Tensor elements (at this stage scattering rates can be calculated) and 
B) Contracting this generated Scattering Tensor to propagate the initialised population distribution in time. 

It is assumed that all the files of the code are in the Working Directory of Matlab and that 2 folders named 'Listfiles' and 'ScatTensfiles' have also been previously created by the user in the Working Directory. 
File named 'parfor_progress.txt' will be automatically created and it can be used to monitor the status of the code when it is running.

A) Generating the Scattering Tensor elements:

- Open the file 'GenerateListAndScatteringTensor.m' with your Matlab Installation. 
- Set the required mesh precision using the appropriate value for Nx. Note that the scaling of the computational cost is ~N^2.5 where N is the total number of elements in the mesh ( Nx*Ny). Use a suitable value of N to keep the computational cost within practical limits.
- Set the appropriate functional form of the dispersion. Although this code requires a functional form of the dispersion it is equally straightforward to use any arbitrary dispersion. The future versions will allow for it.
- Run the Code. The code will first generate, in the folder './Listfiles', the list files for possible combinations of elements which will contribute to the calculation of scattering tensor. At this stage the Scattering Tensor is not constructed yet.
- Since the list calculated in the previous step can be very large and computational cost can be high, the user can decide to split the calculation of the scattering tensor over different runs. This is achieved by setting the values of startelem and endelem. For eg. If there are a total of 400 elements in the mesh, the user may opt to calculate the scattering Tensor for only elements between 1 to 200 first and then for the remaining later. Note that the user need not calculate the list again for this. The code generates the files for the elements of the scattering tensor in the folder './ScatTensfiles' in the Working directory. 
- Once the scattering tensor has been generated, the user can calculate the scattering rates to get an insight into the thermalisation dynamics. The equilibrium population can be chosen based on the values of chemical potential (mu) and the temperature (Temp).


B) Contraction of the generated Scattering Tensor to propagate the initialised population distribution in time:

- Open the file 'AddExcitationAndTimePropagate.m' with your Matlab Installation. 
- The first part of the code, setting mesh precision and dispersion, is a repetition of the previous step. If the objects are already in the Workspace this part of the code is redundant. However, it is recommended to use it to ensure that the correct variables are set for time propagation.
- The user needs to set the relevant chemical potential, mu and temperature, Temp for the equilibrium population. Note that excitations will be added to this equilibrium population later.
- The code allows the user to study either a localised excitation or an excitation in the form of a gaussian distribution. The user can set the appropriate variables eg. Radius of the Gaussian distribution, the amplitude, the coordinates for the localised excitation etc. 
- Set the required initial and final time for time propagation and run the code.
- The remaining section of the code helps you plot the results and also visualise them in the movie format.