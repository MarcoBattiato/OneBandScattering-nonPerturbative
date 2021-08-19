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


Nx=3;     % Set the number of required points in the mesh along X and Y direction
Ny=Nx;
Lx=Nx;     % Set the width of domain in X and Y direction. Mesh element width along x =Lx/Nx=1 in this case
Ly=Nx;
disp=@(x,y) 0.06*sqrt((x-0.5*Lx)^2+(y-0.5*Ly)^2);    % set the required dispersion.
b=Band(Nx,Ny,Lx,Ly);                    % define a band object and get energy coefficients and Energy maximum and minimum in each element
b=b.Energycoefficients(disp);
b=b.Energylimits();

pBS=Plot(100,Lx,Ly);   % This will help to visualize the band structure
pBS=pBS.GetVal(b);
figure();                       
pBS.Plotpop();

 tic
  List=ListScatTensor(b);  % This part calculates the list and writes the List files in the folder. Please ensure that ./Listfiles/ folder exists in Work Directory
 toc
 % Please ensure TimestepStatus.txt and AMatrix.txt are in Work Directory
 tic
startelem=1; endelem=Nx*Ny;
 Scattering=ScatTensorList(b,startelem,endelem); % This will calculate the scattering tensor and store it in the folder ./ScatTensfiles/. Please ensure this folder is present in Work Directory
 toc
 
%% %%%%%%%%%%%%%%% Plot the Scattering Rates to check the generated Scattering Tensor and to get a primary insight into the dynamics%%%%%%%%%%
 
one=@(x,y) 1;p=Plot(100,Lx,Ly);
Ones=b.Energycoefficients(one);
mu=0.1; Temp=700;
f=@(x,y) 1/(1+exp((disp(x,y)-mu)/(8.617*10^-5*Temp)));  % This will introduce a Fermi-Dirac population
EqPopulation=b.Energycoefficients(f);
p=p.GetVal(EqPopulation);
figure();                       
p.Plotpop();              % Plot Equilibrium distribution
Population=EqPopulation;
Lambda=Scattering.ScatteringRates(Scattering.C,Ones.Ecoeff,Scattering.IndicesList,Scattering.N,EqPopulation.Ecoeff);
Population.Ecoeff=Lambda;
p=p.GetVal(Population);
figure();                       
p.Plotpop();         % Plot Scattering rates
