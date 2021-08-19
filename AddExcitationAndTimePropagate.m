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

Nx=20;     % Set the number of required points in the mesh along X and Y direction
Ny=Nx;
Lx=Nx;     % Set the width of domain in X and Y direction. Mesh element width along x =Lx/Nx=1 in this case
Ly=Nx;
disp=@(x,y) 0.06*sqrt((x-0.5*Lx)^2+(y-0.5*Ly)^2);    % set the required dispersion.
b=Band(Nx,Ny,Lx,Ly);                    % define a band object and get energy coefficients and Energy maximum and minimum in each element
b=b.Energycoefficients(disp);

pBS=Plot(100,Lx,Ly);   % This will help to visualize if the correct band structure is being used
pBS=pBS.GetVal(b);
figure();                       
pBS.Plotpop();
%% %%%%%%%%%%%%%%%%%%% Introduce Equilibrium Population distribution and Stabilise it %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=Plot(100,Lx,Ly);
mu=0.1; Temp=700;
 f=@(x,y) 1/(1+exp((disp(x,y)-mu)/(8.617*10^-5*Temp)));  % This will introduce a Fermi-Dirac population and visualize it
Population=b.Energycoefficients(f); 
p=p.GetVal(Population);
figure();                       
p.Plotpop();
 
Scattering.atol=10^-6; Scattering.rtol=10^-6;dense_time=1.0;  % Set the tolerances for Time Propagation
init_dt=0.001; final_time=0.005;init_time=0.0;   % Set starting time step, final and initial time
Scattering=Scattering.RK853(init_dt,final_time,Population,init_time,dense_time); % Propagate the population to get numerical steady state distribution
save('Stabilized FD.mat','Population','-v7.3','-nocompression');  % Save the stabilised population distribution. 

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Add Excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('Stabilized FD.mat'); % Comment this if variables are already in workspace

%% %%%%%%%%%%%%% This introduces a Gaussian excitation
 Radius=6.0;
 Amplitude=0.2;   % 0.2 for Big Excitation and 0.1 for small excitation
 Spread=1.5; % 1.5 for big excitation and 1 for small excitation
 Excfunc=@(x,y) Amplitude*exp(-(sqrt((x-0.5*Lx)^2+(y-0.5*Ly)^2)-Radius)^2/(2*Spread^2));
 Excitation=b.Energycoefficients(Excfunc);
 Population.Ecoeff=Population.Ecoeff+Excitation.Ecoeff;

 %% %%%%%%%%%%%%%%%%%%%%%%%%% This part introduces localised excitation. Comment out accordingly to study one type of excitation%%%%%%%%%%%%%%%%%%%
% BumpAmplitude=0.2;
% BumpSpread=1.0;
% Bx=5;  % Bump Center in x direction
% By=5;  % Bump center in Y direction
% Bumpfn=@(x,y) BumpAmplitude*exp(-((x-Bx)^2+(y-By)^2)/(2*BumpSpread^2));
% Bump=b.Energycoefficients(Bumpfn); 
% Population.Ecoeff=Population.Ecoeff+Bump.Ecoeff;

p=p.GetVal(Population);      % Visualise the introduced excitation and resultant population distribution
figure();
p.Plotpop();                

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Propogate the excitation %%%%%%%%%%%%%%

Scattering.atol=10^-3; Scattering.rtol=10^-3;
init_dt=0.001; final_time=10.0;init_time=0.0;dense=0.0001;
Scattering=Scattering.RK853(init_dt,final_time,Population,init_time,dense);
Population.Ecoeff=Scattering.PopChange;
p=p.GetVal(Population);
figure();
p.Plotpop();
save('Thermalization.mat','Scattering.Movie_dense','-v7.3','-nocompression');% This saves the thermalisation data

%% Plot the thermalisation data and make movie

%%%%%%%%%%%%%Data can also be plotted by reading Population*.mat files stored automatically during time propagation. Here we show data plotting using directing Scattering object.

clearvars Movie;

for i=1:size(Scattering.Movie_dense,3) % change the limit of this for loop as per requirement by inspecting contents of Scattering.Movie_dense. Only part of array as per final time is populated.
 Population.Ecoeff= Scattering.Movie_dense(:,:,i);
 % Population.Ecoeff=Var(:,1:N,i);  %%%% Use this when Population*.mat are read to plot data
p=p.GetVal(Population);
figure();                      
p.Plotpop();
Movie(i)=getframe(gcf);
end
 close all;
% %%%%%%%%%%%%%%%%%%%%%%%%% Make Movie %%%%%%%%%%%%%%%%%%%
 c=date;d=num2str(count);
 filename=sprintf('Movie_%s_%s',c,d);
 writerObj = VideoWriter(filename,'MPEG-4');
  writerObj.FrameRate = 10;
open(writerObj);

for i=1:length(Movie)
    frame = Movie(i) ;    
    writeVideo(writerObj, frame);
end
close(writerObj);
end
count=count+1;

