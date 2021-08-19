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



classdef Integrator2D
    
    properties (Access=private)
    % coefficients of the basis polynomials
    aA0
    aA1
    aA2
    aA3
    aA4
    aA5
    aAp0
    aAp1
    aAp2
    aAp3
    aAp4
    aAp5
    aB0
    aB1
    aB2
    aB3
    aB4
    aB5
    aBp0
    aBp1
    aBp2
    aBp3
    aBp4
    aBp5
    aC0
    aC1
    aC2
    aC3
    aC4
    aC5
    aCp0
    aCp1
    aCp2
    aCp3
    aCp4
    aCp5
    aD0
    aD1
    aD2
    aD3
    aD4
    aD5
    aDp0
    aDp1
    aDp2
    aDp3
    aDp4
    aDp5
    %Geometric properties
    Gx
    Gy
    lAx
    lAy
    hAx
    hAy
    lBx
    lBy
    hBx
    hBy
    lCx
    lCy
    hCx
    hCy
    lDx
    lDy
    hDx
    hDy
    % Energy polynomial coefficients after mapping
    mu0
    muA1
    muA2
    muA3
    muA4
    muA5
    muB1
    muB2
    muB3
    muB4
    muB5
    muC1
    muC2
    muC3
    muC4
    muC5
    muD1
    muD2
    muD3
    muD4
    muD5
    % Energy coefficients before mapping
    epsA0
    epsA1
    epsA2
    epsA3
    epsA4
    epsA5
    epsB0
    epsB1
    epsB2
    epsB3
    epsB4
    epsB5
    epsC0
    epsC1
    epsC2
    epsC3
    epsC4
    epsC5
    epsD0
    epsD1
    epsD2
    epsD3
    epsD4
    epsD5
    % Storage for Monte Carlo Integration
    precision=0.0001            % minimum alpha
    Npts=2000           % number of Monte Carlo points
    randBy
    randCx
    randCy
    randDx
    randDy            % random number array for all the 5 dimensions to be integrated
                  % Coordinates in the integration
    yBsq
    xCsq
    yCsq
    xDsq
    yDsq             %squares of the cordinates
    alpha
    beta
    gamma
    discr
    solxBplussq
    solxAplussq
    solxBminussq     
    solxAminussq    % Calculations and solutions
    kernelplus
    kernelminus
    w               % scattering amplitude
    
    %%%%% Mapping parameters
    T1
    T2
    T3           % Coordinate transformation matrices.
    U            %U matrix
    
    % Mesh properties
    Wx
    Wy
    Cx
    Cy
    Cx1
    Cy1
    Cx2
    Cy2
    Cx3
    Cy3
    
    end
    properties (Access=public)
    Wplus
    Wminus
    solxBplus
    solxAplus
    solxBminus
    solxAminus
    yA
    yB
    yC
    xD
    yD 
    xC
    end
   methods
       function obj=Integrator2D(Gx,Gy,wx,wy,w) % Constructor Integrator object. This object carries out the core Monte Carlo integration.
       
           
           obj.Gx=Gx; obj.Gy=Gy;obj.Wx=wx;obj.Wy=wy; obj.w=w;
                      
           obj.yA=zeros(obj.Npts,1);
           obj.yB=zeros(obj.Npts,1);obj.xC=zeros(obj.Npts,1);obj.yC=zeros(obj.Npts,1);
           obj.xD=zeros(obj.Npts,1);obj.yD=zeros(obj.Npts,1);
           obj.yBsq=zeros(obj.Npts,1);obj.xCsq=zeros(obj.Npts,1);obj.yCsq=zeros(obj.Npts,1);
           obj.xDsq=zeros(obj.Npts,1);obj.yDsq=zeros(obj.Npts,1);
           
           
           obj.beta=zeros(obj.Npts,1);obj.gamma=zeros(obj.Npts,1);
           obj.discr=zeros(obj.Npts,1);obj.solxBplus=zeros(obj.Npts,1);
           obj.solxBminus=zeros(obj.Npts,1);obj.solxBplussq=zeros(obj.Npts,1);
           obj.solxBminussq=zeros(obj.Npts,1);
           obj.solxAplus=zeros(obj.Npts,1);obj.solxAminus=zeros(obj.Npts,1);
           obj.solxAplussq=zeros(obj.Npts,1);obj.solxAminussq=zeros(obj.Npts,1);
           obj.kernelplus=zeros(obj.Npts,1);obj.kernelminus=zeros(obj.Npts,1);
           obj.Wminus=zeros(obj.Npts,1);obj.Wplus=zeros(obj.Npts,1); 
           
          
           rng(0);
           obj.randBy=rand(obj.Npts,1);obj.randCx=rand(obj.Npts,1);obj.randCy=rand(obj.Npts,1);
           obj.randDx=rand(obj.Npts,1);obj.randDy=rand(obj.Npts,1);
           
           obj.T1=eye(6);obj.T2=eye(6);
           obj.T2(2,2)=-1;obj.T2(3,3)=-1;
           obj.T3=obj.T2; obj.T3(2,1)=-Gx;obj.T3(3,1)=-Gy;
           obj.T3(4,1)=Gx*Gy;obj.T3(4,2)=Gy;obj.T3(4,3)=Gx;
           obj.T3(5,1)=Gx*Gx;obj.T3(5,2)=2*Gx;
           obj.T3(6,1)=Gy*Gy;obj.T3(6,3)=2*Gy;                   %Transformation matrices
           
               
       end
       
       function obj=Reseed(obj,seed) % Changes seed for random number generator
           rng(seed);
           obj.randBy=rand(obj.Npts,1);obj.randCx=rand(obj.Npts,1);obj.randCy=rand(obj.Npts,1);
           obj.randDx=rand(obj.Npts,1);obj.randDy=rand(obj.Npts,1);
       end
       
       function obj=Reassign(obj,b,bx,by,cx,cy,dx,dy,ex,ey) % Assigns necessary coefficients for a combination of elements in calculation of scattering tensor elements. 
           
           obj.Cx=b.Cx(1,bx);obj.Cy=b.Cy(1,by);obj.Cx1=b.Cx(1,cx);obj.Cy1=b.Cy(1,cy);
           obj.Cx2=b.Cx(1,dx);obj.Cy2=b.Cy(1,dy);obj.Cx3=b.Cx(1,ex);obj.Cy3=b.Cy(1,ey);
           
           obj.lAx=b.lkx(1,bx);obj.lAy=b.lky(1,by);obj.hAx=b.hkx(1,bx);obj.hAy=b.hky(1,by);
           obj.lBx=b.lkx(1,cx);obj.lBy=b.lky(1,cy);obj.hBx=b.hkx(1,cx);obj.hBy=b.hky(1,cy);
           obj.lCx=-b.lkx(1,dx);obj.lCy=-b.lky(1,dy);obj.hCx=-b.hkx(1,dx);obj.hCy=-b.hky(1,dy);
           obj.lDx=obj.Gx-b.lkx(1,ex);obj.lDy=obj.Gy-b.lky(1,ey);obj.hDx=obj.Gx-b.hkx(1,ex);obj.hDy=obj.Gy-b.hky(1,ey); 
           
           obj.epsA0=b.Ecoeff(1,(bx+(by-1)*b.Nx));obj.epsA1=b.Ecoeff(2,(bx+(by-1)*b.Nx));obj.epsA2=b.Ecoeff(3,(bx+(by-1)*b.Nx));...
           obj.epsA3=b.Ecoeff(4,(bx+(by-1)*b.Nx)); obj.epsA4=b.Ecoeff(5,(bx+(by-1)*b.Nx));obj.epsA5=b.Ecoeff(6,(bx+(by-1)*b.Nx));
       
           obj.epsB0=b.Ecoeff(1,(cx+(cy-1)*b.Nx));obj.epsB1=b.Ecoeff(2,(cx+(cy-1)*b.Nx));obj.epsB2=b.Ecoeff(3,(cx+(cy-1)*b.Nx));...
           obj.epsB3=b.Ecoeff(4,(cx+(cy-1)*b.Nx));obj.epsB4=b.Ecoeff(5,(cx+(cy-1)*b.Nx));obj.epsB5= b.Ecoeff(6,(cx+(cy-1)*b.Nx));
       
           obj.epsC0=b.Ecoeff(1,(dx+(dy-1)*b.Nx));obj.epsC1=b.Ecoeff(2,(dx+(dy-1)*b.Nx));obj.epsC2=b.Ecoeff(3,(dx+(dy-1)*b.Nx));...
           obj.epsC3=b.Ecoeff(4,(dx+(dy-1)*b.Nx));obj.epsC4=b.Ecoeff(5,(dx+(dy-1)*b.Nx));obj.epsC5=b.Ecoeff(6,(dx+(dy-1)*b.Nx));
           
           obj.epsD0=b.Ecoeff(1,(ex+(ey-1)*b.Nx));obj.epsD1=b.Ecoeff(2,(ex+(ey-1)*b.Nx));obj.epsD2=b.Ecoeff(3,(ex+(ey-1)*b.Nx));...
           obj.epsD3=b.Ecoeff(4,(ex+(ey-1)*b.Nx));obj.epsD4=b.Ecoeff(5,(ex+(ey-1)*b.Nx));obj.epsD5=b.Ecoeff(6,(ex+(ey-1)*b.Nx));
               
       end
       
       
       function obj=Mapping(obj,ap,a,b,c,d,index) % Changes the order of element combinations during construction of scattering tensor elements  
           % Assign mapped Basis function coefficients
           
           coeffA=Integrator2D.Umatrix(obj.Wx,obj.Wy,obj.Cx,obj.Cy)*obj.T1;
           obj.aA0=coeffA(a,1);obj.aA1=coeffA(a,2);obj.aA2=coeffA(a,3);obj.aA3=coeffA(a,4);obj.aA4=coeffA(a,5);obj.aA5=coeffA(a,6);
           
           coeffB=Integrator2D.Umatrix(obj.Wx,obj.Wy,obj.Cx1,obj.Cy1)*obj.T1;
           obj.aB0=coeffB(b,1);obj.aB1=coeffB(b,2);obj.aB2=coeffB(b,3);obj.aB3=coeffB(b,4);obj.aB4=coeffB(b,5);obj.aB5=coeffB(b,6);
           
           coeffC=Integrator2D.Umatrix(obj.Wx,obj.Wy,obj.Cx2,obj.Cy2)*obj.T2;
           obj.aC0=coeffC(c,1);obj.aC1=coeffC(c,2);obj.aC2=coeffC(c,3);obj.aC3=coeffC(c,4);obj.aC4=coeffC(c,5);obj.aC5=coeffC(c,6);
           
           coeffD=Integrator2D.Umatrix(obj.Wx,obj.Wy,obj.Cx3,obj.Cy3)*obj.T3;
           obj.aD0=coeffD(d,1);obj.aD1=coeffD(d,2);obj.aD2=coeffD(d,3);obj.aD3=coeffD(d,4);obj.aD4=coeffD(d,5);obj.aD5=coeffD(d,6);
           
           switch index
               case 1
                 obj.aAp0=coeffA(ap,1);obj.aAp1=coeffA(ap,2);obj.aAp2=coeffA(ap,3);obj.aAp3=coeffA(ap,4);obj.aAp4=coeffA(ap,5);obj.aAp5=coeffA(ap,6);
                 obj.aBp0=1;obj.aBp1=0;obj.aBp2=0;obj.aBp3=0;obj.aBp4=0;obj.aBp5=0;
                 obj.aCp0=1;obj.aCp1=0;obj.aCp2=0;obj.aCp3=0;obj.aCp4=0;obj.aCp5=0;
                 obj.aDp0=1;obj.aDp1=0;obj.aDp2=0;obj.aDp3=0;obj.aDp4=0;obj.aDp5=0;
               case 2
                 obj.aBp0=coeffB(ap,1);obj.aBp1=coeffB(ap,2);obj.aBp2=coeffB(ap,3);obj.aBp3=coeffB(ap,4);obj.aBp4=coeffB(ap,5);obj.aBp5=coeffB(ap,6);
                 obj.aAp0=1;obj.aAp1=0;obj.aAp2=0;obj.aAp3=0;obj.aAp4=0;obj.aAp5=0;
                 obj.aCp0=1;obj.aCp1=0;obj.aCp2=0;obj.aCp3=0;obj.aCp4=0;obj.aCp5=0;
                 obj.aDp0=1;obj.aDp1=0;obj.aDp2=0;obj.aDp3=0;obj.aDp4=0;obj.aDp5=0;
               case 3
                 obj.aCp0=coeffC(ap,1);obj.aCp1=coeffC(ap,2);obj.aCp2=coeffC(ap,3);obj.aCp3=coeffC(ap,4);obj.aCp4=coeffC(ap,5);obj.aCp5=coeffC(ap,6);
                 obj.aBp0=1;obj.aBp1=0;obj.aBp2=0;obj.aBp3=0;obj.aBp4=0;obj.aBp5=0;
                 obj.aAp0=1;obj.aAp1=0;obj.aAp2=0;obj.aAp3=0;obj.aAp4=0;obj.aAp5=0;
                 obj.aDp0=1;obj.aDp1=0;obj.aDp2=0;obj.aDp3=0;obj.aDp4=0;obj.aDp5=0;
               case 4
                 obj.aDp0=coeffD(ap,1);obj.aDp1=coeffD(ap,2);obj.aDp2=coeffD(ap,3);obj.aDp3=coeffD(ap,4);obj.aDp4=coeffD(ap,5);obj.aDp5=coeffD(ap,6);
                 obj.aBp0=1;obj.aBp1=0;obj.aBp2=0;obj.aBp3=0;obj.aBp4=0;obj.aBp5=0;
                 obj.aCp0=1;obj.aCp1=0;obj.aCp2=0;obj.aCp3=0;obj.aCp4=0;obj.aCp5=0;
                 obj.aAp0=1;obj.aAp1=0;obj.aAp2=0;obj.aAp3=0;obj.aAp4=0;obj.aAp5=0;
               case 5
                 obj.aAp0=1;obj.aAp1=0;obj.aAp2=0;obj.aAp3=0;obj.aAp4=0;obj.aAp5=0;
                 obj.aBp0=1;obj.aBp1=0;obj.aBp2=0;obj.aBp3=0;obj.aBp4=0;obj.aBp5=0;
                 obj.aCp0=1;obj.aCp1=0;obj.aCp2=0;obj.aCp3=0;obj.aCp4=0;obj.aCp5=0;
                 obj.aDp0=1;obj.aDp1=0;obj.aDp2=0;obj.aDp3=0;obj.aDp4=0;obj.aDp5=0;
           end
          
           % Assign mapped energy coefficients  %% 
           
           [obj.mu0,obj.muA1,obj.muA2,obj.muA3,obj.muA4,obj.muA5,obj.muB1,obj.muB2,obj.muB3,obj.muB4,obj.muB5,...
               obj.muC1,obj.muC2,obj.muC3,obj.muC4,obj.muC5,obj.muD1,obj.muD2,obj.muD3,obj.muD4,obj.muD5]=...
               Integrator2D.Energymap(obj.epsA0,obj.epsA1,obj.epsA2,obj.epsA3,obj.epsA4,obj.epsA5,obj.epsB0,obj.epsB1,obj.epsB2,...
               obj.epsB3,obj.epsB4,obj.epsB5,obj.epsC0,obj.epsC1,obj.epsC2,obj.epsC3,obj.epsC4,obj.epsC5,...
               obj.epsD0,obj.epsD1,obj.epsD2,obj.epsD3,obj.epsD4,obj.epsD5,coeffA,coeffB,coeffC,coeffD);
           
       end
       
       
       function obj=InitializeIntegration(obj) % Initialising step before Monte Carlo Integration
           % Initialize MC points and their squares
           obj.yB=obj.lBy+(obj.hBy-obj.lBy)*obj.randBy;
           obj.xC=(obj.lCx+(obj.hCx-obj.lCx)*obj.randCx);     
           obj.yC=(obj.lCy+(obj.hCy-obj.lCy)*obj.randCy);
           obj.xD=(obj.lDx+(obj.hDx-obj.lDx)*obj.randDx);
           obj.yD=(obj.lDy+(obj.hDy-obj.lDy)*obj.randDy);
           
           obj.yBsq=obj.yB .^2; obj.xCsq=obj.xC .^2;obj.yCsq=obj.yC .^2;
           obj.xDsq=obj.xD .^2;obj.yDsq=obj.yD .^2;
           
           % Calculate the intermediate terms
           obj.alpha=obj.muA4+obj.muB4;
           obj.yA=-obj.yB-obj.yC-obj.yD;
           obj.beta=obj.muB1-obj.muA1+2*obj.muA4*(obj.xC+obj.xD)-obj.muA3*obj.yA+obj.muB3*obj.yB;
           obj.gamma=obj.mu0+obj.muA2*obj.yA+...
               obj.muA5*obj.yA .^2 -obj.muA3*(obj.xC+obj.xD).*obj.yA+...
               obj.muB2*obj.yB+obj.muB5*obj.yB .^2+2*obj.muA4*obj.xC .*obj.xD+...
               (obj.muC1-obj.muA1)*obj.xC+(obj.muC4+obj.muA4)*obj.xC .^2+obj.muC2*obj.yC+...
               obj.muC5*obj.yC .^2+obj.muC3*obj.xC .*obj.yC+...
               (obj.muD1-obj.muA1)*obj.xD+(obj.muD4+obj.muA4)*obj.xD .^2+obj.muD2*obj.yD+...
               obj.muD5*obj.yD .^2+obj.muD3*obj.xD .*obj.yD;
           obj.discr=obj.beta .^2-4*obj.alpha*obj.gamma;

           
           if(abs(obj.alpha)>obj.precision)
           
               obj.solxBplus=(-obj.beta+sqrt(abs(obj.discr)))/(2*obj.alpha);
               obj.solxBminus=(-obj.beta-sqrt(abs(obj.discr)))/(2*obj.alpha);
               obj.solxAplus=-obj.solxBplus-obj.xC-obj.xD;
               obj.solxAminus=-obj.solxBminus-obj.xC-obj.xD;
               
               obj.solxBplussq=obj.solxBplus .^2;
               obj.solxBminussq=obj.solxBminus .^2;
               obj.solxAplussq=obj.solxAplus .^2;
               obj.solxAminussq=obj.solxAminus .^2;
      
           elseif (abs(obj.beta) >obj.precision)
               
               obj.solxBplus=-obj.gamma./obj.beta;
               obj.solxBminus=(obj.lBx-obj.hBx+obj.lBx).*ones(obj.Npts,1);

               obj.solxAplus=-obj.solxBplus-obj.xC-obj.xD;
               obj.solxAminus=(obj.lAx-obj.hAx+obj.lAx).*ones(obj.Npts,1);
               
               obj.solxBplussq=obj.solxBplus .^2;
               obj.solxBminussq=obj.solxBminus .^2;
               obj.solxAplussq=obj.solxAplus .^2;
               obj.solxAminussq=obj.solxAminus .^2;
               
           else
              
               obj.solxBplus=(obj.lBx-obj.hBx+obj.lBx).*ones(obj.Npts,1);
               obj.solxBminus=(obj.lBx-obj.hBx+obj.lBx).*ones(obj.Npts,1);

               obj.solxAplus=(obj.lAx-obj.hAx+obj.lAx).*ones(obj.Npts,1);
               obj.solxAminus=(obj.lAx-obj.hAx+obj.lAx).*ones(obj.Npts,1);
               
               obj.solxBplussq=obj.solxBplus .^2;
               obj.solxBminussq=obj.solxBminus .^2;
               obj.solxAplussq=obj.solxAplus .^2;
               obj.solxAminussq=obj.solxAminus .^2;
              
           end
           
           obj.Wplus=Integrator2D.CalW(obj.solxAplus,obj.yA,obj.solxBplus,obj.yB,-obj.xC,-obj.yC,-obj.xD,-obj.yD,obj.w,obj.Wx,obj.Wy); 
           obj.Wminus=Integrator2D.CalW(obj.solxAminus,obj.yA,obj.solxBminus,obj.yB,-obj.xC,-obj.yC,-obj.xD,-obj.yD,obj.w,obj.Wx,obj.Wy);
           
           obj.kernelplus=((obj.solxAplus>=obj.lAx)&(obj.solxAplus<=obj.hAx)&...
               (obj.solxBplus>=obj.lBx)&(obj.solxBplus<=obj.hBx)&...
               (obj.yA>=obj.lAy)&(obj.yA<=obj.hAy)&(obj.discr>0))./...
               sqrt(abs(obj.discr));
           
           obj.kernelminus=((obj.solxAminus>=obj.lAx)&(obj.solxAminus<=obj.hAx)&...
               (obj.solxBminus>=obj.lBx)&(obj.solxBminus<=obj.hBx)&...
               (obj.yA>=obj.lAy)&(obj.yA<=obj.hAy)&(obj.discr>0))./...
               sqrt(abs(obj.discr));
           

       end
       
       function value=Integrate(obj,prefactor) % Actual Monte Carlo Integration
           Area=(obj.hCx-obj.lCx)*(obj.hCy-obj.lCy)*(obj.hDx-obj.lDx)*(obj.hDy-obj.lDy)*(obj.hBy-obj.lBy);
           value=prefactor*sum(((obj.Wplus.*(obj.aA0+obj.aA1*obj.solxAplus+obj.aA2*obj.yA+obj.aA3*obj.solxAplus .*obj.yA+...
               obj.aA4*obj.solxAplussq+obj.aA5*obj.yA .^2).*(obj.aAp0+obj.aAp1*obj.solxAplus+obj.aAp2*obj.yA+...
               obj.aAp3*obj.solxAplus .*obj.yA+obj.aAp4*obj.solxAplussq+obj.aAp5*obj.yA .^2).*...
               (obj.aB0+obj.aB1*obj.solxBplus+obj.aB2*obj.yB+obj.aB3*obj.solxBplus .*obj.yB+...
               obj.aB4*obj.solxBplussq+obj.aB5*obj.yB .^2).*...
               (obj.aBp0+obj.aBp1*obj.solxBplus+obj.aBp2*obj.yB+obj.aBp3*obj.solxBplus .*obj.yB+...
               obj.aBp4*obj.solxBplussq+obj.aBp5*obj.yB .^2).*obj.kernelplus)+...
               (obj.Wminus.*(obj.aA0+obj.aA1*obj.solxAminus+obj.aA2*obj.yA+obj.aA3*obj.solxAminus .*obj.yA+...
               obj.aA4*obj.solxAminussq+obj.aA5*obj.yA .^2).*(obj.aAp0+obj.aAp1*obj.solxAminus+obj.aAp2*obj.yA+...
               obj.aAp3*obj.solxAminus .*obj.yA+obj.aAp4*obj.solxAminussq+obj.aAp5*obj.yA .^2).*...
               (obj.aB0+obj.aB1*obj.solxBminus+obj.aB2*obj.yB+obj.aB3*obj.solxBminus .*obj.yB+...
               obj.aB4*obj.solxBminussq+obj.aB5*obj.yB .^2).*...
               (obj.aBp0+obj.aBp1*obj.solxBminus+obj.aBp2*obj.yB+obj.aBp3*obj.solxBminus .*obj.yB+...
               obj.aBp4*obj.solxBminussq+obj.aBp5*obj.yB .^2).*obj.kernelminus)).*...
               (obj.aC0+obj.aC1*obj.xC+obj.aC2*obj.yC+obj.aC3*obj.xC .*obj.yC+obj.aC4*obj.xCsq+obj.aC5*obj.yC .^2).*...
               (obj.aCp0+obj.aCp1*obj.xC+obj.aCp2*obj.yC+obj.aCp3*obj.xC .*obj.yC+obj.aCp4*obj.xCsq+obj.aCp5*obj.yC .^2).*(obj.aD0+...
               obj.aD1*obj.xD+obj.aD2*obj.yD+obj.aD3*obj.xD .*obj.yD+obj.aD4*obj.xDsq+obj.aD5*obj.yD .^2).*...
               (obj.aDp0+obj.aDp1*obj.xD+obj.aDp2*obj.yD+obj.aDp3*obj.xD .*obj.yD+obj.aDp4*obj.xDsq+obj.aDp5*obj.yD .^2))*Area/obj.Npts;
       end
      
       function obj=Testmapping(obj,index) % This method was created for testing. Not needed in the program
           
           switch index
               case 1
                   
                   obj.aA0=obj.mu0;obj.aA1=obj.muA1;obj.aA2=obj.muA2;obj.aA3=obj.muA3;obj.aA4=obj.muA4;obj.aA5=obj.muA5;
                   
                   obj.aB0=1;obj.aB1=0;obj.aB2=0;obj.aB3=0;obj.aB4=0;obj.aB5=0;
                   
                   obj.aC0=1;obj.aC1=0;obj.aC2=0;obj.aC3=0;obj.aC4=0;obj.aC5=0;
                   
                   obj.aD0=1;obj.aD1=0;obj.aD2=0;obj.aD3=0;obj.aD4=0;obj.aD5=0;
                   
               case 2
                   
                   obj.aA0=1;obj.aA1=0;obj.aA2=0;obj.aA3=0;obj.aA4=0;obj.aA5=0;
                   
                   obj.aB0=0;obj.aB1=obj.muB1;obj.aB2=obj.muB2;obj.aB3=obj.muB3;obj.aB4=obj.muB4;obj.aB5=obj.muB5;
                   
                   obj.aC0=1;obj.aC1=0;obj.aC2=0;obj.aC3=0;obj.aC4=0;obj.aC5=0;
                   
                   obj.aD0=1;obj.aD1=0;obj.aD2=0;obj.aD3=0;obj.aD4=0;obj.aD5=0;
                   
               case 3
                   
                   obj.aA0=1;obj.aA1=0;obj.aA2=0;obj.aA3=0;obj.aA4=0;obj.aA5=0;
                   
                   obj.aB0=1;obj.aB1=0;obj.aB2=0;obj.aB3=0;obj.aB4=0;obj.aB5=0;
                   
                   obj.aC0=0;obj.aC1=obj.muC1;obj.aC2=obj.muC2;obj.aC3=obj.muC3;obj.aC4=obj.muC4;obj.aC5=obj.muC5;
                   
                   obj.aB0=1;obj.aB1=0;obj.aB2=0;obj.aB3=0;obj.aB4=0;obj.aB5=0;
                   
               case 4
                   
                   obj.aA0=1;obj.aA1=0;obj.aA2=0;obj.aA3=0;obj.aA4=0;obj.aA5=0;
                   
                   obj.aB0=1;obj.aB1=0;obj.aB2=0;obj.aB3=0;obj.aB4=0;obj.aB5=0;
                   
                   obj.aC0=1;obj.aC1=0;obj.aC2=0;obj.aC3=0;obj.aC4=0;obj.aC5=0;
                   
                   obj.aD0=0;obj.aD1=obj.muD1;obj.aD2=obj.muD2;obj.aD3=obj.muD3;obj.aD4=obj.muD4;obj.aD5=obj.muD5;
                   
           end
       end
   end
       
   methods(Static)
        function A=Umatrix(Wx,Wy,Cx,Cy) % Generates the matrix for coordinate mapping 
           A=eye(6);
           A(1,1)=1/sqrt(Wx*Wy);
           A(2,1)=-2*sqrt(3)*Cx/(sqrt(Wy)*Wx^1.5);A(2,2)=2*sqrt(3)/(sqrt(Wy)*Wx^1.5);
           A(3,1)=-2*sqrt(3)*Cy/(sqrt(Wx)*Wy^1.5);A(3,3)=2*sqrt(3)/(sqrt(Wx)*Wy^1.5);
           A(4,1)=12*Cx*Cy/(Wx^1.5*Wy^1.5);A(4,2)=-12*Cy/(Wx^1.5*Wy^1.5);A(4,3)=-12*Cx/(Wx^1.5*Wy^1.5);
           A(4,4)=12/(Wx^1.5*Wy^1.5);
           A(5,1)=(sqrt(5)/(2*sqrt(Wx)*sqrt(Wy)))*((12*Cx^2/Wx^2)-1);A(5,2)=-12*sqrt(5)*Cx/(sqrt(Wx)*sqrt(Wy)*Wx^2);
           A(5,5)=6*sqrt(5)/(sqrt(Wx)*sqrt(Wy)*Wx^2);
           A(6,1)=(sqrt(5)/(2*sqrt(Wx)*sqrt(Wy)))*((12*Cy^2/Wy^2)-1);A(6,3)=-12*sqrt(5)*Cy/(sqrt(Wx)*sqrt(Wy)*Wy^2);
           A(6,6)=6*sqrt(5)/(sqrt(Wx)*sqrt(Wy)*Wy^2);
           
       end
       function [mu0,muA1,muA2,muA3,muA4,muA5,muB1,muB2,muB3,muB4,muB5,...
               muC1,muC2,muC3,muC4,muC5,muD1,muD2,muD3,muD4,muD5]=...
               Energymap(epsA0,epsA1,epsA2,epsA3,epsA4,epsA5,epsB0,epsB1,epsB2,epsB3,epsB4,epsB5,epsC0,epsC1,...
               epsC2,epsC3,epsC4,epsC5,epsD0,epsD1,epsD2,epsD3,epsD4,epsD5,coeffA,coeffB,coeffC,coeffD) % Generates matrix for mapping the coefficients of dispersion.
           epsilonA=[epsA0,epsA1,epsA2,epsA3,epsA4,epsA5]*coeffA;
           epsilonB=[epsB0,epsB1,epsB2,epsB3,epsB4,epsB5]*coeffB;
           epsilonC=-[epsC0,epsC1,epsC2,epsC3,epsC4,epsC5]*coeffC;
           epsilonD=-[epsD0,epsD1,epsD2,epsD3,epsD4,epsD5]*coeffD;
           muA1=epsilonA(1,2);muA2=epsilonA(1,3);muA3=epsilonA(1,4);muA4=epsilonA(1,5);muA5=epsilonA(1,6);
           muB1=epsilonB(1,2);muB2=epsilonB(1,3);muB3=epsilonB(1,4);muB4=epsilonB(1,5);muB5=epsilonB(1,6);
           muC1=epsilonC(1,2);muC2=epsilonC(1,3);muC3=epsilonC(1,4);muC4=epsilonC(1,5);muC5=epsilonC(1,6);
           muD1=epsilonD(1,2);muD2=epsilonD(1,3);muD3=epsilonD(1,4);muD4=epsilonD(1,5);muD5=epsilonD(1,6);
           mu0=epsilonA(1,1)+epsilonB(1,1)+epsilonC(1,1)+epsilonD(1,1);
           
       end
       function W=CalW(x1,y1,x2,y2,x3,y3,x4,y4,w,wx,wy) % Calculates the value of scattering amplitude w
                one=ones(size(x1,1),1);
               temp1=((sqrt(sqrt((x1-x3).^2+(y1-y3).^2)/sqrt(wx^2+wy^2)).*(sqrt((x1-x3).^2+(y1-y3).^2)<=sqrt(wx^2+wy^2)))+...
                         (one.*(sqrt((x1-x3).^2+(y1-y3).^2)>sqrt(wx^2+wy^2))));
               temp2=(sqrt(sqrt((x2-x4).^2+(y2-y4).^2)/sqrt(wx^2+wy^2)).*(sqrt((x2-x4).^2+(y2-y4).^2)<=sqrt(wx^2+wy^2))+...
                         (one.*(sqrt((x2-x4).^2+(y2-y4).^2)>sqrt(wx^2+wy^2))));
               temp3=(sqrt(sqrt((x1-x4).^2+(y1-y4).^2)/sqrt(wx^2+wy^2)).*(sqrt((x1-x4).^2+(y1-y4).^2)<=sqrt(wx^2+wy^2))+...
                         (one.*(sqrt((x1-x4).^2+(y1-y4).^2)>sqrt(wx^2+wy^2))));
               temp4=(sqrt(sqrt((x2-x3).^2+(y2-y3).^2)/sqrt(wx^2+wy^2)).*(sqrt((x2-x3).^2+(y2-y3).^2)<=sqrt(wx^2+wy^2))+...
                         (one.*(sqrt((x2-x3).^2+(y2-y3).^2)>sqrt(wx^2+wy^2))));
           W=w.*temp1.*temp2.*temp3.*temp4;
       end
       
   end
end
