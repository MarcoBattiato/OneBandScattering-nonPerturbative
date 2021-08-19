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



classdef Band
    properties (Access=public)
     % Input parametes
        Nx  % Mesh elements in x direction
        Ny  % Mesh elements in y direction
        Lx  % Width of domain in x direction
        Ly  % Width of domain in y direction
     % Mesh parameters
     Cx      % Array of element centers X coordinates
     Cy      % Array of element centers Y coordinates
     lkx     % Array of element x lower limits
     lky     % Array of element y  lower limits
     hkx     % Array of element x upper limits
     hky     % Array of element y upper limits
     widthx  % element width in x direction
     widthy  % element width in y direction
     % Energy parameters
     Elims
     Ecoeff
     % Matrices for Time Propagation
     Ecoeff_next
     Err3
     Err5
    end
    methods
        function obj=Band(Nx,Ny,Lx,Ly)          % Constructor for the Band object
            obj.Nx=Nx; obj.Ny=Ny; obj.Lx=Lx; obj.Ly=Ly;
            
            obj.widthx=obj.Lx/obj.Nx; obj.widthy=obj.Ly/obj.Ny;
            
            obj.Cx=zeros(obj.Nx,1);obj.Cy=zeros(obj.Ny,1);
            obj.lkx=zeros(obj.Nx,1);obj.lky=zeros(obj.Ny,1);obj.hkx=zeros(obj.Nx,1);obj.hky=zeros(obj.Ny,1);
      
            obj.Cx=obj.widthx*((1:obj.Nx)-0.5);obj.Cy=obj.widthy*((1:obj.Ny)-0.5);
            obj.lkx=obj.widthx*((1:Nx)-1);obj.lky=obj.widthy*((1:Ny)-1);
            obj.hkx=obj.widthx*(1:Nx);obj.hky=obj.widthy*(1:Ny);
            

            obj.Elims=zeros(obj.Nx*obj.Ny,2); obj.Ecoeff=zeros(6,obj.Nx*obj.Ny);
            
        end
        
     function obj=Energycoefficients(obj,disp)  % This method calculates the modal representation of dispersion
            k=0;
            for i=1:obj.Ny
                for j=1:obj.Nx
                    k=k+1;
                    U=Band.Psimat(obj.lkx(1,j),obj.lky(1,i),obj.hkx(1,j),obj.lky(1,i),obj.lkx(1,j),obj.hky(1,i),...
                        obj.hkx(1,j),obj.hky(1,i),obj.lkx(1,j),obj.Cy(1,i),obj.Cx(1,j),obj.hky(1,i),obj.Cx(1,j),obj.Cy(1,i),obj.widthx,obj.widthy);
                    fNatural=[disp(obj.lkx(1,j),obj.lky(1,i));disp(obj.hkx(1,j),obj.lky(1,i));disp(obj.lkx(1,j),obj.hky(1,i));...
                        disp(obj.hkx(1,j),obj.hky(1,i));disp(obj.lkx(1,j),obj.Cy(1,i));disp(obj.Cx(1,j),obj.hky(1,i))];
                    fL=U\fNatural;
                    obj.Ecoeff(:,k)=fL;
                    
                end
            end
                    
     end
     function obj=Energylimits(obj)             % This method calculates minimum and maximum value of energy over each element
             k=0;
             A = [];
             Aeq = [];
             b = [];
             beq = [];
            for i=1:obj.Ny
                for j=1:obj.Nx
                    k=k+1;
                    lb = [obj.lkx(1,j) obj.lky(1,i)];
                    ub = [obj.hkx(1,j) obj.hky(1,i)];
                    y = @(kx,ky) obj.Ecoeff(1,k)*Band.psi0(obj.widthx,obj.widthy)+obj.Ecoeff(2,k)*Band.psi1(kx,obj.Cx(1,j),obj.widthx,obj.widthy)+obj.Ecoeff(3,k)*Band.psi2(ky,obj.Cy(1,i),obj.widthx,obj.widthy)+...
                        obj.Ecoeff(4,k)*Band.psi3(kx,ky,obj.Cx(1,j),obj.Cy(1,i),obj.widthx,obj.widthy)+obj.Ecoeff(5,k)*Band.psi4(kx,obj.Cx(1,j),obj.widthx,obj.widthy)+...
                        obj.Ecoeff(6,k)*Band.psi5(ky,obj.Cy(1,i),obj.widthx,obj.widthy);
                    yx = @(x) y(x(1), x(2));
                    [x,obj.Elims(k,1)]= fmincon(yx, [obj.Cx(1,j), obj.Cy(1,i)], A, b, Aeq, beq, lb, ub);
                    yx = @(x) -y(x(1), x(2));
                    [x1,fval]= fmincon(yx, [obj.Cx(1,j), obj.Cy(1,i)], A, b, Aeq, beq, lb, ub);
                    obj.Elims(k,2)=-fval;
                end
            end
     end
     function A=LagrtoNat(obj,fL)               % This method gives the natural representation (actual values) of a function.
         A=zeros(6,size(fL,2));
        for r=1:size(fL,2)
         for i=1:obj.Ny
             for j=1:obj.Nx
                 U=Band.Psimat(obj.lkx(1,j),obj.lky(1,i),obj.hkx(1,j),obj.lky(1,i),obj.lkx(1,j),obj.hky(1,i),...
                     obj.hkx(1,j),obj.hky(1,i),obj.lkx(1,j),obj.Cy(1,i),obj.Cx(1,j),obj.hky(1,i),obj.Cx(1,j),obj.Cy(1,i),obj.widthx,obj.widthy);
                 A(:,r)=U*fL(:,r);
                 
             end
         end
        end
         
         
     end
    end
       methods(Static)
        function U=Psimat(x,y,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,Cx,Cy,wx,wy)  % Generates coefficient matrix
           U=zeros(6,6);
           U(1,1)=Band.psi0(wx,wy); U(1,2)=Band.psi1(x,Cx,wx,wy);U(1,3)=Band.psi2(y,Cy,wx,wy);U(1,4)=Band.psi3(x,y,Cx,Cy,wx,wy);U(1,5)=Band.psi4(x,Cx,wx,wy);U(1,6)=Band.psi5(y,Cy,wx,wy);
           U(2,1)=Band.psi0(wx,wy); U(2,2)=Band.psi1(x1,Cx,wx,wy);U(2,3)=Band.psi2(y1,Cy,wx,wy);U(2,4)=Band.psi3(x1,y1,Cx,Cy,wx,wy);U(2,5)=Band.psi4(x1,Cx,wx,wy);U(2,6)=Band.psi5(y1,Cy,wx,wy);
           U(3,1)=Band.psi0(wx,wy); U(3,2)=Band.psi1(x2,Cx,wx,wy);U(3,3)=Band.psi2(y2,Cy,wx,wy);U(3,4)=Band.psi3(x2,y2,Cx,Cy,wx,wy);U(3,5)=Band.psi4(x2,Cx,wx,wy);U(3,6)=Band.psi5(y2,Cy,wx,wy);
           U(4,1)=Band.psi0(wx,wy); U(4,2)=Band.psi1(x3,Cx,wx,wy);U(4,3)=Band.psi2(y3,Cy,wx,wy);U(4,4)=Band.psi3(x3,y3,Cx,Cy,wx,wy);U(4,5)=Band.psi4(x3,Cx,wx,wy);U(4,6)=Band.psi5(y3,Cy,wx,wy);
           U(5,1)=Band.psi0(wx,wy); U(5,2)=Band.psi1(x4,Cx,wx,wy);U(5,3)=Band.psi2(y4,Cy,wx,wy);U(5,4)=Band.psi3(x4,y4,Cx,Cy,wx,wy);U(5,5)=Band.psi4(x4,Cx,wx,wy);U(5,6)=Band.psi5(y4,Cy,wx,wy);
           U(6,1)=Band.psi0(wx,wy); U(6,2)=Band.psi1(x5,Cx,wx,wy);U(6,3)=Band.psi2(y5,Cy,wx,wy);U(6,4)=Band.psi3(x5,y5,Cx,Cy,wx,wy);U(6,5)=Band.psi4(x5,Cx,wx,wy);U(6,6)=Band.psi5(y5,Cy,wx,wy);
           
        end
        
        function z=psi0(wx,wy)          % These remaining functions are the basis functions used
            z=1/sqrt(wx*wy);
        end
        function z=psi1(x,Cx,wx,wy)
            z=2*sqrt(3)*(x-Cx)/(sqrt(wx*wy)*wx);
        end
        function z=psi2(y,Cy,wx,wy)
            z=2*sqrt(3)*(y-Cy)/(sqrt(wx*wy)*wy);
        end
        function z=psi3(x,y,Cx,Cy,wx,wy)
            z=12*(x-Cx)*(y-Cy)/(sqrt(wx*wy)*wy*wx);
        end
        function z=psi4(x,Cx,wx,wy)
            z=(sqrt(5)/(2*sqrt(wx*wy)))*((12*((x-Cx)/wx)^2)-1);
        end
        function z=psi5(y,Cy,wx,wy)
            z=(sqrt(5)/(2*sqrt(wx*wy)))*((12*((y-Cy)/wy)^2)-1);
        end
    end
    
end