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



classdef Plot       % This object allows for generating 3 dimensional plots with desried resolution 
    properties
        X             
        Y
        Z
        Xplot
        Yplot
        N
        Wi
        Wj
        Cx
        Cy
    end
    methods
        function obj=Plot(Div,Lx,Ly)    % Constructor for plot object. User provides the divisions and length of domain in X and Y directions 
            obj.N=(Div+1)^2;
            %fprintf('Plot points are %d \n',obj.N);
            obj.X=zeros(obj.N,1); obj.Y=zeros(obj.N,1);
            obj.Wi=Lx/Div;obj.Wj=Ly/Div;
            count=1;
            for j=1:Div+1
                for i=1:Div+1
                    obj.X(count,1)=(i-1)*obj.Wi;
                    obj.Y(count,1)=(j-1)*obj.Wj;
                    count=count+1;
                end
            end
            obj.Cx=obj.Wi*((1:Div)-0.5);obj.Cy=obj.Wj*((1:Div)-0.5);
            obj.Xplot=reshape(obj.X,[sqrt(obj.N),sqrt(obj.N)]);
            obj.Yplot=reshape(obj.Y,[sqrt(obj.N),sqrt(obj.N)]);
            
        end
        
        function obj=GetVal(obj,Pop) % This method gives the interpolated value of function at the plot points
             obj.Z=zeros(obj.N,1);
            for i=1:obj.N
                   
                    [m,n]=Plot.Elemfinder(obj.X(i,1),obj.Y(i,1),Pop.widthx,Pop.widthy);
                    m=m-(m>Pop.Nx);n=n-(n>Pop.Ny);
                    
                    Elem=(n-1)*Pop.Nx+m;
                    
                    obj.Z(i,1)=Plot.Value(obj.X(i,1),obj.Y(i,1),Pop.Ecoeff(:,Elem),Pop.Cx(m),Pop.Cy(n),Pop.widthx,Pop.widthy);
                    
            end
            
            obj.Z=reshape(obj.Z,[sqrt(obj.N),sqrt(obj.N)]);
         
        end
        
        function Plotpop(obj) % This method plots the functionion
            
            surf(obj.Xplot,obj.Yplot,obj.Z);
            zlim([-0.001 inf]);
            view(40,20);
            
            s.EdgeColor = 'none';
        end
    end
    
    methods(Static)
        function [i,j]=Elemfinder(x,y,wx,wy) % These widths are band object widths. This method gives the element in the mesh for the given point
          
            i=floor(x/wx)+(rem(x,wx)~=0)+(x==0);
            j=floor(y/wy)+(rem(y,wy)~=0)+(y==0);
        end
        function A=Value(x,y,fcoeff,Cx,Cy,Wx,Wy) % This method calculates the function value at given point coordinates
            A=fcoeff(1)*Plot.psi0(Wx,Wy)+fcoeff(2)*Plot.psi1(x,Cx,Wx,Wy)+fcoeff(3)*Plot.psi2(y,Cy,Wx,Wy)+...
                fcoeff(4)*Plot.psi3(x,y,Cx,Cy,Wx,Wy)+fcoeff(5)*Plot.psi4(x,Cx,Wx,Wy)+fcoeff(6)*Plot.psi5(y,Cy,Wx,Wy);
        end
        function z=psi0(wx,wy) % These remaining functions are the used basis functions
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
