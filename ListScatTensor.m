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



classdef ListScatTensor < handle
    properties(Access=public)
        Gx=0;
        Gy=0;
        C
        ScatTens
        Nscat
        IndicesList
        Integral
        N
        Onefunct
        PopChange
        kx
        ky
        DelPopulation
        Movie
        BumpMag
        BumpMagcheck
        PrintList
        end
    methods
        function obj=ListScatTensor(b) % Constructor generates the list of possible combination of elements using band dispersion
            obj.N=b.Nx*b.Ny;
            NoElements=obj.N;
           
            widthx=b.widthx;widthy=b.widthy;Gx=obj.Gx;Gy=obj.Gy;
            Nx=b.Nx;Ny=b.Ny;
            Energylimits=b.Elims;
            
           ListScatTensor.parfor_progress(NoElements);
           ListCell=cell(NoElements,1);
            fprintf('Starting loops!\n');
            
            parfor Elementb=1:(NoElements)
                k=0;A=zeros(NoElements^3,8);
              
                for Elementc=1:Elementb
                
                    for Elementd=1:NoElements
                        for Elemente=1:min((Elementd),(NoElements*(Elementb-Elementd)+Elementc))
                            
                           if(((Energylimits(Elementb,2)+Energylimits(Elementc,2))>=(Energylimits(Elementd,1)+Energylimits(Elemente,1)))...
                                    && ((Energylimits(Elementb,1)+Energylimits(Elementc,1))<=(Energylimits(Elementd,2)+Energylimits(Elemente,2))))
                                
                                
                                bx=rem(Elementb,Nx)+Nx*(rem(Elementb,Nx)==0);
                                by=floor(Elementb/Nx)+(1-(rem(Elementb,Nx)==0));
                                
                                cx=rem(Elementc,Nx)+Nx*(rem(Elementc,Nx)==0);
                                cy=floor(Elementc/Nx)+(1-(rem(Elementc,Nx)==0));
                                
                                dx=rem(Elementd,Nx)+Nx*(rem(Elementd,Nx)==0);
                                dy=floor(Elementd/Nx)+(1-(rem(Elementd,Nx)==0));
                                
                                ex=rem(Elemente,Nx)+Nx*(rem(Elemente,Nx)==0);
                                ey=floor(Elemente/Nx)+(1-(rem(Elemente,Nx)==0));
                                

                               if((abs(bx+cx-dx-ex)<=1)&& (abs(by+cy-dy-ey)<=1))
                                
                                Val=0;
                                a=Integrator2D(Gx,Gy,widthx,widthy,1);
                                
                                
                                a=a.Reseed(0);
                                a=a.Reassign(b,bx,by,cx,cy,dx,dy,ex,ey);
                                
                                a=a.Mapping(1,1,1,1,1,1);
                                
                                a=a.InitializeIntegration();
                                
                                Val=a.Integrate(1);
                                
                                if(Val>0)
                                    k=k+1;
                                    
                                    A(k,1)=bx;A(k,2)=by;A(k,3)=cx;
                                    A(k,4)=cy;A(k,5)=dx;A(k,6)=dy;
                                    A(k,7)=ex;A(k,8)=ey;
                                    
                                end
                                end
                           end
                        end
                    end
                end
                A(all(~A,2), : ) = [];
                ListCell{Elementb}=A;
                A=[];
                ListScatTensor.parfor_progress;
            end
             ListCell=ListCell(cellfun(@(x) ~isequal(x,0), ListCell));
    
            fprintf('Writing the list \n');
            elements=0;
          for i=1:size(ListCell,1)
                D=ListCell{i};
                elements=elements+size(D,1);
                filename=sprintf('./Listfiles/List%d.mat',i);
                save(filename,'D','-v7.3','-nocompression');
                clearvars D;
          end
            fprintf('The number of scatterings are %d\n',elements);
          
        end
         
        function A=Testintegral(obj,b) % This method was constructed to carry out some tests. It is not required now.
            
            a=Integrator2D(obj.Gx,obj.Gy,b.widthx,b.widthy);
            
            bx=1;by=1;cx=2;cy=1;dx=1;dy=1;ex=1;ey=1;
            
            a=a.Reassign(b,bx,by,cx,cy,dx,dy,ex,ey);
            
            a=a.Mapping(1,1,1,1,1,5);
            
            a=a.Testmapping(1);
            a=a.InitializeIntegration();
            Val1=a.Integrate()
            
            a=a.Testmapping(2);
            Val2=a.Integrate()
            
            a=a.Testmapping(3);
            Val3=a.Integrate()
            
            a=a.Testmapping(4);
            Val4=a.Integrate()
            
            A=Val1+Val2+Val3+Val4;
          
        end
    end
    methods (Static)
        function A=Integration(Onefunct,Delpop) % This method was constructed to carry out some tests. It is not required now.
            
            C=Onefunct.*Delpop;
            A=sum(C(:));
        end
        function percent=parfor_progress(N)
            error(nargchk(0, 1,nargin, 'struct'));
            if nargin < 1
                N = -1;
            end
            percent = 0;
            w = 50; % Width of progress bar
            if N > 0
                f = fopen('parfor_progress.txt', 'w');
                if f<0
                    error('Do you have write permissions for %s?', pwd);
                end
                fprintf(f, '%d\n', N); % Save N at the top of progress.txt
                fclose(f);
                
                if nargout == 0
                    disp(['  0%[>', repmat(' ', 1, w), ']']);
                end
            elseif N == 0
                delete('parfor_progress.txt');
                percent = 100;
                
                if nargout == 0
                    disp([repmat(char(8), 1, (w+9)), char(10), '100%[', repmat('=', 1, w+1), ']']);
                end
            else
                if ~exist('parfor_progress.txt', 'file')
                    error('parfor_progress.txt not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.txt.');
                end
                
                f = fopen('parfor_progress.txt', 'a');
                fprintf(f, '1\n');
                fclose(f);
                
                f = fopen('parfor_progress.txt', 'r');
                progress = fscanf(f, '%d');
                fclose(f);
                percent = (length(progress)-1)/progress(1)*100;
                
                if nargout == 0
                    perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
                    disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
                end
            end
        end % This function mointors the progress of parfor. It has been directly taken without modifications from: https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fwww.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2Fsubmissions%2F32101%2Fversions%2F3%2Fcontents%2Fparfor_progress.m&embed=web.
    end
end
