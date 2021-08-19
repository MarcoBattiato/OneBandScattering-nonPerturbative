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




classdef ScatTensorList < handle
    properties(Access=public)
        path;
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
        Movie_dense
        BumpMag
        BumpMagcheck
        atol
        rtol
        reject
        old_StpSize
        Stp_size
        err_old
        a=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;5.26001519587677318785587544488e-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1.97250569845378994544595329183e-2, 5.91751709536136983633785987549e-2,0,0,...
            0,0,0,0,0,0,0,0,0,0,0;2.95875854768068491816892993775e-2, 0,  8.87627564304205475450678981324e-2,0,0,0,0,0,0,0,0,0,0,0,0;2.41365134159266685502369798665e-1, 0,...
            -8.84549479328286085344864962717e-1, 9.24834003261792003115737966543e-1,0,0,0,0,0,0,0,0,0,0,0;3.7037037037037037037037037037e-2, 0, 0, 1.70828608729473871279604482173e-1,...
            1.25467687566822425016691814123e-1,0,0,0,0,0,0,0,0,0,0;3.7109375e-2, 0, 0, 1.70252211019544039314978060272e-1,6.02165389804559606850219397283e-2, -1.7578125e-2,0,0,0,0,0,0,0,0,0;...
            3.70920001185047927108779319836e-2, 0, 0, 1.70383925712239993810214054705e-1, 1.07262030446373284651809199168e-1,-1.53194377486244017527936158236e-2, 8.27378916381402288758473766002e-3,...
            0,0,0,0,0,0,0,0;6.24110958716075717114429577812e-1, 0, 0, -3.36089262944694129406857109825e+0,-8.68219346841726006818189891453e-1, 2.75920996994467083049415600797e+1, ...
            2.01540675504778934086186788979e+1, -4.34898841810699588477366255144e+1,0,0,0,0,0,0,0;4.77662536438264365890433908527e-1, 0, 0, -2.48811461997166764192642586468e+0,...
            -5.90290826836842996371446475743e-1, 2.12300514481811942347288949897e+1,1.52792336328824235832596922938e+1, -3.32882109689848629194453265587e+1, -2.03312017085086261358222928593e-2,...
            0,0,0,0,0,0;-9.3714243008598732571704021658e-1, 0, 0,5.18637242884406370830023853209e+0, 1.09143734899672957818500254654e+0, -8.14978701074692612513997267357e+0, -1.85200656599969598641566180701e+1,...
            2.27394870993505042818970056734e+1, 2.49360555267965238987089396762e+0, -3.0467644718982195003823669022e+0,0,0,0,0,0;2.27331014751653820792359768449e+0, 0, 0,-1.05344954667372501984066689879e+1,...
            -2.00087205822486249909675718444e+0, -1.79589318631187989172765950534e+1, 2.79488845294199600508499808837e+1,-2.85899827713502369474065508674e+0, -8.87285693353062954433549289258e+0, ...
            1.23605671757943030647266201528e+1, 6.43392746015763530355970484046e-1,0,0,0,0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0;5.61675022830479523392909219681e-2, 0, 0, 0, 0, 0,...
            2.53500210216624811088794765333e-1, -2.46239037470802489917441475441e-1, -1.24191423263816360469010140626e-1,1.5329179827876569731206322685e-1, 8.20105229563468988491666602057e-3,...
            7.56789766054569976138603589584e-3, -8.298e-3,0,0;3.18346481635021405060768473261e-2, 0, 0, 0, 0, 2.83009096723667755288322961402e-2, 5.35419883074385676223797384372e-2,...
            -5.49237485713909884646569340306e-2, 0, 0,-1.08347328697249322858509316994e-4, 3.82571090835658412954920192323e-4, -3.40465008687404560802977114492e-4, 1.41312443674632500278074618366e-1,0;...
            -4.28896301583791923408573538692e-1, 0, 0, 0, 0, -4.69762141536116384314449447206e+0, 7.68342119606259904184240953878e+0, 4.06898981839711007970213554331e+0,3.56727187455281109270669543021e-1,...
            0, 0, 0, -1.39902416515901462129418009734e-3, 2.9475147891527723389556272149e+0, -9.15095847217987001081870187138e+0];
        
        b=[5.42937341165687622380535766363e-2, 0, 0, 0, 0, 4.45031289275240888144113950566e+0,1.89151789931450038304281599044e+0,-5.8012039600105847814672114227e+0, 3.1116436695781989440891606237e-1,...
            -1.52160949662516078556178806805e-1, 2.01365400804030348374776537501e-1,4.47106157277725905176885569043e-2];

        bhh =[0.244094488188976377952755905512e+00, 0.733846688281611857341361741547e+00, 0.220588235294117647058823529412e-01];

        er =[0.1312004499419488073250102996e-01, 0, 0, 0, 0, -0.1225156446376204440720569753e+01,-0.4957589496572501915214079952e+00, 0.1664377182454986536961530415e+01, -0.3503288487499736816886487290e+00,...
            0.3341791187130174790297318841e+00,0.8192320648511571246570742613e-01, -0.2235530786388629525884427845e-01];
       
        end
    methods
        function obj=ScatTensorList(b,startelem,endelem) % Constructor reads the list of possible element combinations generated and creates elements of scattering tensor
            obj.N=b.Nx*b.Ny;SizeScatTensor=0;
          
            obj.IndicesList=dlmread('AMatrix.txt');
            
             obj.C=size(obj.IndicesList,1);
             
            for read=startelem:endelem
                
                filename=sprintf('./Listfiles/List%d.mat',read);
                list=load(filename);
                list=list.D;
                lpsize=size(list,1);
                SizeScatTensor=SizeScatTensor+lpsize;
                Indices=obj.IndicesList;
                NoIndices=obj.C;
                
                ScatCell={};
                widthx=b.widthx;widthy=b.widthy;Gx=obj.Gx;Gy=obj.Gy;
                Nx=b.Nx;
               
                ScatTensorList.parfor_progress(lpsize);
                fprintf('Starting parallel loops for element %d\n',read);
                Bx=list(:,1);
                By=list(:,2);
                Cx=list(:,3);
                Cy=list(:,4);
                Dx=list(:,5);
                Dy=list(:,6);
                Ex=list(:,7);
                Ey=list(:,8);
                
                Count=1;
                
                parfor Element=1:lpsize
                    
                    ScatTens=struct('ElemList',0,'A',zeros(6,4,NoIndices));
                    bx=Bx(Element);
                    by=By(Element);
                    cx=Cx(Element);
                    cy=Cy(Element);
                    dx=Dx(Element);
                    dy=Dy(Element);
                    ex=Ex(Element);
                    ey=Ey(Element);
                    
                    a=Integrator2D(Gx,Gy,widthx,widthy,1);
                    Val=0;
                    a=a.Reseed(0);
                    a=a.Reassign(b,bx,by,cx,cy,dx,dy,ex,ey);
                    
                    a=a.Mapping(1,1,1,1,1,1);
                    
                    a=a.InitializeIntegration();
                    
                    prefactor= 1.0*(1-0.5*((bx==cx)&&(by==cy)))*(1-0.5*((dx==ex)&&(dy==ey)))*(1-0.5*((bx==dx)&&(by==dy))*((cx==ex)&&(cy==ey)));
                    
                    Val=a.Integrate(prefactor);
                    
                    ScatTens.ElemList=[(bx+(by-1)*Nx),(cx+(cy-1)*Nx),(dx+(dy-1)*Nx),(ex+(ey-1)*Nx)];
                    ScatTens.A(1,1,1)=Val;      %% This stores the first row of all matrices
                    ScatTens.A(1,2,1)=Val;
                    ScatTens.A(1,3,1)=Val;
                    ScatTens.A(1,4,1)=Val;
                    %This for loop Stores all orders for 'abcd' as '0000'
                    
                    for f=2:6
                        a=a.Mapping(f,1,1,1,1,1);
                        ScatTens.A(f,1,1)=a.Integrate(prefactor);
                        
                        a=a.Mapping(f,1,1,1,1,2);
                        ScatTens.A(f,2,1)=a.Integrate(prefactor);
                        
                        a=a.Mapping(f,1,1,1,1,3);
                        ScatTens.A(f,3,1)=a.Integrate(prefactor);
                        
                        a=a.Mapping(f,1,1,1,1,4);
                        ScatTens.A(f,4,1)=a.Integrate(prefactor);
                    end
                    
                    for i=2:NoIndices
                        for j=1:6
                            a=a.Mapping(j,Indices(i,1),Indices(i,2),Indices(i,3),Indices(i,4),1);
                            ScatTens.A(j,1,i)=a.Integrate(prefactor);
                            
                            a=a.Mapping(j,Indices(i,1),Indices(i,2),Indices(i,3),Indices(i,4),2);
                            ScatTens.A(j,2,i)=a.Integrate(prefactor);
                            
                            a=a.Mapping(j,Indices(i,1),Indices(i,2),Indices(i,3),Indices(i,4),3);
                            ScatTens.A(j,3,i)=a.Integrate(prefactor);
                            
                            a=a.Mapping(j,Indices(i,1),Indices(i,2),Indices(i,3),Indices(i,4),4);
                            ScatTens.A(j,4,i)=a.Integrate(prefactor);
                        end
                    end
                    
                    ScatCell{Element}=ScatTens;
                    ScatTensorList.parfor_progress;
                    
                end
               
                
                ScatCell=ScatCell(cellfun(@(x) ~isequal(x,0), ScatCell));
                D=cat(2,ScatCell{:});
                filename=sprintf('./ScatTensfiles/Tens%d.mat',read);
                save(filename,'D','-v7.3','-nocompression');
                clearvars D;
            end
            
            fprintf('The size of scattering tensor is %d\n',SizeScatTensor);
            obj.Nscat=size(obj.ScatTens,2);
        end
       
   

        function obj=Timepropagation(obj,dt,time,Pop,Ex,Ey)%,b)% b is required to check if integration conserves energy. This was initial Time Propagation using RK4 and not in use now.
            tic
            k=1;

             Nsteps=round(time/dt);
            ElementBump=(Ey-1)*Pop.Nx+Ex;   % Element for the center of the localised excitation introduced
            ElementCenter=(Pop.Ny-1)*Pop.Nx+(Pop.Nx);
            ElementPt1=(Pop.Ny-Ey-1)*Pop.Nx+(Pop.Nx-Ex);
            ElementPt2=(Ey-1)*Pop.Nx+(Pop.Nx-Ex);
            ElementPt3=(Pop.Ny-Ey-1)*Pop.Nx+Ex;
            obj.Movie=zeros(6,(Pop.Nx*Pop.Ny),(Nsteps+1));
            one=@(x,y) 1;
            Kx=@(x,y) x;
            Ky=@(x,y) y;
            obj.Onefunct=Pop.Energycoefficients(one);
            obj.kx=Pop.Energycoefficients(Kx);
            obj.ky=Pop.Energycoefficients(Ky);
            %obj.Movie(:,:,1)=Pop.Ecoeff;
            obj.BumpMag=zeros(Nsteps+1,5);
            obj.PopChange=Pop.Ecoeff;
            PlotPop=Pop;
            pscat=Plot(100,Pop.Lx,Pop.Ly);
            pscat=pscat.GetVal(PlotPop);
%             obj.BumpMag(k,1)=Plot.Value(Ex,Ey,PlotPop.Ecoeff(:,ElementBump),PlotPop.Cx(Ex),PlotPop.Cy(Ey),PlotPop.widthx,PlotPop.widthy);
%             obj.BumpMag(k,2)=Plot.Value((Pop.Nx),(Pop.Ny-1),PlotPop.Ecoeff(:,ElementCenter),PlotPop.Cx(Pop.Nx),PlotPop.Cy(Pop.Ny-1),PlotPop.widthx,PlotPop.widthy);
%             obj.BumpMag(k,3)=Plot.Value((Pop.Nx-Ex),(Pop.Ny-Ey-1),PlotPop.Ecoeff(:,ElementPt1),PlotPop.Cx(Pop.Nx-Ex),PlotPop.Cy(Pop.Ny-Ey-1),PlotPop.widthx,PlotPop.widthy);
%             obj.BumpMag(k,4)=Plot.Value((Pop.Nx-Ex),(Ey-1),PlotPop.Ecoeff(:,ElementPt2),PlotPop.Cx((Pop.Nx-Ex)),PlotPop.Cy((Ey-1)),PlotPop.widthx,PlotPop.widthy);
%             obj.BumpMag(k,5)=Plot.Value(Ex,(Pop.Ny-Ey-1),PlotPop.Ecoeff(:,ElementPt3),PlotPop.Cx(Ex),PlotPop.Cy((Pop.Ny-Ey-1)),PlotPop.widthx,PlotPop.widthy);
%             %obj.BumpMagcheck(k)=pscat.Max;
            [Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem]=ScatTensorList.Extractobject(obj);
            output=fopen('TimestepStatus.txt','a');
            for t=1:Nsteps 
                k=k+1;
                obj.PopChange=ScatTensorList.RK4(Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem,dt,obj.PopChange);%,t,b);
                PlotPop.Ecoeff=obj.PopChange;
                Var=obj.PopChange;
                pscat=pscat.GetVal(PlotPop);
%                 figure();                       
%                 pscat.Plotpop();
%             obj.Movie(:,:,(t+1))=obj.PopChange;
%             obj.BumpMag(k,1)=Plot.Value(Ex,Ey,PlotPop.Ecoeff(:,ElementBump),PlotPop.Cx(Ex),PlotPop.Cy(Ey),PlotPop.widthx,PlotPop.widthy);
%             obj.BumpMag(k,2)=Plot.Value((Pop.Nx),(Pop.Ny-1),PlotPop.Ecoeff(:,ElementCenter),PlotPop.Cx(Pop.Nx),PlotPop.Cy(Pop.Ny-1),PlotPop.widthx,PlotPop.widthy);
%             obj.BumpMag(k,3)=Plot.Value((Pop.Nx-Ex),(Pop.Ny-Ey-1),PlotPop.Ecoeff(:,ElementPt1),PlotPop.Cx(Pop.Nx-Ex),PlotPop.Cy(Pop.Ny-Ey-1),PlotPop.widthx,PlotPop.widthy);
%             obj.BumpMag(k,4)=Plot.Value((Pop.Nx-Ex),(Ey-1),PlotPop.Ecoeff(:,ElementPt2),PlotPop.Cx((Pop.Nx-Ex)),PlotPop.Cy((Ey-1)),PlotPop.widthx,PlotPop.widthy);
%             obj.BumpMag(k,5)=Plot.Value(Ex,(Pop.Ny-Ey-1),PlotPop.Ecoeff(:,ElementPt3),PlotPop.Cx(Ex),PlotPop.Cy((Pop.Ny-Ey-1)),PlotPop.widthx,PlotPop.widthy);

                fprintf(output,'The completed timestep is %d\n',t);
                if(rem(t,10)==0)
                    file=sprintf('Population%d.mat',t);
                    save(file,'Var','-v7.3','-nocompression');  
                end
            end
            fclose(output);
            toc
        end
        
        function obj=RK4Adaptive(obj,init_dt,final_time,Pop,init_time) % This method was initial implementation of adaptive time stepping.No longer in use. 
            
            % This section is for Checking Integration conservation
            disp=@(x,y)0.11718*sqrt((x-0.5*Pop.Lx)^2+(y-0.5*Pop.Ly)^2);
            one=@(x,y) 1;
            Kx=@(x,y) x;
            Ky=@(x,y) y;
            obj.Onefunct=Pop.Energycoefficients(one);
            obj.kx=Pop.Energycoefficients(Kx);
            obj.ky=Pop.Energycoefficients(Ky);
            bsample=Pop.Energycoefficients(disp);
            
        steps=0;
        b5 = [35.0/384; 0; 500.0/1113; 125.0/192; -2187.0/6784; 11.0/84; 0];
        b5star=[5179.0/57600;0.0;7571.0/16695;393.0/640;-92097.0/339200;187.0/2100;1.0/40];
        c = [0; 0.2; 0.3; 4./5.; 8./9.; 1.; 1.];
        a = [0,0,0,0,0; 1.0/5,0,0,0,0; 3.0/40, 9.0/40,0,0,0;44.0/45, -56.0/15, 32.0/9,0,0; 19372.0/6561, -25360.0/2187,...
             64448.0/6561, -212.0/729,0;9017.0/3168, -355.0/33, 46732.0/5247, 49.0/176, -5103.0/18656];
         
         obj.Stp_size=init_dt; obj.old_StpSize=init_dt; obj.err_old=1e-4; Stp_min=1e-10;
         
         [Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem]=ScatTensor.Extractobject(obj);
         output=fopen('TimestepStatus.txt','a');

         while(init_time<=final_time)

        k1=obj.Stp_size*ScatTensor.CalDelPop(Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem,Pop.Ecoeff);

        k2=obj.Stp_size*ScatTensor.CalDelPop(Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem,(Pop.Ecoeff+a(2,1)*k1));

        k3=obj.Stp_size*ScatTensor.CalDelPop(Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem,(Pop.Ecoeff+a(3,1)*k1+a(3,2)*k2));

        k4=obj.Stp_size*ScatTensor.CalDelPop(Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem,(Pop.Ecoeff+a(4,1)*k1+a(4,2)*k2...
                                                                                                   +a(4,3)*k3));
        k5=obj.Stp_size*ScatTensor.CalDelPop(Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem,(Pop.Ecoeff+a(5,1)*k1...
                                                                                       +a(5,2)*k2+a(5,3)*k3+a(5,4)*k4));
        k6=obj.Stp_size*ScatTensor.CalDelPop(Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem,(Pop.Ecoeff+a(6,1)*k1...
                                                                                  +a(6,2)*k2+a(6,3)*k3+a(6,4)*k4+a(6,5)*k5));                                                                           
        
        Pop.Ecoeff_next=Pop.Ecoeff+b5(1)*k1+b5(2)*k2+b5(3)*k3+b5(4)*k4+b5(5)*k5+b5(6)*k6;
        Pop.Ecoeff_nextStar=Pop.Ecoeff+ b5star(1)*k1+ b5star(2)*k2+ b5star(3)*k3+ b5star(4)*k4+ b5star(5)*k5+ b5star(6)*k6;
        
        err=ScatTensorList.compute_error(obj.atol,obj.rtol,Pop);
      
        if(adapt_Stp_siz(obj,err))
            steps=steps+1;
            fprintf(output,'The completed timestep is %d\n',steps);
            if(err<0.8)
            obj.Stp_size=1.5*obj.Stp_size;
            end
            init_time=init_time+obj.old_StpSize;
           Pop.Ecoeff=Pop.Ecoeff_next;
           Var=Pop.Ecoeff;
          % Check Integration and Plot
%            ScatTensorList.CheckInt(obj.Onefunct.Ecoeff,obj.kx.Ecoeff,obj.ky.Ecoeff,bsample.Ecoeff,k1);
%            p=Plot(100,Pop.Lx,Pop.Ly);
%            p=p.GetVal(Pop);
%            figure();
%            p.Plotpop();
             if(rem(steps,10)==0)
                 file=sprintf('Population%d.mat',steps);
                 save(file,'Var');
             end
          
        elseif(obj.Stp_size<Stp_min)
                fprint("Step Size is too low!!! Exiting\n");
                break;
        end
        
        
        end
        end
        function obj=RK853(obj,init_dt,final_time,Pop,init_time,dense)
            one=@(x,y) 1;
            obj.Onefunct=Pop.Energycoefficients(one);
            dense_time=init_time+dense;
            % This section is for Checking Integration conservation
            %             disp=@(x,y)0.11718*sqrt((x-0.5*Pop.Lx)^2+(y-0.5*Pop.Ly)^2);
%                         Kx=@(x,y) x;
%                         Ky=@(x,y) y;
%                         obj.kx=Pop.Energycoefficients(Kx);
%                         obj.ky=Pop.Energycoefficients(Ky);
            %             bsample=Pop.Energycoefficients(disp);
            
            steps=1;%n=30;error=zeros(n,2);dt=0.01;

            obj.Movie=zeros(6,(Pop.Nx*Pop.Ny)+1,1000); obj.Movie(:,1:(Pop.Nx*Pop.Ny),1)=Pop.Ecoeff;count=1;obj.Movie(1,(Pop.Nx*Pop.Ny)+1,1)=init_time;
            obj.Movie_dense=zeros(6,(Pop.Nx*Pop.Ny)+1,1000); obj.Movie_dense(:,1:(Pop.Nx*Pop.Ny),1)=Pop.Ecoeff;count_dense=1;
            obj.Movie_dense(1,(Pop.Nx*Pop.Ny)+1,1)=init_time;

            obj.Stp_size=init_dt; obj.old_StpSize=init_dt; obj.err_old=1e-4; Stp_min=1e-10;
            
            %global ElemStruct;
            [Combinations,Onefunctcoeff,Indices,NoScat,Nelem]=ScatTensorList.Extractobject(obj);
              
            output=fopen('TimestepStatus.txt','a');
            
            while(init_time<=final_time)
                tic
                k=zeros(6,Nelem,12);
            %while(steps<=n)
                %obj.Stp_size=init_dt+((steps-1)*dt);
                k(:,:,1)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,Pop.Ecoeff);
                
                k(:,:,2)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,(Pop.Ecoeff+(obj.a(2,1)*obj.Stp_size*k(:,:,1))));
                
                k(:,:,3)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,(Pop.Ecoeff+obj.Stp_size*(obj.a(3,1)*k(:,:,1)+obj.a(3,2)*k(:,:,2))));
                
                k(:,:,4)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,(Pop.Ecoeff+obj.Stp_size*(obj.a(4,1)*k(:,:,1)+obj.a(4,2)*k(:,:,2)...
                    +obj.a(4,3)*k(:,:,3))));
                k(:,:,5)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,(Pop.Ecoeff+obj.Stp_size*(obj.a(5,1)*k(:,:,1)...
                    +obj.a(5,2)*k(:,:,2)+obj.a(5,3)*k(:,:,3)+obj.a(5,4)*k(:,:,4))));
                k(:,:,6)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,(Pop.Ecoeff+obj.Stp_size*(obj.a(6,1)*k(:,:,1)...
                    +obj.a(6,2)*k(:,:,2)+obj.a(6,3)*k(:,:,3)+obj.a(6,4)*k(:,:,4)+obj.a(6,5)*k(:,:,5))));
                k(:,:,7)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,(Pop.Ecoeff+obj.Stp_size*(obj.a(7,1)*k(:,:,1)...
                    +obj.a(7,2)*k(:,:,2)+obj.a(7,3)*k(:,:,3)+obj.a(7,4)*k(:,:,4)+obj.a(7,5)*k(:,:,5)+obj.a(7,6)*k(:,:,6))));
                k(:,:,8)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,(Pop.Ecoeff+obj.Stp_size*(obj.a(8,1)*k(:,:,1)...
                    +obj.a(8,2)*k(:,:,2)+obj.a(8,3)*k(:,:,3)+obj.a(8,4)*k(:,:,4)+obj.a(8,5)*k(:,:,5)+obj.a(8,6)*k(:,:,6)+obj.a(8,7)*k(:,:,7))));
                k(:,:,9)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,(Pop.Ecoeff+obj.Stp_size*(obj.a(9,1)*k(:,:,1)...
                    +obj.a(9,2)*k(:,:,2)+obj.a(9,3)*k(:,:,3)+obj.a(9,4)*k(:,:,4)+obj.a(9,5)*k(:,:,5)+obj.a(9,6)*k(:,:,6)+obj.a(9,7)*k(:,:,7)+obj.a(9,8)*k(:,:,8))));
                k(:,:,10)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,(Pop.Ecoeff+obj.Stp_size*(obj.a(10,1)*k(:,:,1)...
                    +obj.a(10,2)*k(:,:,2)+obj.a(10,3)*k(:,:,3)+obj.a(10,4)*k(:,:,4)+obj.a(10,5)*k(:,:,5)+obj.a(10,6)*k(:,:,6)+obj.a(10,7)*k(:,:,7)+obj.a(10,8)*k(:,:,8)+...
                    obj.a(10,9)*k(:,:,9))));
                k(:,:,11)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,(Pop.Ecoeff+obj.Stp_size*(obj.a(11,1)*k(:,:,1)...
                    +obj.a(11,2)*k(:,:,2)+obj.a(11,3)*k(:,:,3)+obj.a(11,4)*k(:,:,4)+obj.a(11,5)*k(:,:,5)+obj.a(11,6)*k(:,:,6)+obj.a(11,7)*k(:,:,7)+obj.a(11,8)*k(:,:,8)+...
                    obj.a(11,9)*k(:,:,9)+obj.a(11,10)*k(:,:,10))));
                k(:,:,12)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,(Pop.Ecoeff+obj.Stp_size*(obj.a(12,1)*k(:,:,1)...
                    +obj.a(12,2)*k(:,:,2)+obj.a(12,3)*k(:,:,3)+obj.a(12,4)*k(:,:,4)+obj.a(12,5)*k(:,:,5)+obj.a(12,6)*k(:,:,6)+obj.a(12,7)*k(:,:,7)+...
                    obj.a(12,8)*k(:,:,8)+obj.a(12,9)*k(:,:,9)+obj.a(12,10)*k(:,:,10)+obj.a(12,11)*k(:,:,11))));
                
                
                Pop8orderSlope=obj.b(1,1)*k(:,:,1)+obj.b(1,2)*k(:,:,2)+obj.b(1,3)*k(:,:,3)+obj.b(1,4)*k(:,:,4)+obj.b(1,5)*k(:,:,5)+obj.b(1,6)*k(:,:,6)...
                    +obj.b(1,7)*k(:,:,7)+obj.b(1,8)*k(:,:,8)+obj.b(1,9)*k(:,:,9)+obj.b(1,10)*k(:,:,10)+obj.b(1,11)*k(:,:,11)+obj.b(1,12)*k(:,:,12);
                
                Pop.Err3=Pop8orderSlope-obj.bhh(1,1)*k(:,:,1)-obj.bhh(1,2)*k(:,:,9)-obj.bhh(1,3)*k(:,:,12);
                Pop.Err5=obj.er(1,1)*k(:,:,1)+obj.er(1,2)*k(:,:,2)+obj.er(1,3)*k(:,:,3)+obj.er(1,4)*k(:,:,4)+obj.er(1,5)*k(:,:,5)+obj.er(1,6)*k(:,:,6)+obj.er(1,7)*k(:,:,7)...
                    +obj.er(1,8)*k(:,:,8)+obj.er(1,9)*k(:,:,9)+obj.er(1,10)*k(:,:,10)+obj.er(1,11)*k(:,:,11)+obj.er(1,12)*k(:,:,12);
                
                Pop.Ecoeff_next=Pop.Ecoeff+obj.Stp_size*Pop8orderSlope;
                
                err=ScatTensorList.compute_error853(obj.atol,obj.rtol,Pop);
%                 
%                 error(steps,1)=obj.Stp_size; error(steps,2)=err;
%                 steps=steps+1; 
%                 
%                 fprintf("%f \t                   %f\n",obj.Stp_size,err);
                
               if(adapt_Stp_siz(obj,err))
                    fprintf(output,'The completed timestep is %d\n',steps);
                    fprintf("The calculated step size is %f\n",obj.Stp_size);
                    init_time=init_time+obj.old_StpSize;
                    fprintf("The Time is now :%f\n",init_time);
                    count_dense=count_dense+1;count=count+1;
                    if(init_time>=dense_time)
                        rcont=ScatTensorList.PrepareDense(Pop.Ecoeff,Pop.Ecoeff_next,k,obj.old_StpSize,Combinations,Onefunctcoeff,Indices,NoScat,Nelem,obj.a);
                         while((init_time>dense_time) && (dense_time <= final_time))
                               obj.Movie_dense(:,1:(Pop.Nx*Pop.Ny),count_dense)=ScatTensorList.Dense_out853(rcont, init_time, dense_time,obj.old_StpSize);
                               obj.Movie_dense(1,(Pop.Nx*Pop.Ny)+1,count_dense)=dense_time;
                                dense_time= dense_time+dense;
                                count_dense=count_dense+1;
                         end
                    end
                    Pop.Ecoeff=Pop.Ecoeff_next;
                   % p=Plot(100,Pop.Lx,Pop.Ly);
                   % p=p.GetVal(Pop);
                    %figure();
                    % p.Plotpop();
                    %Var=Pop.Ecoeff;
                    obj.Movie(:,1:(Pop.Nx*Pop.Ny),count)=Pop.Ecoeff;obj.Movie(1,(Pop.Nx*Pop.Ny)+1,count)=init_time;
                    obj.Movie_dense(:,1:(Pop.Nx*Pop.Ny),count_dense)=Pop.Ecoeff;obj.Movie_dense(1,(Pop.Nx*Pop.Ny)+1,count_dense)=init_time;
                    Var=obj.Movie;Var_dense=obj.Movie_dense;
                    steps=steps+1; 
                    % Check Integration and Plot
                    %            ScatTensorList.CheckInt(obj.Onefunct.Ecoeff,obj.kx.Ecoeff,obj.ky.Ecoeff,bsample.Ecoeff,k1);
                    %            p=Plot(100,Pop.Lx,Pop.Ly);
                    %            p=p.GetVal(Pop);
                    %            figure();
                    %            p.Plotpop();
                    
                    if(rem(steps,10)==0)
                        file=sprintf('Population%d.mat',steps);
                        save(file,'Var','Var_dense','-v7.3','-nocompression');
                    end
                    
                elseif(obj.Stp_size<Stp_min)
                    fprint("Step Size is too low!!! Exiting\n");
                    break;
                end
                toc
                
            end
        end % This methode uses our implementation of adaptive time stepping with DP853. We use this mostly for time propagation.
    
        function output=adapt_Stp_siz(obj,err) % This method calculates the next timestep based on error for the proposed time step.
        safe_fac = 0.9;scale=0.0; minscale = 0.333333;maxscale = 6.0;beta=0.1;alpha =1.0/8.0-beta*0.2; 
if(err<=1)
    obj.old_StpSize= obj.Stp_size;
  if(err == 0)
    scale = maxscale;
  else
   scale = (safe_fac*(obj.err_old^beta)*(err^-alpha));
   scale = min(max(scale, minscale),maxscale);
  end
fprintf("Accepted Step::Scale is %f \n",scale);
 if(obj.reject) 
obj.Stp_size= obj.Stp_size* min(1.0, scale);
 else 
obj.Stp_size=obj.Stp_size*scale;
 end
obj.err_old = max(err, 1.0e-4);
obj.reject=false;
output=true;

else
scale = max(safe_fac*(err^-alpha), minscale);
fprintf("Rejected Step :: Scale is %f \n",scale);
obj.Stp_size=obj.Stp_size*scale;
obj.reject=true;
output=false;
end
  end
        
            
        function A=Testintegral(obj,b) % This method was used for testing purposes. Not needed in actual code.
            
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
        function A=Integration(Onefunct,Delpop)
            
            C=Onefunct.*Delpop;
            A=sum(C(:));
        end
        function A=RK4(Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem,h,Population)%,t,b) % This methode carries our RK4 step.
            %fprintf('The value of Delpop is\n');
            DelPop=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem,Population);
            %obj.DelPopulation=DelPop;
           % disp(DelPop);
%             obj.Integral(t,1)=ScatTensor.Integration(obj.Onefunct.Ecoeff,DelPop);
%             obj.Integral(t,2)=ScatTensor.Integration(b.Ecoeff,DelPop);
%             obj.Integral(t,3)=ScatTensor.Integration(obj.kx.Ecoeff,DelPop);
%             obj.Integral(t,4)=ScatTensor.Integration(obj.ky.Ecoeff,DelPop);
            k1=h*DelPop;
            
            k2=h*ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem,(Population+0.5*k1));
          
            k3=h*ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem,(Population+0.5*k2));
           
            k4=h*ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,ElemStruct,NoScat,Nelem,(Population+k3));
            A=Population+(1/6)*(k1+2*k2+2*k3+k4);
        end
        function DelPop=CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,Population) % This is the actual contraction of the scattering tensor
         
          DelPop=zeros(6,Nelem);
              parfor l=1:Nelem
                   name=strcat('./ScatTensfiles/Tens',int2str(l),'.mat');
                   ElemStruct=load(name);
                   ElemStruct=ElemStruct.D;
                   
                   DelPoploc=zeros(6,Nelem);
                   Z=size(num2cell(ElemStruct),2);
                   for i=1:Z
                     k=[ElemStruct(i).ElemList(1,1),ElemStruct(i).ElemList(1,2),ElemStruct(i).ElemList(1,3),ElemStruct(i).ElemList(1,4)];   
                 
                    set=1;
             
                     if(set)
                         
                         for j=1:Combinations
                             Stat=((Onefunctcoeff(Indices(j,1),k(1))-Population(Indices(j,1),k(1)))*...
                                 (Onefunctcoeff(Indices(j,2),k(2))-Population(Indices(j,2),k(2)))*...
                                 (Population(Indices(j,3),k(3)))*(Population(Indices(j,4),k(4))))-...
                                 ((Population(Indices(j,1),k(1)))*(Population(Indices(j,2),k(2)))*...
                                 (Onefunctcoeff(Indices(j,3),k(3))-Population(Indices(j,3),k(3)))*...
                                 (Onefunctcoeff(Indices(j,4),k(4))-Population(Indices(j,4),k(4))));
                             
                             
                             DelPoploc(:,k(1))=DelPoploc(:,k(1))+(Stat*ElemStruct(i).A(1:6,1,j));
                             DelPoploc(:,k(2))=DelPoploc(:,k(2))+(Stat*ElemStruct(i).A(1:6,2,j));
                             DelPoploc(:,k(3))=DelPoploc(:,k(3))-(Stat*ElemStruct(i).A(1:6,3,j));
                             DelPoploc(:,k(4))=DelPoploc(:,k(4))-(Stat*ElemStruct(i).A(1:6,4,j));
                        
                         end
                       
                      end
              
                   end
                  DelPop(:,:,l)=DelPoploc;
                   DelPoploc=[];
                   k=[];
                   Stat=[];
                   ElemStruct=[];                                
              end
              DelPop=sum(DelPop,3);
              clearvars ElemStruct;
        end
        function Lambda=ScatteringRates(Combinations,Onefunctcoeff,Indices,Nelem,Population) % We calculate scattering rates here
           Lambda=zeros(6,Nelem);
                
              parfor l=1:Nelem
                   name=strcat('./ScatTensfiles/Tens',int2str(l),'.mat');
                   ElemStruct=load(name);
                   ElemStruct=ElemStruct.D;    
                      DelPoploc=zeros(6,Nelem);
                   Z=size(num2cell(ElemStruct),2);
                   fprintf('Starting loop for element %d\n',l);
                   for i=1:Z
                     k=[ElemStruct(i).ElemList(1,1),ElemStruct(i).ElemList(1,2),ElemStruct(i).ElemList(1,3),ElemStruct(i).ElemList(1,4)];   
                     set=1;
                    

                     if(set)
                             
                         for j=1:Combinations
                             StatB=(Onefunctcoeff(Indices(j,1),k(1))*...
                                 (Onefunctcoeff(Indices(j,2),k(2))-Population(Indices(j,2),k(2)))*...
                                 (Population(Indices(j,3),k(3)))*(Population(Indices(j,4),k(4))))+...
                                 (Onefunctcoeff(Indices(j,1),k(1))*(Population(Indices(j,2),k(2)))*...
                                 (Onefunctcoeff(Indices(j,3),k(3))-Population(Indices(j,3),k(3)))*...
                                 (Onefunctcoeff(Indices(j,4),k(4))-Population(Indices(j,4),k(4))));
                             
                             StatC=(Onefunctcoeff(Indices(j,2),k(2))*...
                                 (Onefunctcoeff(Indices(j,1),k(1))-Population(Indices(j,1),k(1)))*...
                                 (Population(Indices(j,3),k(3)))*(Population(Indices(j,4),k(4))))+...
                                 (Onefunctcoeff(Indices(j,2),k(2))*(Population(Indices(j,1),k(1)))*...
                                 (Onefunctcoeff(Indices(j,3),k(3))-Population(Indices(j,3),k(3)))*...
                                 (Onefunctcoeff(Indices(j,4),k(4))-Population(Indices(j,4),k(4))));
                             
                             StatD=(Onefunctcoeff(Indices(j,3),k(3))*...
                                 (Onefunctcoeff(Indices(j,4),k(4))-Population(Indices(j,4),k(4)))*...
                                 (Population(Indices(j,1),k(1)))*(Population(Indices(j,2),k(2))))+...
                                 (Onefunctcoeff(Indices(j,3),k(3))*(Population(Indices(j,4),k(4)))*...
                                 (Onefunctcoeff(Indices(j,1),k(1))-Population(Indices(j,1),k(1)))*...
                                 (Onefunctcoeff(Indices(j,2),k(2))-Population(Indices(j,2),k(2))));
                             
                             StatE=(Onefunctcoeff(Indices(j,4),k(4))*...
                                 (Onefunctcoeff(Indices(j,3),k(3))-Population(Indices(j,3),k(3)))*...
                                 (Population(Indices(j,1),k(1)))*(Population(Indices(j,2),k(2))))+...
                                 (Onefunctcoeff(Indices(j,4),k(4))*(Population(Indices(j,3),k(3)))*...
                                 (Onefunctcoeff(Indices(j,1),k(1))-Population(Indices(j,1),k(1)))*...
                                 (Onefunctcoeff(Indices(j,2),k(2))-Population(Indices(j,2),k(2))));
                             
                            
                              DelPoploc(:,k(1))=DelPoploc(:,k(1))+(StatB*ElemStruct(i).A(1:6,1,j));
                              DelPoploc(:,k(2))=DelPoploc(:,k(2))+(StatC*ElemStruct(i).A(1:6,2,j));
                             DelPoploc(:,k(3))=DelPoploc(:,k(3))+(StatD*ElemStruct(i).A(1:6,3,j));
                             DelPoploc(:,k(4))=DelPoploc(:,k(4))+(StatE*ElemStruct(i).A(1:6,4,j));

                                
                         end
                        
                      end
                   end
                   Lambda(:,:,l)=DelPoploc;
                   fprintf('Done Combinations\n');
                   DelPoploc=[];
                   k=[];
                   Stat=[];
                   ElemStruct=[];     
              end
              Lambda=sum(Lambda,3);
             
        end
        function [Combinations,Onefunctcoeff,Indices,NoScat,Nelem]=Extractobject(obj) % Time propgation reads the Scattering Tensor object and extracts these details for future use.
            Combinations=obj.C;
            Onefunctcoeff=obj.Onefunct.Ecoeff;
            Indices=obj.IndicesList;
            NoScat=obj.Nscat;
            Nelem=obj.N;
            
        end
        function percent=parfor_progress(N) % This function mointors the progress of parfor. It has been directly taken without modifications from: https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fwww.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2Fsubmissions%2F32101%2Fversions%2F3%2Fcontents%2Fparfor_progress.m&embed=web.
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
        end 
        function err=compute_error(atol,rtol,Pop) % This method calculates the error for the proposed time step for first adaptive time stepping implementation. No longer in use.
            err=0.0;
            Delta=Pop.Ecoeff_next-Pop.Ecoeff_nextStar;
           
            for i=1:size(Pop.Ecoeff_next,2)
                  Val=Pop.LagrtoNat(Pop.Ecoeff_next(:,i));
                  ValDelta=Pop.LagrtoNat(Delta(:,i));
                 scale=atol+(rtol*min(abs(min(Val)),abs(1-max(Val))));
                 Del_OverScale=max(abs(max(ValDelta)),abs(min(ValDelta)))/scale;
                 err=err+Del_OverScale^2;
            end
            
            err=err/size(Pop.Ecoeff_next,2);
        end
        
        function err=compute_error853(atol,rtol,Pop) % This method calculates the error for the proposed time step. This methode is in use.
            err3=0.0;err5=0.0;
           
            for i=1:size(Pop.Ecoeff_next,2)
                  Val=Pop.LagrtoNat(Pop.Ecoeff_next(:,i));
                  Val3=Pop.LagrtoNat(Pop.Err3(:,i));Val5=Pop.LagrtoNat(Pop.Err5(:,i));
                 scale=atol+(rtol*min(abs(min(Val)),abs(1-max(Val))));
                 Del_OverScale3=max(abs(max(Val3)),abs(min(Val3)))/scale;
                 Del_OverScale5=max(abs(max(Val5)),abs(min(Val5)))/scale;
                 err3=err3+Del_OverScale3^2;
                 err5=err5+Del_OverScale5^2;
            end
            err3=sqrt(err3)/size(Pop.Ecoeff_next,2);
            err5=sqrt(err5)/size(Pop.Ecoeff_next,2);
            deno=(err3*err3)+(0.01*err5*err5);
            if(deno==0.0)
                deno=1.0;
            end
            err=(err5*err5)/(sqrt(deno));
        end
        function CheckInt(OnefunctEcoeff,kxEcoeff,kyEcoeff,bandEcoeff,k1) % This method was created for testing purposes only.
            fprintf("Particle number change is %f \n",ScatTensorList.Integration(OnefunctEcoeff,k1));
            fprintf("kx momentum number change is %f \n",ScatTensorList.Integration(kxEcoeff,k1));
            fprintf("ky momentum number change is %f \n",ScatTensorList.Integration(kyEcoeff,k1));
            fprintf("Energy change is %f \n",ScatTensorList.Integration(bandEcoeff,k1));
        end
        function rcont=PrepareDense(Ecoeff,Ecoeff_next,k,old_StpSize,Combinations,Onefunctcoeff,Indices,NoScat,Nelem,a) % In dense output this method prepares the required coefficents.
            d41 = -0.84289382761090128651353491142e+1; d46 = 0.56671495351937776962531783590e+0; d47 = -0.30689499459498916912797304727e+1; d48 = 0.23846676565120698287728149680e+1;
            d49 = 0.21170345824450282767155149946e+1;d410 = -0.87139158377797299206789907490e+0;d411 = 0.22404374302607882758541771650e+1;d412 = 0.63157877876946881815570249290e+0;
            d413 = -0.88990336451333310820698117400e-1;d414 = 0.18148505520854727256656404962e+2;d415 = -0.91946323924783554000451984436e+1;d416 = -0.44360363875948939664310572000e+1;
    
            d51 = 0.10427508642579134603413151009e+2;d56 = 0.24228349177525818288430175319e+3;d57 = 0.16520045171727028198505394887e+3;d58 = -0.37454675472269020279518312152e+3;
            d59 = -0.22113666853125306036270938578e+2;d510 = 0.77334326684722638389603898808e+1;d511 = -0.30674084731089398182061213626e+2;d512 = -0.93321305264302278729567221706e+1;
            d513 = 0.15697238121770843886131091075e+2;d514 = -0.31139403219565177677282850411e+2;d515 = -0.93529243588444783865713862664e+1;d516 = 0.35816841486394083752465898540e+2;
    
            d61 = 0.19985053242002433820987653617e+2;d66 = -0.38703730874935176555105901742e+3;d67 = -0.18917813819516756882830838328e+3;d68 = 0.52780815920542364900561016686e+3;
            d69 = -0.11573902539959630126141871134e+2;d610 = 0.68812326946963000169666922661e+1;d611 = -0.10006050966910838403183860980e+1;d612 = 0.77771377980534432092869265740e+0;
            d613 = -0.27782057523535084065932004339e+1;d614 = -0.60196695231264120758267380846e+2;d615 = 0.84320405506677161018159903784e+2;d616 = 0.11992291136182789328035130030e+2;
    
            d71 = -0.25693933462703749003312586129e+2;d76 = -0.15418974869023643374053993627e+3;d77 = -0.23152937917604549567536039109e+3;d78 = 0.35763911791061412378285349910e+3;
            d79 = 0.93405324183624310003907691704e+2;d710 = -0.37458323136451633156875139351e+2;d711 = 0.10409964950896230045147246184e+3;d712 = 0.29840293426660503123344363579e+2;
            d713 = -0.43533456590011143754432175058e+2;d714 = 0.96324553959188282948394950600e+2;d715 = -0.39177261675615439165231486172e+2;d716 = -0.14972683625798562581422125276e+3;

             NewSlope= ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,Ecoeff_next);
                    
    rcont(:,:,1) = Ecoeff;
    rcont(:,:,2) = Ecoeff_next - Ecoeff;
    rcont(:,:,3) = old_StpSize * k(1,1) - rcont(:,:,2);
    rcont(:,:,4) = rcont(:,:,2) - old_StpSize * NewSlope - rcont(:,:,3);
    rcont(:,:,5) = d41 * k(:,:,1)+ d46 * k(:,:,6) + d47 * k(:,:,7) + d48 * k(:,:,8) + d49 * k(:,:,9) + d410 * k(:,:,10)+ d411 * k(:,:,11)+ d412 * k(:,:,12);
    rcont(:,:,6)= d51 * k(:,:,1) + d56 * k(:,:,6) + d57 * k(:,:,7) + d58 * k(:,:,8) + d59 * k(:,:,9) + d510 * k(:,:,10)+ d511 * k(:,:,11)+ d512 * k(:,:,12);
    rcont(:,:,7) = d61 * k(:,:,1) + d66 * k(:,:,6) + d67 * k(:,:,7) + d68 * k(:,:,8) + d69 * k(:,:,9) + d610 * k(:,:,10)+ d611 * k(:,:,11)+ d612 * k(:,:,12);
    rcont(:,:,8) = d71 * k(:,:,1) + d76 * k(:,:,6) + d77 * k(:,:,7) + d78 * k(:,:,8) + d79 * k(:,:,9) + d710 * k(:,:,10)+ d711 * k(:,:,11)+ d712 * k(:,:,12);
    
    
    temp=Ecoeff+old_StpSize*(a(14,1)*k(:,:,1)+a(14,7)*k(:,:,7)+a(14,8)*k(:,:,8)+a(14,9)*k(:,:,9)+a(14,10)*k(:,:,10)+...
         a(14,11)*k(:,:,11)+a(14,12)*k(:,:,12)+a(14,13)*NewSlope);
    k(:,:,10)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,temp);
    
    temp=[];
      temp=Ecoeff+old_StpSize*(a(15,1)*k(:,:,1)+a(15,6)*k(:,:,6)+a(15,7)*k(:,:,7)+a(15,8)*k(:,:,8)+a(15,11)*k(:,:,11)+...
           a(15,12)*k(:,:,12)+a(15,13)*NewSlope+a(15,14)*k(:,:,10));
    k(:,:,11)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,temp);
    
    temp=[];  
    temp=Ecoeff+old_StpSize*(a(16,1)*k(:,:,1)+a(16,6)*k(:,:,6)+a(16,7)*k(:,:,7)+a(16,8)*k(:,:,8)+a(16,9)*k(:,:,9)+...
           a(16,13)*NewSlope+a(16,14)*k(:,:,10)+a(16,15)*k(:,:,11));
    k(:,:,12)=ScatTensorList.CalDelPop(Combinations,Onefunctcoeff,Indices,NoScat,Nelem,temp);
    
    rcont(:,:,5)=old_StpSize*(rcont(:,:,5)+d413*NewSlope+d414*k(:,:,10)+d415*k(:,:,11)+d416*k(:,:,12));
    rcont(:,:,6)=old_StpSize*(rcont(:,:,6)+d513*NewSlope+d514*k(:,:,10)+d515*k(:,:,11)+d516*k(:,:,12));
    rcont(:,:,7)=old_StpSize*(rcont(:,:,7)+d613*NewSlope+d614*k(:,:,10)+d615*k(:,:,11)+d616*k(:,:,12)); 
    rcont(:,:,8)=old_StpSize*(rcont(:,:,8)+d713*NewSlope+d714*k(:,:,10)+d715*k(:,:,11)+d716*k(:,:,12));
        end
        function DenseOutput=Dense_out853(rcont, init_time, dense_time,old_StpSize) % This method writes the dense output.
            s=(dense_time-(init_time-old_StpSize))/old_StpSize;
            s1=1.0-s;
            DenseOutput=rcont(:,:,1)+s*(rcont(:,:,2)+s1*(rcont(:,:,3)+s*(rcont(:,:,4)+s1*(rcont(:,:,5)+s*(rcont(:,:,6)+s1*(rcont(:,:,7)+s*rcont(:,:,8)))))));
        end
end
end
