%Bio-geography based Krill Herd Migration 
%----------------------------------------
clc
% clear
% close all
qw=1;
% Problem Statement
% Npar = 3;
% VarLow=[-5.12 -5.12 -5.12];
% VarHigh = [5.12 5.12 5.12];
%obj1--
% Npar = 7;
% VarLow=[0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 2 2 2];
%obj4--
% Npar = 13;
% VarLow=[0 0 0 0 0 0 0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 1 1 1 1 1 100 100 100 1];
%obj3--
% Npar = 10;
% VarLow=[0 0 0 0 0 0 0 0 0 0];
% VarHigh = [10 10 10 10 10 10 10 10 10 10];
%obj5--
% Npar = 8;
% VarLow=[100 1000 1000 10 10 10 10 10];
% VarHigh = [10000 10000 10000 1000 1000 1000 1000 1000];
%obj6--
% Npar = 7;
% VarLow=[-10 -10 -10 -10 -10 -10 -10];
% VarHigh = [10 10 10 10 10 10 10];
%obj2--
% Npar = 5;
% VarLow=[78 33 27 27 27];
% VarHigh = [102 45 45 45 45];
%obj8--
% Npar = 5;
% VarLow=[0 0 0 0 0];
% VarHigh = [1 1 1 10 10];
[Npar VarLow VarHigh]=bounds()

while qw<201


% parameters
NP=100;       %number of krills
MaxIter=100;  %number of iterations
Vf=0.2;       %foraging velocity
Dmax=0.005;   %max diffusion
Nmax=0.01;   
pmod=0.5;
wf = 0.3;     %inertia for foraging
wn= 0.3;      %inertia for movement
mu=0.3;

XBest = VarLow;
XBestFit = fitnessFunc(XBest);
GB=XBestFit;
t = cputime;

%Initialization only for memory allocation
Krill=repmat(struct('Position',zeros(1,Npar),'Velocity',zeros(1,Npar)),NP,1);

%Random Population
for ii = 1:NP
    Krill(ii).Position = rand(1,Npar).* (VarHigh - VarLow) + VarLow;
    Krill(ii).N=rand(1,Npar);
    Krill(ii).F=rand(1,Npar);
    Krill(ii).Fit = fitnessFunc(Krill(ii).Position);
end

%main loop
for kk = 1:MaxIter
    
   %Sorting the population according to their fitness
   [~,sortind]=sort([Krill.Fit]);
   Krill = Krill(sortind);

   %store the best krill
   if Krill(1).Fit<XBestFit
      XBest=Krill(1).Position;
      XBestFit=Krill(1).Fit;
   end
  
   for ii=1:NP
 
       %movement
       Krill(ii).N=Nmax+wn*(Krill(ii).N);

       %foraging
       beta=XBest-Krill(ii).Position;
       Krill(ii).F=Vf*beta+ wf*(Krill(ii).F);

       %diffusion
       D=Dmax*rand(1,Npar);

       %new position
       Krill(ii).Position=Krill(ii).Position+Krill(ii).N+D+Krill(ii).F;

       %Krill Migration operator
       if rand<pmod
           for jj=1:Npar
              j=randperm(NP,1);
              if rand<mu
                 Krill(ii).Position(jj)=Krill(j).Position(jj);
              end
           end
       end

       %limit within boundaries
       Krill(ii).Positon=limiter(Krill(ii).Position,VarLow,VarHigh);

       %find the fitness value
       Krill(ii).Fit =fitnessFunc(Krill(ii).Position);
   end
   
   % store the best value in each iteration
   GB=[GB XBestFit];
end
t1=cputime;

fprintf('The time taken is %3.2f seconds \n',t1-t);
fprintf('The best value is :');
XBest
XBestFit
executiontime=t1-t;

% Convergence Plot
% figure(1)
% plot(0:MaxIter,GB, 'linewidth',1.2);
% title('Convergence');
% xlabel('Iterations');
% ylabel('Objective Function (Cost)');
% grid('on')

[fitness fitness1 y] = fitnessFunc(XBest)

finalfitness(qw)=fitness;
finaly(qw,:)=y;
finalx(qw,:)=XBest;
finalexecutiontime(qw)=executiontime;
qw=qw+1
end
bestfitness=min(finalfitness)

for i=1:200
   if(finalfitness(i)==bestfitness)
       optimumX=finalx(i,:);
       optimumY=finaly(i,:);
       optimumFitness=(1/finalfitness(i))-1;
       reqTime=finalexecutiontime(i);
       
   end
end


result=[optimumFitness optimumX optimumY reqTime];

       