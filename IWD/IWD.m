% Intelligent Water Drops
%------------------------
clc
% clear
% close all
qw=1;
% Problem Statement

%obj4--
% Npar = 13;
% VarLow=[0 0 0 0 0 0 0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 1 1 1 1 1 100 100 100 1];
%obj6--
% Npar = 7;
% VarLow=[-10 -10 -10 -10 -10 -10 -10];
% VarHigh = [10 10 10 10 10 10 10];
%obj1--
% Npar = 7;
% VarLow=[0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 2 2 2];
%obj2--
% Npar = 5;
% VarLow=[78 33 27 27 27];
% VarHigh = [102 45 45 45 45];
%obj 5--
% Npar = 8;
% VarLow=[100 1000 1000 10 10 10 10 10];
% VarHigh = [10000 10000 10000 1000 1000 1000 1000 1000];
% obj3--
% Npar = 10;
% VarLow=[0 0 0 0 0 0 0 0 0 0];
% VarHigh = [10 10 10 10 10 10 10 10 10 10];
% Npar = 5;
% VarLow=[0 0 0 0 0];
% VarHigh = [1 1 1 10 10];
[Npar VarLow VarHigh]=bounds()
while qw<201


% parameters
N=100;        %number of candidates
MaxIter=100;  %number of iterations
av=1;
bv=0.01;
cv=1;
as=1;
bs=0.02;
cs=1;

% initialize a random value as best value
SBest = rand(1,Npar).* (VarHigh - VarLow) + VarLow;
TB=fitnessFunc(SBest);
GB=TB;
t = cputime;

%random initialization of water drops and memory
soil = zeros(N, Npar);
gsoil = zeros(N, Npar);
fsoil = zeros(N, Npar);

F = zeros(N, 1);
for ii = 1:N
    soil(ii,:) = rand(1,Npar).* (VarHigh - VarLow) + VarLow;
    
    % calculate the fitness
    F(ii) = fitnessFunc(soil(ii,:));
end

%Main loop
for it=1:MaxIter
    for ii=1:N
        
      %calculate gsoil
      for j=1:Npar
       if min(soil(:,j))>0
          gsoil(ii,j)=soil(ii,j)-min(soil(:,j));
       else
          gsoil(ii,j)=soil(ii,j);
       end
      end 
      
      %calculate fsoil
      es=0.001;
      for j=1:Npar
         fsoil(ii,j)=1/(es + gsoil(ii,j));
      end
      
      %calculate probability
      for j=1:Npar
         p(j)=fsoil(ii,j)/sum(fsoil(ii,:)); 
      end
      
      %calculate Velocity
      for j=1:Npar
        Viwd(j)=av/(bv +cv +(soil(ii,j))^2); 
      end
      
      HUD=abs(mean(soil)-soil(ii,:));
      
      time=HUD./Viwd;
      
      %calculate delsoil
      delsoil=as./(bs+(cs*time));
      
      %update the soil postion
      soil(ii,:)=soil(ii,:) + delsoil;
      
      %check constraints
      soil(ii,:)=limiter(soil(ii,:),VarHigh,VarLow);
      
      %find fitness
      F(ii)=fitnessFunc(soil(ii,:));
      
      %check for optimality
      if F(ii)<TB
         TB=F(ii);
         SBest=soil(ii,:);
      end
      
    end
    
    % store the best value in each iteration
    GB=[GB TB];
end

t1=cputime;

fprintf('The time taken is %3.2f seconds \n',t1-t);
fprintf('The best value is :');
XBest = SBest
TB

executiontime=t1-t;
% Convergence Plot
%figure(1)
%plot(0:MaxIter,GB, 'linewidth',1.2);
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
bestfitness=max(finalfitness)

for i=1:200
   if(finalfitness(i)==bestfitness)
       optimumX=finalx(i,:);
       optimumY=finaly(i,:);
       optimumFitness=(1/finalfitness(i))-1;
       reqTime=finalexecutiontime(i);
       
   end
end


result=[optimumFitness optimumX optimumY reqTime];





