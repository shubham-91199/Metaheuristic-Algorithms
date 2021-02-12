% Magnetotatic Bacteria Optimizer
%--------------------------------
% clc
% clear
% close all
qw=1;
% Problem Statement

% obj4--
% Npar = 13;
% VarLow=[0 0 0 0 0 0 0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 1 1 1 1 1 100 100 100 1];
%obj6-
% Npar = 7;
% VarLow=[-10 -10 -10 -10 -10 -10 -10];
% VarHigh =[10 10 10 10 10 10 10];
%obj1-
% Npar = 7;
% VarLow=[0 0 0 0 0 0 0];
% VarHigh =[1 1 1 1 2 2 2];
%OBJ2-
% Npar = 5;
% VarLow=[78 33 27 27 27];
% VarHigh = [102 45 45 45 45];
% obj5--
% Npar = 8;
% VarLow=[100 1000 1000 10 10 10 10 10];
% VarHigh = [10000 10000 10000 1000 1000 1000 1000 1000];
% obj3--
% Npar = 10;
% VarLow=[0 0 0 0 0 0 0 0 0 0];
% VarHigh = [10 10 10 10 10 10 10 10 10 10];
%obj8--
% Npar = 5;
% VarLow=[0 0 0 0 0];
% VarHigh = [1 1 1 10 10];
[Npar VarLow VarHigh]=bounds()

while qw<201

% parameters
MaxIter=100;         %Maximum iteration
N=50;          %Initial population
c1 = 20;
c2 =0.003;
mp =0.6;
B = 1;

% initialize a random value as best value
XBest = rand(1,Npar).* (VarHigh - VarLow) + VarLow;
XBestCost = fitnessFunc(XBest);
GB = XBestCost;

t=cputime;

% Initialize of magnetotatic bacteria and memory
X = zeros(N,Npar);
F = zeros(N, 1);
Fnorm = F;

for ii=1:N       
    % initialize with a random solution
    X(ii,:)=rand(1,Npar) .* (VarHigh - VarLow) + VarLow;

    % find the fitness of the solution
    F(ii)=fitnessFunc(X(ii,:));
end

% Main Loop
for it=1:MaxIter
   
   %normalized fitnesss
   for ii= 1:N
      Fnorm(ii) = (F(ii) - min(F))/(max(F) - min(F));  
   end
    
   %interaction distance
   for ii = 1:N
       for jj = 1:N
          for kk = 1:Npar
             D(ii,jj,kk) = (X(ii,kk) - X(jj,kk));  
          end
       end
   end
   
   %Moments
   for ii = 1:N
       for jj=1:N
         for kk = 1:Npar
           p = round(rand*N +0.5);
           q = round(rand*N +0.5);
           for dd = 1:Npar
             DD(dd) = D(ii,jj,dd); 
           end
           E(ii,jj,kk) = (D(ii,jj,kk)/(1+c1*norm(DD) + c2*D(p,q)))^3;
         end
       end
   end
   M = E/B;
   
   %moment regulation
   for ii=1:N
     if (rand<mp)
       X(ii,:) = X(ii,:) + (XBest - X(ii,:)).*rand(1,Npar);
     else
        for jj = 1:Npar
           r = round(rand*N +0.5);
           X(ii,jj) =  X(r,jj);
        end
     end
     
     % limit solution within search space
     X(ii,:) = limiter(X(ii,:),VarHigh,VarLow);
     
     % find new fitness
     F(ii) = fitnessFunc(X(ii,:));
   end
   
   %MTS Replacement
   [F, sortOrder] = sort(F);
   X = X(sortOrder,:);
   
   for ii= round(N/2):N
       if rand<0.5
         for jj = 1:Npar
           r = round(rand*N +0.5);
           X(ii,jj) = X(ii,jj) + rand*M(ii,r,jj);
         end
       else
          X(ii,:)=rand(1,Npar) .* (VarHigh - VarLow) + VarLow;
       end
       
       % limit solution within search space
       X(ii,:) = limiter(X(ii,:),VarHigh,VarLow);
       
       % find new fitness
       F(ii) = fitnessFunc(X(ii,:));
   end
   
   [F, sortOrder] = sort(F);
   X = X(sortOrder,:);
   
   if (F(1) < XBestCost)
      XBestCost = F(1);
      XBest = X(1,:);
   end
   
    % store the best value in each iteration
    GB=[GB XBestCost];
end

t1=cputime;

fprintf('The time taken is %3.2f seconds \n',t1-t);
fprintf('The best value is :');
XBest
XBestCost

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






