%Firefly algorithm
%-----------------
 clc
% clear
% close all
% warning off
qw=1;

%obj1--
% Npar = 7;
% VarLow=[0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 2 2 2];
%obj2--
% Npar = 5;
% VarLow=[78 33 27 27 27];
% VarHigh = [102 45 45 45 45];
%obj3--
% Npar = 10;
% VarLow=[0 0 0 0 0 0 0 0 0 0];
% VarHigh = [10 10 10 10 10 10 10 10 10 10];
%obj4--
% Npar = 13;
% VarLow=[0 0 0 0 0 0 0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 1 1 1 1 1 100 100 100 1];
% obj5--
% Npar = 8;
% VarLow=[100 1000 1000 10 10 10 10 10];
% VarHigh = [10000 10000 10000 1000 1000 1000 1000 1000];
%obj6--
% Npar = 7;
% VarLow=[-10 -10 -10 -10 -10 -10 -10];
% VarHigh = [10 10 10 10 10 10 10];
%obj8--
% Npar = 5;
% VarLow=[0 0 0 0 0];
% VarHigh = [1 1 1 10 10];
[Npar VarLow VarHigh]=bounds()
while qw<201
% parameters
n=100;         %number of fireflies
Gen = 100;
alpha=0.3;      % Randomness 0--1 (highly random)
gamma=1.0;      % Absorption coefficient
delta=0.97;      % Randomness reduction (similar to 
                % an annealing schedule)
betamin=0.2;    %minimum value of beta

% initialize a random value as best value
nbest=rand(1,Npar) .* (VarHigh - VarLow) + VarLow;
Lightbest = fitnessFunc(nbest);
GB = Lightbest;

t=cputime;

% Initialization of fireflies and memory
ns=zeros(n,Npar);
Lightn=zeros(n,1);
for ii = 1:n
    
    ns(ii,:)= rand(1,Npar) .* (VarHigh - VarLow) + VarLow;
    Lightn(ii)= fitnessFunc(ns(ii,:));  
    
    % save the best value
    if (Lightn(ii)<Lightbest)
       nbest=ns(ii,:);
       Lightbest = Lightn(ii);
    end

end


for k=1:Gen
    % Reducing Alpha
    alpha=alpha_new(alpha,Gen);

    % Move all fireflies to the better locations
    ns=move(n,ns,Lightn,alpha,betamin,gamma,VarLow,VarHigh);

    % limit the solutions within the search space
    for ii=1:n
       ns(ii,:)=limiter(ns(ii,:),VarHigh,VarLow); 
    end
    
    % Evaluate new solutions (for all n fireflies)
    for ii=1:n
       Lightn(ii) = fitnessFunc(ns(ii,:));
       if (Lightn(ii)<Lightbest)
           nbest=ns(ii,:);
           Lightbest=Lightn(ii);
       end
    end
     
    % store the best value in each iteration
    GB = [GB Lightbest];
end   

t1=cputime;

fprintf('The time taken is %3.2f seconds \n',t1-t);
fprintf('The best value is :');
XBest = nbest
Lightbest

executiontime=t1-t;
% Convergence Plot
% figure(1)
% plot(0:Gen,GB, 'linewidth',1.2);
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










