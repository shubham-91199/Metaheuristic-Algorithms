%Grey Wolf Optimization
% ---------------------
 clc
% clear
% close all
qw=1;
% Problem Statement

% obj3--
% Npar = 10;
% VarLow=[0 0 0 0 0 0 0 0 0 0];
% VarHigh = [10 10 10 10 10 10 10 10 10 10];
%obj5--
% Npar = 8;
% VarLow=[100 1000 1000 10 10 10 10 10];
% VarHigh = [10000 10000 10000 1000 1000 1000 1000 1000];
%obj4--
% Npar = 13;
% VarLow=[0 0 0 0 0 0 0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 1 1 1 1 1 100 100 100 1];
% %obj6--
% Npar = 7;
% VarLow=[-10 -10 -10 -10 -10 -10 -10];
% VarHigh = [10 10 10 10 10 10 10];
% %obj1--
% Npar = 7;
% VarLow=[0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 2 2 2];
% %obj2--
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
WolfSize = 100;
MaxIter = 100;

% initialize a random value as best value
gamma.Position=rand(1,Npar).* (VarHigh - VarLow) + VarLow;
gamma.Fit = fitnessFunc(gamma.Position);
GB=gamma.Fit;

t=cputime;

% Initialization of memory
X = repmat(struct('Position',zeros(1,Npar),'Fit',zeros(1,Npar)),WolfSize,1);
X(1) = gamma; % first wolf is gamma

for ii = 2:WolfSize
    %initial position of wolves
    X(ii).Position= rand(1,Npar).* (VarHigh - VarLow) + VarLow;
    
    % fitness of the solutions
    X(ii).Fit=fitnessFunc(X(ii).Position);
    
    % evaluate gamma
    if (X(ii).Fit < gamma.Fit)
        gamma.Position = X(ii).Position;
        gamma.Fit = X(ii).Fit;
    end
end
% sort wolves according to fitness
[~, sortind] = sort([X.Fit]);
X = X(sortind);

%initialize gamma beta and delta
gamma = X(1);
Beta = X(2);
delta = X(3);

%initialize A and C vectors
a=2*ones(1,Npar);
A1=2*rand(1,Npar).*a -a;
A2=2*rand(1,Npar).*a -a;
A3=2*rand(1,Npar).*a -a;
C1=2*rand(1,Npar);
C2=2*rand(1,Npar);
C3=2*rand(1,Npar);

% Main Loop
for jj = 1:MaxIter
    for ii = 1:WolfSize
        
        %compute encircling vectors
        Da=abs(C1.*gamma.Position-X(ii).Position);
        Db=abs(C2.*Beta.Position-X(ii).Position);
        Dd=abs(C3.*delta.Position-X(ii).Position);
        
        %update wolf position
        X1=gamma.Position - A1.*Da;
        X2=Beta.Position - A2.*Db;
        X3=delta.Position - A3.*Dd;
        X(ii).Position = (X1 + X2 + X3)/3;
        
        %maintian constraints
        X(ii).Position=limiter(X(ii).Position,VarHigh,VarLow);
        
        %update fitness
        X(ii).Fit=fitnessFunc(X(ii).Position);
        
        % update gamma if better
        if X(ii).Fit < gamma.Fit
            gamma = X(ii);
        end
    end
       
    % sort wolves according to fitness
    [~, sortind] = sort([X.Fit]);
    X = X(sortind);
    
    %update beta and delta
    Beta = X(2);
    delta = X(3);
    
    %update A and C vectors
    a=2*(1-jj/MaxIter);
    A1=2*rand(1,Npar).*a -a;
    A2=2*rand(1,Npar).*a -a;
    A3=2*rand(1,Npar).*a -a;
    C1=2*rand(1,Npar);
    C2=2*rand(1,Npar);
    C3=2*rand(1,Npar);
    
    % store the best value in each iteration
    GB = [GB gamma.Fit];
end

t1=cputime;

fprintf('The time taken is %3.2f seconds \n',t1-t);
fprintf('The best value is :');
gamma.Position
gamma.Fit
XBest=gamma.Position;

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






