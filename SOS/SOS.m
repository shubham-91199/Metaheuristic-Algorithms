%Symbiotic organism Search 
clc%-------------------------
qw=1;
% clc
% clear
% close all
% warning off
%obj3--
% Npar = 10;
% VarLow=[0 0 0 0 0 0 0 0 0 0];
% VarHigh = [10 10 10 10 10 10 10 10 10 10];

% obj5--
% Npar = 8;
% VarLow=[100 1000 1000 10 10 10 10 10];
% VarHigh = [10000 10000 10000 1000 1000 1000 1000 1000];
%obj4--
% Npar = 13;
% VarLow=[0 0 0 0 0 0 0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 1 1 1 1 1 100 100 100 1];
% obj6--
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
%obj8--
% Npar = 5;
% VarLow=[0 0 0 0 0];
% VarHigh = [1 1 1 10 10];
[Npar VarLow VarHigh]=bounds()
while qw<201
%Problem Statement

% parameters
n=50;
MaxIter=50;

XBest = rand(1,Npar).* (VarHigh - VarLow) + VarLow;
XBestFit = fitnessFunc(XBest);
GB = XBestFit;

%Initialization of memory
X=repmat(struct('Position',zeros(1,Npar),'Velocity',zeros(1,Npar)),n,1);

t=cputime;

%Initialize Random ecosystem
for ii = 1:n
    X(ii).Position = rand(1,Npar).* (VarHigh - VarLow) + VarLow;
    X(ii).Position=limiter(X(ii).Position,VarHigh,VarLow);
    X(ii).Fit = fitnessFunc(X(ii).Position);

    if X(ii).Fit < XBestFit
        XBest = X(ii).Position;
        XBestFit = X(ii).Fit;
    end
end


%Main Loop
for jj=1:MaxIter
    for ii=1:n
        %MUTUALISM PHASE
        j=randperm(n,1);   %random organism to interact
        Mutual_Vector=(X(ii).Position+X(j).Position)/2; % mutual vector
        
        % new solutions
        Xinew = X(ii).Position+rand(1,Npar).*(XBest-Mutual_Vector);
        Xjnew = X(j).Position+rand(1,Npar).*(XBest-Mutual_Vector);
        
        % limit solutions within search space
        Xinew=limiter(Xinew,VarHigh,VarLow);
        Xjnew=limiter(Xjnew,VarHigh,VarLow);
        
        % find new fitness values
        Xinewfit=fitnessFunc(Xinew);
        Xjnewfit=fitnessFunc(Xjnew);

        % update organisms if better
        if Xinewfit<X(ii).Fit
            X(ii).Position=Xinew;
            X(ii).Fit=Xinewfit;
            if X(ii).Fit < XBestFit
              XBest = X(ii).Position;
              XBestFit = X(ii).Fit;
            end
        end
        if Xjnewfit<X(j).Fit
            X(j).Position=Xjnew;
            X(j).Fit=Xjnewfit;
            if X(j).Fit < XBestFit
              XBest = X(j).Position;
              XBestFit = X(j).Fit;
            end
        end

        %COMMENSALISM PHASE
        j=randperm(n,1);   %random organism to interact
        Xinew = X(ii).Position+(rand).*(XBest-X(ii).Position);
        Xinew=limiter(Xinew,VarHigh,VarLow); % limit within search space
        Xinewfit=fitnessFunc(Xinew);  % update fitness
        
        % update organism if better
        if Xinewfit<X(ii).Fit
            X(ii).Position=Xinew;
            X(ii).Fit=Xinewfit;
            if X(ii).Fit < XBestFit
              XBest = X(ii).Position;
              XBestFit = X(ii).Fit;
            end
        end

        %PARASITISM PHASE
        j=randperm(n,1);   %random organism to become parasite
        Xpar=X(j).Position;
        Xm=randperm(Npar,1); %random vector dimension
        Xpar(Xm) = rand* (VarHigh(Xm) - VarLow(Xm)) + VarLow(Xm);
        
        % limit within search space
        Xpar=limiter(Xpar,VarHigh,VarLow);
        Xparfit=fitnessFunc(Xpar);  % update fitness

        % update organism if better
        if Xparfit<X(ii).Fit
            X(ii).Position=Xpar;
            X(ii).Fit=Xparfit;
            if X(ii).Fit < XBestFit
              XBest = X(ii).Position;
              XBestFit = X(ii).Fit;
            end
        end

    end
    
    % store the best value in each iteration
    GB = [GB XBestFit];
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






