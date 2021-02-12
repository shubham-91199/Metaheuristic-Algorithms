%Big Bang Big Crunch Algorithm
%-----------------------------
 clc
% clear
% close all
qw=1;
% Problem Statement
[Npar VarLow VarHigh]=bounds()
%obj4--
% Npar = 13;
% VarLow=[0 0 0 0 0 0 0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 1 1 1 1 1 100 100 100 1];
% obj3--
% Npar = 10;
% VarLow=[0 0 0 0 0 0 0 0 0 0];
% VarHigh = [10 10 10 10 10 10 10 10 10 10];
% obj5--
% Npar = 8;
% VarLow=[100 1000 1000 10 10 10 10 10];
% VarHigh = [10000 10000 10000 1000 1000 1000 1000 1000];
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
% %obj7--
% Npar = 8;
% VarLow=[0 1.2 20 9 6.5 0 0 0];
% VarHigh = [5 2.4 60 9.3 7 294000 294000 277200];
% obj8--
% Npar = 5;
% VarLow=[0 0 0 0 0];
% VarHigh = [1 1 1 10 10];
while qw<201

%BBBC parameters
N=100;       %number of candidates
MaxIter=100;  %number of iterations

% initialize a random value as best value
XBest = rand(1,Npar).* (VarHigh - VarLow) + VarLow;
FBest=fitnessFunc(XBest);
GB=FBest;
t = cputime;

%intialize solutions and memory
X = zeros(N, Npar);
F = zeros(N, 1);

for ii = 1:N
    X(ii,:) = rand(1,Npar).* (VarHigh - VarLow) + VarLow;
    
    % calculate the fitness of solutions
    F(ii) = fitnessFunc(X(ii,:));
end

%Main Loop
for it=1:MaxIter

    %Find the centre of mass 
    %-----------------------

    %numerator term
    num=zeros(1,Npar);
    for ii=1:N
        for jj=1:Npar
            num(jj)=num(jj)+(X(ii,jj)/F(ii));
        end
    end

    %denominator term
    den=sum(1./F);

    %centre of mass
    Xc=num/den; 

    %generate new solutions
    %----------------------
    for ii=1:N

        %new solution from centre of mass
        for jj=1:Npar      
            New=X(ii,:);
            New(jj)=Xc(jj)+((VarHigh(jj)*rand)/it^2);
        end

        %boundary constraints
        New=limiter(New,VarHigh,VarLow);
        %new fitness
        newFit=fitnessFunc(New);

        %check whether the solution is better than previous solution
        if newFit<F(ii)
            X(ii,:)=New;
            F(ii)=newFit;
            if F(ii)<FBest
                XBest=X(ii,:);
                FBest=F(ii);   
            end
        end

    end

    % store the best value in each iteration
    GB=[GB FBest];
end

t1=cputime;

fprintf('The time taken is %3.2f seconds \n',t1-t);
fprintf('The best value is :');
XBest;
FBest;

executiontime=t1-t;
% Convergence Plot
% figure(1)
% plot(0:MaxIter,GB, 'linewidth',1.2);
% title('Convergence');
% xlabel('Iterations');
% ylabel('Objective Function (Cost)');
% grid('on')

[fitness fitness1 y] = fitnessFunc(XBest);

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


result=[optimumFitness optimumX reqTime];




