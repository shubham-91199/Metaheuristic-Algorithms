%Satin Bower bird algorithm
%---------------------------
 clc
% clear
% close all
qw=1;
%obj5--
% Npar = 8;
% VarLow=[100 1000 1000 10 10 10 10 10];
% VarHigh = [10000 10000 10000 1000 1000 1000 1000 1000];
%obj3--
% Npar = 10;
% VarLow=[0 0 0 0 0 0 0 0 0 0];
% VarHigh = [10 10 10 10 10 10 10 10 10 10];
%obj4--
% Npar = 13;
% VarLow=[0 0 0 0 0 0 0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 1 1 1 1 1 100 100 100 1];
%obj6--
% Npar = 7;
% VarLow=[-10 -10 -10 -10 -10 -10 -10];
% VarHigh = [10 10 10 10 10 10 10];
%obj2--
% Npar = 5;
% VarLow=[78 33 27 27 27];
% VarHigh = [102 45 45 45 45];
%obj1--
% Npar = 7;
% VarLow=[0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 2 2 2];
%obj8--
% Npar = 5;
% VarLow=[0 0 0 0 0];
% VarHigh = [1 1 1 10 10];
[Npar VarLow VarHigh]=bounds()
while qw<201
% Problem Statement
%Parameters
nb=100;  %number of bowerbirds
maxIter=100;  %maximum number of iterations
alpha = 0.1; %greatest step size
pMut = 0.05;  %mutation probability
Z = 0.02;    % spread factor

% initialize a random value as best value
BestSol= rand(1,Npar).* (VarHigh - VarLow) + VarLow; %elite 
BestCost= fitnessFunc(BestSol);
GB = BestCost;
t = cputime;

%intialize of solutions and memory
pop = zeros(nb,Npar);
Fit = zeros(nb, 1);
F = zeros(nb, 1);

for ii = 1:nb
    % initialize bowers with random values
    pop(ii,:) = rand(1,Npar).* (VarHigh - VarLow) + VarLow;

    % find the fitness for of the bower
    Fit(ii) = fitnessFunc(pop(ii,:));
    
    % store the best value
    if Fit(ii)<BestCost
        BestSol = pop(ii,:);
        BestCost= Fit(ii);
    end
end

for jj=1:maxIter
    
    %Calculate F from Fitness
    for ii=1:nb
        if Fit(ii) > 0
            F(ii) = 1/(1+Fit(ii));  
        else
            F(ii) = 1 + abs(Fit(ii));
        end
    end
    Fsum=sum(F); % sum of F
   
    % calculate probability of each bower
    prob = F/sum(F);
    
    oldpop = pop;
    oldfit = Fit;

    %deterministic changes
    for ii=1:nb
        for k = 1:Npar
            %select a target bower
            target = Roulette(prob);
            
            % calculate step size
            lambda = alpha/(1+prob(ii));

            % update bower
            pop(ii,k) = pop(ii,k) + lambda*(((pop(target,k) + BestSol(k))/2) - pop(ii,k));

            %mutation
            if rand<pMut
                pop(ii,k) = pop(ii,k) + Z*randn*(VarHigh(k) - VarLow(k));
            end
      
        end
        
        % limit solution within searchspace
        pop(ii,:)=limiter(pop(ii,:),VarHigh, VarLow);

        % caclulate the new fitness
        Fit(ii) = fitnessFunc(pop(ii,:));
    end
    
    %Elitism
    pop = [pop
          oldpop];
    Fit = [Fit
           oldfit];
    
    [Fit, index]=sort(Fit); % sort pop according to fitness
    Fit = Fit(1:nb);
    pop = pop(index, :);
    pop = pop(1:nb, :);
    
    % store the best value
    BestSol = pop(1, :);
    BestCost = Fit(1);
    
    % store the best value in each iteration
    GB= [GB BestCost];
end
t1=cputime;

fprintf('The time taken is %3.2f seconds \n',t1-t);
fprintf('The best value is :');
XBest = BestSol
BestCost
executiontime=t1-t;

% Convergence Plot
% figure(1)
% plot(0:maxIter,GB, 'linewidth',1.2);
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







