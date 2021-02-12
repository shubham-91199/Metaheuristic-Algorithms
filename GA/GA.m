% Genetic Algorithm
%------------------
% clc
% clear
% close all
qw=1;

% % Problem Statement
% obj5--
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
%obj2--
% Npar = 5;
% VarLow=[78 33 27 27 27];
% VarHigh = [102 45 45 45 45];
%obj3--
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
N = 100;          % population size
G = 100;          % Number of generations
cx = 0.5;         %Crossover length
Pmut = 0.12;      %probability of mutation
res = 3;          %resolution of the problem (10^-res)

% calculate bias if the lower range is less than 0
off = zeros(1,Npar);
for ii = 1:Npar
    if VarLow(ii) <0
        off(ii) = 0 - VarLow(ii);
    end
end

% calculate the ideal genome length
gg = zeros(1,Npar);
for ii = 1:Npar
   gg(ii) = length(dec2bin((VarHigh(ii)+off(ii))*(10^res))); 
end
genl = max(gg);

% initialize a random value as best value
best = rand(1,Npar).* (VarHigh - VarLow) + VarLow;
bestVal = fitnessFunc(best);
GB = bestVal;
t = cputime;


% Initialization of a random population and memory
pop = zeros(N, Npar);
fitness = zeros(N, 1);
for ii=1:N
    pop(ii,:) = rand(1,Npar).*(VarHigh - VarLow) + VarLow;
    
    % calulate the fitness of the population
    fitness(ii,:) = fitnessFunc(pop(ii,:));
end

%Sorting the population according to their fitness
[fitness,sortind] = sort(fitness);
pop = pop(sortind, :);

%Generations Loop
for k = 1:G

    % crossover
    for i = 1:N
        if rand>(i/N)  % higher fitness get better chance to mate
            j = randperm(N,1);  % choose a random member

            % perfrom crossover
            offspring = crossover(pop(i,:),pop(j,:),genl,cx,res,off);

            % limit the range within constraints
            offspring = limiter(offspring,VarHigh,VarLow);

            % calculate the fitness
            fitnessoff = fitnessFunc(offspring);

            % add the offspring to the population
            pop=cat(1,pop,offspring);
            fitness=cat(1,fitness,fitnessoff);
        end
    end

    % Mutation
    for i = 1:N
        if rand<Pmut         
            % perfrom mutation
            mutant = mutate(pop(i,:),genl,res,off);

            % limit the range within constraints
            mutant = limiter(mutant,VarHigh,VarLow);

            % calculate the fitness
            fitnessoff = fitnessFunc(mutant);

            % add the offspring to the population
            pop=cat(1,pop,mutant);
            fitness=cat(1,fitness,fitnessoff);
        end
    end 

    % Elimination
    %Sorting the population according to their fitness
    [fitness,sortind] = sort(fitness);
    pop = pop(sortind, :);

    % remove weaker populations
    pop = pop(1:N, :);
    fitness = fitness(1:N);

    % Storing best solution
    bestVal = fitness(1);
    best = pop(1,:);

    % store the best value in each iteration
    GB = [GB bestVal];
   
end

t1=cputime;

fprintf('The time taken is %3.2f seconds \n',t1-t);
fprintf('The best value is :');
XBest = best
bestVal

executiontime=t1-t;
% Convergence Plot
% figure(1)
% plot(0:G,GB, 'linewidth',1.2);
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

result=[optimumFitness optimumX optimumY reqTime];







