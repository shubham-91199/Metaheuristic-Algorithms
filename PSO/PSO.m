% Particle Swarm Optimization
%-----------------------------
clc
% clear
% close all
qw=1;
% Problem Statement
% obj4--
% Npar = 13;
% VarLow=[0 0 0 0 0 0 0 0 0 0 0 0 0];
% VarHigh = [1 1 1 1 1 1 1 1 1 100 100 100 1];
%obj5--
% Npar = 8;
% VarLow=[100 1000 1000 10 10 10 10 10];
% VarHigh = [10000 10000 10000 1000 1000 1000 1000 1000];
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

% Parameters
C1 = 2;
C2 = 4-C1;
Inertia = .3;
DampRatio = .95;
ParticleSize = 100;
MaxIter = 10000;

% initialize a random value as best value
GlobalBest = rand(1,Npar).* (VarHigh - VarLow) + VarLow;
GlobalBestCost = fitnessFunc(GlobalBest);
GB = GlobalBestCost;
t = cputime;

% Initialization of particles and memory
Particle=repmat(struct('Position',zeros(1,Npar),'Velocity',zeros(1,Npar)),ParticleSize,1);

for ii = 1:ParticleSize
    
    % initialize with a random position
    Particle(ii).Position = rand(1,Npar).* (VarHigh - VarLow) + VarLow;

    % find the cost for that position
    Particle(ii).Cost = fitnessFunc(Particle(ii).Position);

    % initialize velocity as random
    Particle(ii).Velocity = rand(1,Npar);
    
    % store current position and cost as localbest
    Particle(ii).LocalBest = Particle(ii).Position;
    Particle(ii).LocalBestCost = Particle(ii).Cost;
    
    % update globalbest cost
    if Particle(ii).Cost < GlobalBestCost
        GlobalBest = Particle(ii).Position;
        GlobalBestCost = Particle(ii).Cost;
    end
end

% Main Loop
for jj = 1:MaxIter
    for ii = 1:ParticleSize
        Inertia = Inertia * DampRatio;
        
        % update velocity 
        Particle(ii).Velocity = rand * Inertia * Particle(ii).Velocity + C1 * rand * (Particle(ii).LocalBest - Particle(ii).Position) +  C2 * rand * (GlobalBest - Particle(ii).Position);
        
        % update to new position
        Particle(ii).Position = Particle(ii).Position + Particle(ii).Velocity;

        % limit position within search space
        Particle(ii).Position = limiter(Particle(ii).Position,VarHigh,VarLow);
        
        % calculate the new cost
        Particle(ii).Cost = fitnessFunc(Particle(ii).Position);

        % update local best and global best
        if Particle(ii).Cost < Particle(ii).LocalBestCost
            Particle(ii).LocalBest = Particle(ii).Position;
            Particle(ii).LocalBestCost = Particle(ii).Cost;

            if Particle(ii).Cost < GlobalBestCost
                GlobalBest = Particle(ii).Position;
                GlobalBestCost = Particle(ii).Cost;
            end
        end
    end
    
    % store the best value in each iteration
    GB = [GB GlobalBestCost];
end

t1=cputime;

fprintf('The time taken is %3.2f seconds \n',t1-t);
fprintf('The best value is :');
XBest = GlobalBest
GlobalBestCost

executiontime=t1-t;

% Convergence Plot
% figure(1)
% %plot(0:MaxIter,GB, 'linewidth',1.2);
% title('Convergence');
% xlabel('Iterations');
% ylabel('Objective Function (Cost)');
% grid('on')
%User has to specify

% objectivefunction =(x(1)-10)^2+5*(x(2)-12)^2+x(3)^4+3*(x(4)-11)^2+10*x(5)^6+7*x(6)^2+x(7)^4-4*x(6)*x(7)-10*x(6)-8*x(7);
% y(1) = 127-2*x(1)^2-3*x(2)^4-x(3)-4*x(4)^2-5*x(5);
% y(2) = 282-7*x(1)-3*x(2)-10*x(3)^2-x(4)+x(5);
% y(3) = 196-23*x(1)-x(2)^2-6*x(6)^2+8*x(7);
% y(4) = -4*x(1)^2-x(2)^2+3*x(1)*x(2)-2*x(3)^2-5*x(6)+11*x(7);

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






