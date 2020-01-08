%-------------------------------------------------------------------------
% Sample Matlab Code for the Modified  Wind Driven Optimization.
% Optimization of the Sphere Function in the range of [-5, 5].
% by Dr. Boulesnane Abdennour - abdennour.boulesnane@gmail.com
%-------------------------------------------------------------------------
%
% Please refer to the following journal article in your research papers:
% Boulesnane A, Meshoul S (2014) A Modified Wind Driven Optimization Model for Global Continuous Optimization.
%  Hybrid Artificial Intelligent Systems. HAIS 2015. Lecture Notes in Computer Science, vol 9121, pp 294-304
%-------------------------------------------------------------------------

tic;
clc
clear;
close all;


global param
%-------------------------------------------------------------

% User defined MWDO parameters:
param.popsize = 100;		% population size.
param.npar = 2;			% Dimension of the problem.
param.maxit = 1000;		% Maximum number of iterations.
param.RT = 5;			% RT coefficient.
param.g = 0;			% gravitational constant.
param.alp = 0.1;		% constants in the update eq.
param.c = 0.2;			% coriolis effect.
maxV = 4;			% maximum allowed speed.
dimMin =-5;			% Lower dimension boundary.
dimMax= 5;			% Upper dimension boundary.
Worstfitness=-Inf;
valuebestfitness=Inf;
interV=maxV;
count=0;
sommeold=0;
sommme=[];
%---------------------------------------------------------------
% Initialize MWDO population, position and velocity:
% Randomize population in the range of [-dimMin, dimMax]:
for K=1:param.popsize
    for z=1:param.npar
        pos(K,z)=randvalue(dimMin,dimMax);
    end
end
% Randomize velocity:
for K=1:param.popsize
    for z=1:param.npar
        vel(K,z)=randvalue(-maxV,maxV); 
    end
end
%---------------------------------------------------------------
% Evaluate initial population
for K=1:param.popsize
    fitness(K)=evaluate(pos(K,:));
end
[valuebestfitness,indx]=min(fitness);
gbest=pos(indx,:);


%----------------------------------------------------------------
somme=0;
for j=1:(param.popsize)
    somme=somme+(fitness(j)-max(max(fitness),valuebestfitness));
end
somme=somme+(valuebestfitness-max(max(fitness),valuebestfitness));

for j=1:(param.popsize)
    pression(j)=exp(-param.npar*((fitness(j)-max(max(fitness),valuebestfitness))/somme));
end
Local_Pression=exp(-param.npar*((valuebestfitness-max(max(fitness),valuebestfitness))/somme));


bestold=Inf;
%-----------------------------------------------------------------

% Start iterations :
for it = 2:param.maxit
    
    % Update the velocity:
    for i=1:param.popsize
        % choose random dimensions:
        a = randperm(param.npar);
        % choose velocity based on random dimension:
        velot(i,:) = vel(i,a);
        vel(i,:) = (1-param.alp)*vel(i,:)-(param.g*pos(i,:))+ ...
            abs(1-Local_Pression/pression(i))*((gbest-pos(i,:)).*param.RT)+ ...
            (param.c*velot(i,:)/pression(i));
    end
    
    % Check velocity:
    vel = min(vel, interV);
    vel = max(vel, -interV);
    % Update air parcel positions:
    pos = pos + vel;
    pos = min(pos, dimMax);
    pos = max(pos, dimMin);
    % Evaluate population: (Pressure)
    for K=1:param.popsize
        fitness(K)=evaluate(pos(K,:));
    end
    
%% LOCAL PRESSURE
    somme=0;
    for h=1:(param.popsize)
        somme=somme+(fitness(h)-max(max(fitness),valuebestfitness));
    end
    somme=somme+(valuebestfitness-max(max(fitness),valuebestfitness));
    
    if(sommeold==somme)
        param.g=rand;
        interV=maxV;
    else
        param.g=0;
        for j=1:(param.popsize)
            pression(j)=exp(-param.npar*((fitness(j)-max(max(fitness),valuebestfitness))/somme));
        end
        Local_Pression=exp(-param.npar*((valuebestfitness-max(max(fitness),valuebestfitness))/somme));
    end

    sommeold=somme;
    [maxfitnss,indx]=min(fitness);
    minpres=pression(indx);
    minpos=pos(indx,:); % min location for this iteration

    if maxfitnss<valuebestfitness
        gbest=minpos;
        Local_Pression=minpres;
        valuebestfitness=maxfitnss;
        
    end

%% Velocite intervalle
    if(bestold==valuebestfitness)
        count=count+1;
    else
        count=0;
    end
    
    if(count==5)
        interV= interV/2;
        % param.alp=param.alp+0.5;
        count=0;
    end
    bestold=valuebestfitness;
end
gbest
valuebestfitness
toc


