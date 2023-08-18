function [bestZ1, bestZ2 , variables] =DE_min(parameters,model)
% DE for min , we should have z1=total costs , so delta is 0
delta=0;eps=0;
%% DE Parameters
Maxit=parameters.Maxit;
npop=parameters.npop;
F=parameters.F;
cr=parameters.cr;
%% Initialization
individual.x=[];
pop=repmat(individual,npop,1);
for i=1:npop
    pop(i).x=CreateRandomSolution(model);
    [pop(i).Z1]=Cost2(pop(i).x,model,eps,delta);
end
%% DE Main loop
it=0;
while it<Maxit   
    it=it+1;
    % Mutation
    for i =1:npop
        popfM=pop;  %copy of pop for mutation
        popfM(i)=[];
        Ps=randsample(popfM,3);
        pop(i).T=Mutation(Ps(1).x,Ps(2).x,Ps(3).x,F);
    end
    % Crossover
    for i =1:npop
        for j=1:model.nVar
            r=rand;
            if r<cr
                pop(i).T(j)=pop(i).x(j);
            end
        end
        pop(i).Z1T=Cost2(pop(i).T,model,eps,delta);
    end
    % Fitness and Selection 
    for i=1:npop
        if pop(i).Z1T<=pop(i).Z1
            pop(i).x=pop(i).T;
            pop(i).Z1=pop(i).Z1T;
        end
    end    
end
[~ , best]=min([pop.Z1]);
[bestZ1,bestZ2,variables]=Cost2(pop(best).x,model,eps,delta);
end

