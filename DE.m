clc;clear;close;
%% Insert Data
model=SSModel_7s();
TP='7s';
%% Parameters 
npar=40;
Maxit=200;
npop=55;
F=0.95;
cr=0.3;
%% one objective runs
tic
parameters.Maxit=Maxit;parameters.npop=npop;
parameters.F=F;parameters.cr=cr;

[minZ1,nadirZ2]=DE_min(parameters,model);
[nadirZ1,maxZ2]=DE_max(parameters,model);
epsilons=linspace(nadirZ2,maxZ2,npar);
epsilons=epsilons(2:end-1);
Pareto_DE=zeros(npar,3);
Pareto_DE(1,:)=[nadirZ2 , minZ1 , nadirZ2];
Pareto_DE(end,:)=[maxZ2 , nadirZ1 , maxZ2];
delta=0.0001;
%% Pareto loop
itp=2; %iteration counter for pareto loop
for eps=epsilons
    UNABLE=0;
    %% Initialization
    individual.x=[];
    pop=repmat(individual,npop,1);
    %BestZ2=zeros(Maxit,1);
    for i=1:npop
        pop(i).x=CreateRandomSolution(model);
        [~ ,pop(i).Z2]=Cost2(pop(i).x,model,eps,delta);
    end
    %% DE Main loop 1
    it=1;
    while it<=Maxit && max([pop.Z2])<eps  
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
            [~,pop(i).Z2T]=Cost2(pop(i).T,model,eps,delta);
        end
        % Fitness and Selection 
        for i=1:npop
            if pop(i).Z2T>=pop(i).Z2
                pop(i).x=pop(i).T;
                pop(i).Z2=pop(i).Z2T;
            end
        end
        it=it+1;
        if it>Maxit && pop(1).Z2<eps
            UNABLE=1;
            break
        end        
    end
    if UNABLE==0
        for i=1:npop
            pop(i).Z1=Cost2(pop(i).x,model,eps,delta);
        end
        %pop=rmfield(pop,'Z2T');
        %% part 2, minimizing z1(eps. cons. objective function) while ensuring 
        % the epsilon constraint in satisfied
        it=0;
        while it<=Maxit   
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
                [pop(i).Z1T,pop(i).Z2T] =Cost2(pop(i).T,model,eps,delta);
            end
            % Fitness and Selection 
            for i=1:npop
                if pop(i).Z1T<=pop(i).Z1 && pop(i).Z2T>=eps
                    pop(i).x=pop(i).T;
                    pop(i).Z1=pop(i).Z1T;
                end
            end
        end
        [~,best]=min([pop.Z1]);
        [bestZ1,bestZ2,variables]=Cost2(pop(best).x,model,eps,delta);
        bestZ1=bestZ1+delta*variables.slack;  %transforming z1 to total cost
        Pareto_DE(itp,:)=[eps,bestZ1,bestZ2];
        format shortG
        disp([itp , Pareto_DE(itp,1) , Pareto_DE(itp,2) , Pareto_DE(itp,3)]);
        itp=itp+1;
    end
end
Pareto_DE(Pareto_DE(:,1)==0,:)=[];
time=toc;
disp('------------------------------------------')
disp(['Time = ' num2str(time)])
%% visualizing pareto frontiers
Cplex_Pareto=readmatrix('Pareto.xlsx','Sheet',TP,'Range','E3:F42');
[p_min,p_max]=pareto_dominant_minmax(Pareto_DE(:,2),Pareto_DE(:,3));
plot(p_min,p_max,'r*');
hold on
[cp_min,cp_max]=pareto_dominant_minmax(Cplex_Pareto(:,1),Cplex_Pareto(:,2));
plot(cp_min , cp_max,'bo')