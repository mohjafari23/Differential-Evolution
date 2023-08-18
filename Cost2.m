function [z1, z2, variables , epsc]=Cost2(X,model,eps,delta)
% cost based on two catagory of random varables:
% X_p: for priority and X_a: for allocation
I=model.I;
S=model.S;
E=model.E;
L=model.L;
F=model.F;
D=model.D;
D0=D;
CAP=model.CAP;
CAP0=CAP;
CAPb=model.CAPb;
CAPi=model.CAPi;
PA=model.PA;
SA=model.SA;
G=model.G;
FC=model.FC;
FCB=model.FCB;
CI=model.CI;
CF=model.CF;
CP=model.CP;
Q=model.Q;
teta=model.teta;
beta=model.beta;
TC=model.TC;
TCw=model.TCw;
CPB=model.CPB;
CE=model.CE;
CW=model.CW;
CB=model.CB;
ppf=model.ppf;
SS=model.SS;
moq=model.moq;
bud=model.bud;
budf=model.budf;

BS=3; % number of strategies for covering shortage

Xp=reshape(X(1:I*S),I,S); %for xp, supplier priority for first stage alocation
BS_PR=reshape(X(I*S+1:I*S+I*E*BS),I,E,BS); %i e BS , for determining the priority of bs
BS_sr=reshape(X(I*S+I*E*BS+1:I*S+I*E*BS+I*BS),I,BS); % i  bs , for selecting supplier for bs
fr=X(I*S+I*E*BS+I*BS+1:I*S+I*E*BS+I*BS+1); % just one , for fortification
 
xr=reshape(X(I*S+I*E*BS+I*BS+1+1:I*S*2+I*E*BS+I*BS+1),I,S); %i s
%or=reshape(X(I*S*2+I*E*BS*2+1+1:I*S*3+I*E*BS*2+1),I,S); % i s
%br=reshape(X(I*S*3+I*E*BS*2+1+1:I*S*3+I*E*BS*2+1+I*E),I,E); %i e
bsxr=reshape(X(I*S*2+I*E*BS+I*BS+1+1:I*S*2+I*E*BS+I*BS+1+I*S*BS),I,S,BS); % i s  bs , for allocating eu , n , w
br=X(I*S*2+I*E*BS+I*BS+1+I*S*BS+1:I*S*2+I*E*BS+I*BS+1+I*S*BS+I);
% order allocation (Xis , z , m) -----------------------------------------
x=zeros(I,S);
CAP=floor(CAP.*(1./(ones(I,S)+G)));
[~, xp]= sort(Xp,2); %determining pririty of allocation
xpc=cell(I,1);
for i=1:I
    xpc{i}=xp(i,:);
end

% there is a xr(i,s) random matrix for quantity allocation
% CAP=CAP.*PA;
for i = 1:I
    vs=zeros(1,S);
    % first round of allocation
    for s=xpc{i}
        vs(s)=1; % visited suppliers 
        if D(i)>=moq(i,s)
            if PA(i,s) > 0
                if CAP(i,s)>=D(i) || nnz(vs)==length(vs)
                    x(i,s)=round(moq(i,s)+xr(i,s)*(D(i)-moq(i,s)));
                    D(i)=D(i)-x(i,s);
                    CAP(i,s)=CAP(i,s)-x(i,s);
                else
                    x(i,s)=round(moq(i,s)+xr(i,s)*(min(CAP(i,s),(D(i)-max(moq(i,vs==0))))-moq(i,s)));
                    D(i)=D(i)-x(i,s);
                    CAP(i,s)=CAP(i,s)-x(i,s);
                end
            end
            if D(i)==0
                break
            end
        else
           xpc{i}(xpc{i}==s)=-1;
        end
    end
    xpc{i}(xpc{i}==-1)=[];
    % second (or more) round of allocation
    while D(i)>0
        for s=xpc{i}
            if PA(i,s) > 0
                xx = ceil(xr(i,s)*(min(D(i),CAP(i,s))));
                x(i,s)= x(i,s) + xx;
                D(i)=D(i)-xx;
                CAP(i,s)=CAP(i,s)-xx;
            end
            if D(i)== 0
                break
            end
        end
    end
end

z=zeros(1,S);
for s = 1:S
    for i =1:I
        if x(i,s)>0
           z(1,s)=1;
           break
        end
    end
end
m=zeros(I,S,L);
for s = 1:S
    for i =1:I
        for l =1:L
            if x(i,s)<=Q(i,s,l)
                m(i,s,l)=1;
                break
            end
        end
    end
end
x_l=reshape(repmat(x,1,L),I,S,L);
for i=1:I
    for s=1:S
        x_l(i,s,m(i,s,:)==0)=0;
    end
end

% pre-positioning (w) ---------------------------------------------------
% CAPi(:,z==0)=0;
% o=round(or.*(CAPi.*PA));       

% second stage variables -------------------------------------------------
% calculating amount of recieved in each scenario (uu), before any strategy
x_e=reshape(repmat(x,1,E),I,S,E);
SA_e=permute(reshape(repelem(SA,1,I),S,I,E),[2 1 3]);
uu=x_e.*SA_e;
b=reshape(repmat(reshape(D0,I,1),1,E),I,E)-reshape(sum(uu,2),I,E); %shortage without strategies

% compensetaing strategies for shortage:
% BS=1 : eu | BS=2 : n | BS=3 : w
% capacities for covering shortage:
CAPfb=zeros(I,S,E,BS);

G_e=reshape(repmat(G,1,E),I,S,E);
%CAP_e=reshape(repmat(CAP,1,E),I,S,E);
CAPfb(:,:,:,1)=floor(uu.*G_e);

CAPb(:,z==1,:)=0;
CAPb=CAPb.*PA;
CAPfb(:,:,:,2)=reshape(repmat(CAPb,1,E),I,S,E);
CAPfb(:,:,:,2)=CAPfb(:,:,:,2).*SA_e;

CAPi(:,z==0)=0;
CAPi=CAPi.*PA;
CAPfb(:,:,:,3)=reshape(repmat(CAPi,1,E),I,S,E);

% CAPfb_s=reshape(sum(sum(CAPfb,2),4),I,E);

% b=round((bb-min(bb,CAPfb_s))+br.*min(bb,CAPfb_s)); %shortage
bsx=zeros(I,S,E,BS);  %amount of allocation to each shortage covering strategy
BS_pr=zeros(I,E,BS); % for detrmining priority of each bs for each i in each e 
for i=1:I
    for e=1:E
        [~ , BS_pr(i,e,:)]=sort(BS_PR(i,e,:));
    end
end

for i=1:I
    for e=1:E
        if b(i,e)>0
            for bs=BS_pr(i,e,:)
                if nnz(CAPfb(i,:,e,bs))>0
                    nzs=find(CAPfb(i,:,e,bs)>0); %nzs: suppliers with none zero capacity of bs
                    BS_ss=round(1+BS_sr(i,bs)*(length(nzs)-1));
                    if bsxr(i,nzs(BS_ss),bs)>=br(i)
                        bsx(i,nzs(BS_ss),e,bs)= min(CAPfb(i,nzs(BS_ss),e,bs),b(i,e));
                        CAPfb(i,nzs(BS_ss),e,bs)=CAPfb(i,nzs(BS_ss),e,bs)-bsx(i,nzs(BS_ss),e,bs);
                        b(i,e)=b(i,e)-bsx(i,nzs(BS_ss),e,bs);
                    end
                    if b(i,e)==0
                        break
                    end
                end
            end
        end
    end
end
eu=bsx(:,:,:,1);
n=bsx(:,:,:,2);        
w=bsx(:,:,:,3);
u=uu+eu;

y=zeros(1,S);
for i=1:I
    for e=1:E
        for s=1:S    
            if n(i,s,e)>0
                y(1,s)=1;
                break
            end
        end
    end
end
o=max(w,[],3);
% fortification -----------------------------
fs_f=find(CF<=budf);  %feasable fortification strategies
fs=fs_f(round(1+fr*(length(fs_f)-1)));
% Costs ----------------------------------------------------------------
FiC=sum(FC.*z)+sum(FCB.*y);
PPC=sum(CI.*o,'all');
FoC=CF(fs);

teta=reshape(repmat(repmat(reshape(teta,I,1),1,S),1,L),I,S,L);
OC1=CP.*x_l+teta.*CP.*x_l;
OC1=reshape(repmat(sum(sum(OC1,3),1),1,E),S,E);
OC1=sum(OC1.*SA,1);

beta=reshape(repmat(repmat(reshape(beta,I,1),1,S),1,L),I,S,L);
OC2=beta.*(CP.*x_l);
OC2=reshape(repmat(sum(sum(OC2,3),1),1,E),S,E);
OC2=sum(OC2.*(1-SA),1);

OC=OC1+OC2;

TC=reshape(repmat(TC,1,E),I,S,E);
TCw=reshape(repmat(TCw,1,E),I,S,E);
TrC=reshape(sum(sum((TC.*(n+u)+TCw.*w),1),2),1,E);

CPB=reshape(repmat(CPB,1,E),I,S,E);
CE=reshape(repmat(CE,1,E),I,S,E);
CW=reshape(repmat(CW,1,E),I,S,E);
CB=reshape(repmat(CB,1,E),I,E);
RC=reshape(sum(sum((CPB.*n+CE.*eu+CW.*w),1),2),1,E)+sum(CB.*b,1);

p=ppf(:,fs);
TotalCost=FiC+PPC+FoC+(OC+TrC+RC)*p;

SS=reshape(repmat(repmat(reshape(SS,1,S),I,1),1,E),I,S,E);
z2=(reshape(sum(sum((SS.*(u+n+w)),1),2),1,E))*p;

% the epsilon constraint -----------------------------------------------
% when z2>=eps there is a slack>0 such that z2-slack=eps
% so having z2>=eps satisfies the epsilon constraint
epsc=z2>=eps; % true if is satisfied
slack=z2-eps;
z1=TotalCost-delta*slack;

variables.x=x;
variables.u=u;
variables.eu=eu;
variables.n=n;
variables.w=w;
variables.b=b;
variables.o=o;
variables.z=z;
variables.y=y;
variables.fs=fs;
variables.slack=slack;




