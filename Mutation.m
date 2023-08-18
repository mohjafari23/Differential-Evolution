function Ans=Mutation(p1,p2,p3,F)
    Move=F*(p2-p3);
    if norm(Move)>0
        Move=Move/norm(Move);
    end
    Ans=zeros(1,length(p1));
    for i=1:length(p1)
        if Move(i)>0   
            Ans(i)=p1(i)+Move(i)*(1-p1(i));
        else
            Ans(i)=p1(i)+p1(i).*Move(i);
        end
    end       
 end
