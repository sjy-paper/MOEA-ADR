function [OffMask,OffDec] = Operator(ParentMask,ParentDec,Site,rbm,allOne,other,Lower,Upper,REAL)
    
    Parent1Mask = ParentMask(1:floor(end/2),:);
    Parent2Mask = ParentMask(floor(end/2)+1:end,:);
    Parent1Dec  = ParentDec(1:floor(end/2),:);
    Parent2Dec  = ParentDec(floor(end/2)+1:end,:);
    
    if any(Site)
        
        OffTemp = BinaryCrossover(rbm.reduce(Parent1Mask(Site,other)),rbm.reduce(Parent2Mask(Site,other)));
        OffTemp = rbm.recover(OffTemp);
        
        OffMaskT1 = false(size(OffTemp,1),size(Parent1Mask,2));
        OffMaskT1(:,other)  = OffTemp;
        OffMaskT1(:,allOne) = true;
        
    else
        OffMaskT1 = [];
    end
    OffMaskT2 = false(sum(~Site)*2,size(ParentMask,2));
    OffMaskT2(:,allOne) = true;
    OffMaskT2(:,other) = BinaryCrossover(Parent1Mask(~Site,other),Parent2Mask(~Site,other));
    
    
    OffMask = [OffMaskT1;OffMaskT2];
    OffMask = BinaryMutation(OffMask);
    
    if any(Site)
        if REAL

            ReduceMask = [Parent1Mask(Site,other);Parent2Mask(Site,other)];
            
            reducePos = (ReduceMask & OffTemp);
            allZero2 = all(~reducePos,1);
            allOne2  = all(reducePos,1);
            other2   = ~allZero2 & ~allOne2;
            
            OffDecTemp1 = Parent1Dec(Site,other);
            OffDecTemp2 = Parent2Dec(Site,other);
            L = Lower(:,other);
            U = Upper(:,other);
            
            OT1 = zeros(sum(Site)*2,sum(other));

            OT1(:,other2) = RealVariation(OffDecTemp1(:,other2),OffDecTemp2(:,other2),L(:,other2),U(:,other2));
            OT1(:,~other2) = RealVariation(OffDecTemp1(:,~other2),OffDecTemp2(:,~other2),L(:,~other2),U(:,~other2));

            
            OffDecT1 = zeros(sum(Site)*2,size(ParentDec,2));
            
            OffDecT1(:,other) = OT1;
            OffDecT1(:,~other) = RealVariation(Parent1Dec(Site,~other),Parent2Dec(Site,~other),Lower(:,~other),Upper(:,~other));

        end
    else
        OffDecT1 = [];
    end
   
    if REAL
        pos = (ParentMask&OffMask);
        allZero3 = all(~pos,1);
        allOne3  = all(pos,1);
        other3   = ~allZero3 & ~allOne3;

        OffDecT2 = zeros(sum(~Site)*2,size(ParentDec,2));
        OffDecT2(:,other3) = RealVariation(Parent1Dec(~Site,other3),Parent2Dec(~Site,other3),Lower(:,other3),Upper(:,other3));
        OffDecT2(:,~other3) = RealVariation(Parent1Dec(~Site,~other3),Parent2Dec(~Site,~other3),Lower(:,~other3),Upper(:,~other3));

        OffDec = [OffDecT1;OffDecT2];
    else
        OffDec = ones(size(OffMask));
    end
end





function Offspring = BinaryCrossover(Parent1,Parent2)
    [N,D] = size(Parent1);
    if N==0 || D ==0
        Offspring = [Parent1;Parent2];
        return;
    end
    
    k = repmat(randi(D,N,1),1,D) < repmat(1:D,N,1);
    
    Offspring1    = Parent1;
    Offspring2    = Parent2;
    Offspring1(k) = Parent2(k);
    Offspring2(k) = Parent1(k);
    Offspring     = [Offspring1;Offspring2];
end

function Offspring = BinaryMutation(Offspring)

    Site = rand(size(Offspring)) < 1/size(Offspring,2);
    Offspring(Site) = ~Offspring(Site);
end


function Offspring = RealVariation(Parent1,Parent2,Lower,Upper)

    [proC,disC,proM,disM] = deal(1,20,1,20);
    [N,D] = size(Parent1);
    beta  = zeros(N,D);
    mu    = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                 (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
    
    Lower = repmat(Lower,2*N,1);
    Upper = repmat(Upper,2*N,1);
    
    Site  = rand(2*N,D) < proM/D;
    mu    = rand(2*N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring       = min(max(Offspring,Lower),Upper);
end


