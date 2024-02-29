function [Population,Dec,Mask,FrontNo,UpdateRatio] = EnvironmentalSelection(Population,Dec,Mask,N,len,num)
    success = false(1,length(Population));

    
    [~,uni] = unique(Population.objs,'rows');
    if length(uni) == 1
        [~,uni] = unique(Population.decs,'rows');
    end
    Population = Population(uni);
    Dec        = Dec(uni,:);
    Mask       = Mask(uni,:);
    N          = min(N,length(Population));

    PopObj = Population.objs;
    [FrontNo,MaxFNo] = NDSort(PopObj,Population.cons,N);

    Next = FrontNo <= MaxFNo;

    Last     = find(FrontNo==MaxFNo);
    Del = Truncation(PopObj(Last,:),sum(Next)-N);
    Next(Last(Del)) = false;
    
    success(uni(Next)) = true;
    s1     = sum(success(len+1:min(len+num,size(success,2))));
    s2     = sum(success(len+num+1:end));

    
    UpdateRatio = (s1+s2+1e-6)/(size(success,2)-len+1e-6);
    
    if num > 1e-6  && s1 > 1e-6
        UpdateRatio = 0.5 + s1/(s1+s2);
        
    end
    
    
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    Dec        = Dec(Next,:);
    Mask       = Mask(Next,:);
    
end


function Del = Truncation(PopObj,K)

    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end
