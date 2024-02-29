classdef MOEAADR < ALGORITHM
% <multi> <real/binary> <large/none> <constrained/none> <sparse>

    methods
        function main(Algorithm,Problem)
            
            
            REAL    = ~strcmp(Problem.encoding,'binary');
            
            if REAL
                Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            else
                Dec = ones(Problem.N,Problem.D);
            end
            
            Mask = false(size(Dec));
            for i = 1 : Problem.N
                Mask(i,randperm(end,randi(end))) = true;
            end
            Population = SOLUTION(Dec.*Mask);   
            
            rho = 0.5;
            
            
            [Population,Dec,Mask,FrontNo] = ...
                    EnvironmentalSelection(Population,Dec,Mask,Problem.N,0,0);

            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo);
                
                [rbm,allOne,other,Site] = trainRBM(Mask,Dec,rho,Problem.N);

                [OffMask,OffDec] = Operator(Mask(MatingPool,:),Dec(MatingPool,:),...
                                    Site,rbm,allOne,other,Problem.lower,Problem.upper,REAL);

                Offspring = SOLUTION(OffDec.*OffMask);
                
                
                [Population,Dec,Mask,FrontNo,UpdateRatio] = ...
                    EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N,length(Population),sum(Site));
                
                
                rho = (rho+UpdateRatio)/2;
            end
        end
    end
end