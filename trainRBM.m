function [rbm,allOne,other,Site] = trainRBM(Mask,Dec,UpdateRatio,N)

    allZero = all(~Mask,1);
    allOne  = all(Mask,1);
    other   = ~allZero & ~allOne;

    Site = false(1,floor(N/2));
    
    rbm = [];
    
    if rand < UpdateRatio
        
        Site(randperm(floor(N/2),randi(ceil(N/2*0.1)))) = true;
        
        K = ceil(sum(mean(abs(Mask(:,other).*Dec(:,other))>1e-6,1)));
        
        K = min(K,ceil(size(Mask,2)*0.2));
        
        rbm = RBM(sum(other),K,8,1,0,0.4,0.15);
        
        rbm.train(Mask(:,other));

    end
end

