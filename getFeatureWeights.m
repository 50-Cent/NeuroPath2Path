function fcoeff = getFeatureWeights(interClassF, intraClassF)

fDim = 6;
no_class = size(intraClassF,3);                 %number of classes
lagrange_coeff = ones(no_class,1);              %lagrange coefficient per class
barrier_coeff = 1;                              %log barrier coefficient
prob_coeff    = 1;                              %probability estimate coefficient   

rng(100);
rnum = randi([1,10],fDim,1);
tmp_coeff = rnum/sum(rnum);
curr_coeff = ones(fDim,1)/fDim;                 %current coefficient estimate

stp = 0.0001;                                      %step length of gradient descent              
tol = 0.01;                                    %convergence criteria tolerance

count = 1;
    while count < 2
        %for iter = 1:10
            while (norm(tmp_coeff-curr_coeff,2) > tol) 
                sm = 0;
                for k = 1:no_class
                    sm = sm+lagrange_coeff(k)*intraClassF(:,:,k);
                end
                df_dcoeff = -interClassF*tmp_coeff + sm*tmp_coeff ...
                            -barrier_coeff*(1./(tmp_coeff+0.0001));

                curr_coeff = tmp_coeff;
                curr_coeff = curr_coeff./sum(curr_coeff);
                tmp_coeff = tmp_coeff-stp*df_dcoeff;   
                tmp_coeff = tmp_coeff./sum(tmp_coeff)        
            end
            
        % end    
        count = count+1
        barrier_coeff = barrier_coeff/2;
        prob_coeff = prob_coeff*2;
        lagrange_coeff = 1.1*lagrange_coeff;               % more intraclass compaction
        tol = tol/2;
    end

fcoeff = tmp_coeff;

end