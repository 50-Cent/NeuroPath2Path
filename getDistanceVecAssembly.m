function [M1, M2] = getDistanceVecAssembly(listnP2P, indxx)

fDim = 6;
M1 = zeros(fDim,fDim);
M2 = zeros(fDim,fDim,length(indxx));
%% Interclass [for parallel processing use parfor]
for k = 1:length(indxx)-1
    dataC1 = listnP2P{k};
    for l = 1:length(dataC1)              
        data1 = dataC1{l}; 
        for m = k+1:length(indxx)
            dataC2 = listnP2P{m};
            for h = 1:length(dataC2)
            %[l h]
               data2 = dataC2{h};
               [~,dstVec,~] = computeDST(data1,data2);
               for p = 1:size(dstVec,2)
                    M1 = M1+dstVec(:,p)*dstVec(:,p)';
               end
            end                
        end
    end
end

for k = 1:length(indxx)
    dataC1 = listnP2P{k};
    for m = 1:length(dataC1)-1
        data1 = dataC1{m};
        for n = m+1:length(dataC1)
            data2 = dataC1{n};
            [~,dstVec,~] = computeDST(data1,data2);
            for p = 1:size(dstVec,2)
                    M2(:,:,k) = M2(:,:,k)+dstVec(:,p)*dstVec(:,p)';
            end
        end
    end        
end

clear dstVec k m n p h data1 data2 l dataC1 dataC2 indxx 

end