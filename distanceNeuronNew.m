function [dsT,pairList] = distanceNeuronNew(distanceMat,data1,data2)

dsT = 0;
[M,N] = size(distanceMat);
statusSwap = 0;
if M>N
    distanceMat = distanceMat';    
    statusSwap = 1;
    tmp = data1;
    data1 = data2;
    data2 = tmp;
end
[M,N] = size(distanceMat);          % M < N always : |data1| = M, |data2| = N (after swap if applied)



kSimilarity = floor(N/M);

pathIndx = [1:N]';
dstSum = 0;
dstMatcopy = distanceMat;
listPairwVal = [];
%% append zero rows (M workers + (N-M) dummy :: N jobs)
for k = 1:kSimilarity    
    dstMatcopy = [dstMatcopy;zeros(length(pathIndx)-M,length(pathIndx))];
    [pairL, dst] = munkres(dstMatcopy); 
    partiallyAssigned = find(pairL(1:M)==0);                     %Munkres can sometimes give assigment 0 (partial assignment)
    if ~isempty(partiallyAssigned)
        unAssigned = [pairL(M+1:end) setdiff(1:length(pathIndx),pairL)];
        for ss = 1:length(partiallyAssigned)
            rw = distanceMat(partiallyAssigned(ss),:);
            [~,candidateP] = min(rw(unAssigned));               % minimum from the unassigned path
            pairL(partiallyAssigned(ss)) = candidateP(1);   % if there are two minimum, take the first one 
        end
    end
       
    %dstSum = dstSum+dst;                            % hold this distance for further analysis
    
    pairIndx = pathIndx(pairL(1:M));
    inD = sub2ind(size(distanceMat),(1:M)',pairIndx);
    dstSum = dstSum+sum(distanceMat(inD));
    listPairwVal = [listPairwVal;(1:M)' pairIndx distanceMat(inD)];
    pathIndx = setdiff(pathIndx,pairIndx);          % these paths in neuron2 NOT assigned in current iter
    
    dstMatcopy = distanceMat(:,pathIndx);    
end
clear dst dstMatcopy k pairL inD pairIndx

%% remaining (N-M) paths assignment + incomplete kSimilarity (fractional)
% we have pathIndx as a set of indices that are left out from neuron2.
% as N-kSimilarity*M < M ==> flip dstMat becuase neuron2 (N-kSimilarity*M) is smaller than
% neuron 1 (M). 

if M > 1
  dstMatcopy = distanceMat'; 
  dstMatcopy = dstMatcopy(pathIndx,:);
  dstMatcopy = [dstMatcopy; zeros(M-length(pathIndx),M)];
  [pairL, dst] = munkres(dstMatcopy);     % this pairL is for neuron1 
  dstSum = dstSum+dst;
  pairIndx = pairL(1:length(pathIndx));
  inD = sub2ind(size(distanceMat),pairIndx',pathIndx);
  listPairwVal = [listPairwVal; pairIndx' pathIndx distanceMat(inD)];
end

clear dst inD pairIndx pathIndx



% dstSum/(M*kSimilarity)

%% Check the degree of each path in the neuron with no of paths N.
[degNeuron2,~] = histcounts(listPairwVal(:,2),N);
[degNeuron1,~] = histcounts(listPairwVal(:,1),M);

%% check the individual values per pair of paths from listPairwVal M<N 
% we take skewness measure

misalignBucket = [];

mD = median(listPairwVal(:,3));
sD = std(listPairwVal(:,3));
sk = skewness(listPairwVal(:,3));
IX = [];
if sk > 0
    IX = find(listPairwVal(:,3) >= mD+sD) ;      
end

listPairwValn = listPairwVal;
listIX = [];
if ~isempty(IX)
    for k = 1:length(IX)
       pvt = listPairwVal(IX(k),1:2);
       hier1 = data1{pvt(1),5};
       hier2 = data2{pvt(2),5};
       if abs(length(hier1)-length(hier2)) < max(length(hier1),length(hier2))/2   %hierarchy mismatch
           if degNeuron1(pvt(1)) > 1
               misalignBucket = [misalignBucket;pvt(2)];
               dstSum = dstSum-listPairwVal(IX(k),3);                                 %delete the distances of pair
               listIX = [listIX;IX(k)];
           end
       end
    end
end


if length(listIX)~=0
    listPairwValn(listIX,:) = []; 
end



clear k pvt hier1 hier2 
%% reassignment to misaligned pair

if ~isempty(misalignBucket)
    for k = 1:length(misalignBucket)
        pvt = misalignBucket(k);
        dstCol = distanceMat(:,pvt);
        [num, id] = min(dstCol);
        listPairwValn = [listPairwValn; id pvt num];
        dstSum = dstSum+num; 
    end
end

dsT = dstSum/N;

clear misalignBucket id pvt num dstCol
%% Final pairlist
%size(pairList)
if statusSwap == 0
    pairList = listPairwValn(:,1:2);
else
    pairList = [listPairwValn(:,2) listPairwValn(:,1)];
end

%% extension if we allow noncommutative measure 

misalignBucket = setdiff([1:M]',listPairwValn(:,1));


%% For best match only [not part of neuropath2path]

% lstX = [];
% for kk = 1: size(distanceMat,2)
%     [~,ZZ] = min(distanceMat(:,kk));
%     lstX = [lstX;ZZ(1) kk]; 
% end
% 
% lstX
end
