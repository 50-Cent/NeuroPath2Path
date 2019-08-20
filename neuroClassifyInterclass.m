function [] = neuroClassifyInterclass(testProp,delS,delA,majorVote,spn)

%% Description
% testprop = ratio of the test data. testprop = 0.1 means that 10% of the
% data are tested
% delS = binwidth of the distribution of paths. 
% delA = search step to find nearest candidates of a test sample (nearest in terms of paths)
% majorVote = majority voting parameter to determine the class
% spn = limiting factor to stop searching for candidates. 
%
% Ex: neuroClassifyInterClass(0.1,15,5,7,5). In majority of the
% experiments, delS, spn, and delA are kept fixed. User can opt other
% values by trials



%tempPath = strsplit(pwd,'/');
%newPath = strcat(strjoin(tempPath(1:end-1),'/'),'/','Dataset','/','Cell_type');
fList = dir('______Enter Data Path Here_______');

maxorderList = zeros(length(fList)-2,1);
for k = 3:length(fList)    
    folder_name = fList(k).name
    basePath = strcat(newPath,'/',folder_name);
    neuronData = readInput(basePath);    
    sz = length(neuronData); 
    pathCount = zeros(sz,1);                    % path length
    avLength = zeros(sz,1);                     % average path length 
    pathorder = [];
    ctt = 1;
    for kk = 1:sz   
        [fT,num,avL,pL] = computenP2P(cell2mat(neuronData{kk}));        
        if ~isempty(fT)                         %% delete files containing axon only
            listnP2P{k-2}{ctt} = fT;
            dataStore{k-2}{ctt} = cell2mat(neuronData{kk});
            pathCount(ctt) = num;
            avLength(ctt) = avL;
            pathorder = [pathorder;pL];
            ctt = ctt+1;
        end        
    end 
    IXI = find(pathCount==0);
    pathCount(IXI)=[];
    avLength(IXI) = [];
    pathList{k-2}       = pathCount;                %storing number of paths in each neuron
    avList{k-2}         = avLength;                 %
    maxorderList(k-2) = max(pathorder);
end
clear count fList  k kk sz folder_name pathCount IXX fT num avLength avL pL pathlength ctt
clear basePath newPath tmpPath currPath neuronData

%% dataset partition
%testProp = 0.1;                          %testing parameter
nmm = randi([1 10000000],1);             %to partition the data randomly
rng(nmm)
no_class = length(listnP2P);
%delS = 15;                               %bin width of the distribution of path sizes: user parameter
trainIDXwL = cell(no_class,1);           
testIDXwL  = cell(no_class,1);
 
testdataperclass = zeros(no_class,1);
for k = 1:no_class
   noData = length(listnP2P{k}); 
   pathL = pathList{k};
   idx = [1:length(pathL)]';   
   [histL, histEdge] = histcounts(pathL,'BinWidth',delS);   
   histLL = histL/length(pathL);                         % probability (skewed and heavy-tailed)in delS bins    
   traindataprop = floor(histLL*(1-testProp)*noData) ;      % divided based on the distribution of paths 
   for kk = 1:length(traindataprop)
      if histL(kk)==1
          traindataprop(kk) = 1;
      end
   end  
   trainIDX = [];
   testIDX = [];
   for kk = 1:length(traindataprop)
       if traindataprop(kk)~=0
          id = (pathL>=histEdge(kk))&(pathL<histEdge(kk+1));      
          id = idx(id);          
          rndP = randperm(length(id));          
          trainIDX = [trainIDX;id(rndP(1:traindataprop(kk)))];          
          testIDX = [testIDX; id(rndP(traindataprop(kk)+1:end))];          
       end
   end
   testdataperclass(k) = length(testIDX);
   trainIDXwL{k} = [trainIDX pathL(trainIDX) avList{k}(trainIDX)];  %per class training set indices with no of paths
   testIDXwL{k}  = [testIDX pathL(testIDX) avList{k}(testIDX)];    %per class test set indices with no of paths
end

clear id trainIDX testIDX kk k mn noData pathL rndP notraindata traindataprop histL histEdge idx id 

disp('<----- Unsupervised Classification  :: nearest neighbor-->')
%% classification [UNSUPERVISED]
%delA = 5;
predictedLabel=cell(no_class,1);
trueLabel=cell(no_class,1);
%majorVote = 5;
for k=1:no_class
    predictedLabel{k}= zeros(testdataperclass(k),1);
    trueLabel{k}=k*ones(testdataperclass(k),1);
end
% 
for k = 1:no_class
    k
    tmp = testIDXwL{k};
    dataA = listnP2P{k};    
    for m = 1:size(tmp,1)       
         pvt = tmp(m,1);
         pvtData = dataA{pvt};
         sZ = tmp(m,2);            %no of paths
         nL = tmp(m,3);
         dstArray = [];
         labelArray = [];
         for n = 1:no_class
             dataC = listnP2P{n};      
             trPathL = trainIDXwL{n}(:,2);                            % only the number of paths
             noCandidate = 0;                     
             fct = 1;
             while noCandidate < majorVote
                 noCandidate = sum((trPathL > sZ-fct*delA) & (trPathL < sZ+fct*delA));
                 fct = fct+1;
                 if fct == spn
                     break;
                 end
             end
             fct = fct-1; 
             trIndx = trainIDXwL{n}(:,1);
             candidateNeuron = trIndx((trPathL > sZ-fct*delA) & (trPathL < sZ+fct*delA));    
             
             if ~isempty(candidateNeuron)             
               for p = 1:length(candidateNeuron)
                  [dsT,~,~] = computeDSTInter(dataC{candidateNeuron(p)},pvtData); 
                  dstArray = [dstArray;dsT];
                  labelArray = [labelArray;n];
               end
             end
             
         end         
         [~,IX] = sort(dstArray);
         predictedLabel{k}(m,1) = mode(labelArray(IX(1:majorVote)));
    end 
    disp(strcat('_____ END :: class',num2str(k),'-','test labeling ______'))
end

clear  candidateNeuron dataA dsT labelArray  IX k m n p pathListpvt pvtData sZ sZfromtraining dataC
clear tmp 
%% performance
tLabel = [];
pLabel = [];
for k = 1:no_class
   tLabel = [tLabel;trueLabel{k}];
   pLabel = [pLabel;predictedLabel{k}];
end

CP = confusionmat(tLabel,pLabel)
accuracY = trace(CP)/sum(sum(CP))*100

end

