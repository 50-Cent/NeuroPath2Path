function []=jobIW()

%% This script is used to obtain the importance weight of each feature.


tempPath = strsplit(pwd,'/');
newPath = strcat(strjoin(tempPath(1:end-1),'/'),'/','Dataset','/','Cell_overall');
fList = dir(newPath);


maxorderList = zeros(length(fList)-2,1);
for k = 3:length(fList)    
    folder_name = fList(k).name
    basePath = strcat(newPath,'/',folder_name);
    neuronData = readInput(basePath);    
    sz = length(neuronData) 
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

%% Optimization
nameNeuron = ["Ganglion","Granule","Motor","Purkinje","Pyramidal"];
noClass = length(listnP2P);

parfor k = 1:noClass-1
    disp(nameNeuron(k))
    for m = k+1:noClass
    disp(nameNeuron(m))
        indxx = [k;m];
        [interClassF, intraClassF] = getDistanceVecAssembly(listnP2P, indxx);
        fcoeff = getFeatureWeights(interClassF, intraClassF);
        disp(strcat(nameNeuron(k),'----',nameNeuron(m)))
        fcoeff
    end
end

clear k m noClass nameNeuron fcoeff intraClassF interClassF indxx 

end


