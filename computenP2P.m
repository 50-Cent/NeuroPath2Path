function [opFeature, noNeuronPaths, avL, pL] = computenP2P(ipData)

%% Initialization
opFeature = [];
noNeuronPaths = 0;
avL = 0;
pL = 0;

%% total sampled locations
noNodes = size(ipData,1);

%% translating the root to (0,0,0)
rootN = ipData(1,3:5);
if norm(rootN,2) == 0
    ipData(:,3:5) = ipData(:,3:5)-rootN;
end
%% alignment of dimensions based on maximum range
rangeX = max(ipData(:,3))-min(ipData(:,3));
rangeY = max(ipData(:,4))-min(ipData(:,4));
rangeZ = max(ipData(:,5))-min(ipData(:,5));

indx = [3 4 5];
[~,IX] = sort([rangeX;rangeY;rangeZ]);
indx = indx(IX);
ipData(:,3:5) = ipData(:,indx);

clear rangeX rangeY rangeZ IX indx 

%% find degree of nodes
degNeuron = zeros(noNodes,1);
for k = 2:noNodes
   degNeuron(ipData(k,1)) = degNeuron(ipData(k,1))+1;
   degNeuron(ipData(k,7)) = degNeuron(ipData(k,7))+1;
end
degNeuron(1) = sum(ipData(:,7)==1);

clear k

%% list of bifurcation nodes
nodeDendrite = (ipData(:,2)==3)+ (ipData(:,2)==4); 
nodeDegUnity = (degNeuron==1);
dendriteEnd = find((nodeDegUnity.*nodeDendrite)==1);
if isempty(dendriteEnd)                                 % It is axon file
    return
end
if dendriteEnd(1)==1
    dendriteEnd(1)=[];
end
noNeuronPaths = length(dendriteEnd);

clear nodeDendrite nodeDegUnity 

%% find bifurcations and multifurcations
bifurList = find(degNeuron(2:end)>=3)+1; 

%% path concurrence
concurList = zeros(noNodes,1);
concurList(dendriteEnd)=1;
for k = 1:size(dendriteEnd,1)
    pvt = dendriteEnd(k);
    concurList = pathConcurrence(ipData,concurList,pvt);
end
concurList(1)= noNeuronPaths;

clear pvt k

%% find bifurcation angle + partition asymmetry
bifurAngle = zeros(length(bifurList),1);
cauL = zeros(length(bifurList),1);

for k = 1:length(bifurList)
   pvt = bifurList(k);
   loc = find(ipData(:,7)==pvt);
   maxAngle = 0;
   maxAsym = 0;
   for m =1:length(loc)-1
       uu = ipData(loc(m),3:5)-ipData(pvt,3:5);
       gg = concurList(loc(m));
       for n=m+1:length(loc)
            vv = ipData(loc(n),3:5)-ipData(pvt,3:5);
            %thetaD = abs(atan2d(norm(cross(uu,vv)),dot(uu,vv)));
            thetaD = abs(acos(dot(uu,vv)/(norm(uu,2)*norm(vv,2))));
            hh = concurList(loc(n));
            asym = abs(gg-hh)/(gg+hh);
            if thetaD > maxAngle
                maxAngle = thetaD;
            end
            if asym > maxAsym
                maxAsym = asym;
            end
       end
   end
   bifurAngle(k) = maxAngle;
   cauL(k)=maxAsym;
end

% find bifurcation at the root.
pvt = 1;
loc = find(ipData(:,7)==pvt);
maxAngle = 0;
maxAsym = 0;
for m =1:length(loc)-1
   uu = ipData(loc(m),3:5)-ipData(pvt,3:5);
   gg = concurList(loc(m));
   for n=m+1:length(loc)
        vv = ipData(loc(n),3:5)-ipData(pvt,3:5);
        thetaD = abs(acos(dot(uu,vv)/(norm(uu,2)*norm(vv,2))));
        hh = concurList(loc(n));
        asym = abs(gg-hh)/(gg+hh);
        if thetaD > maxAngle
            maxAngle = thetaD;
        end
        if asym > maxAsym
            maxAsym = asym;
        end
   end
end
bifurRoot = maxAngle;
bifurcauL = maxAsym;

clear m k uu vv thetaD maxAngle k pvt loc asym gg hh maxAsym 


%% hierarchy can be computed from Concurrence

%% List of paths with location indices
tmp = [];
enumeratePath = cell(noNeuronPaths,1);
for k = 1:noNeuronPaths
    pvt = dendriteEnd(k);
    while pvt ~=1
        tmp = [tmp; pvt];
        pvt = ipData(pvt,7);
    end
    tmp = [tmp; 1];   
    enumeratePath{k} = flipud(tmp);         % each path in enumerate starts with node 1 at top    
    tmp = [];
end

clear k pvt tmp 

%% competitive proximity
delT = 0.5;
comProx = zeros(length(bifurList),1);
for k = 1:length(bifurList)
   pvt = bifurList(k);
   tmpData = ipData(:,3:5);
   
   %%%%%%%%%%%%%%%% standardization parameters%%%%%%%%%
   muC = mean(tmpData,1);
   sdC = std(tmpData,1); 
   tmpData(:,1) = (tmpData(:,1)-muC(1))./sdC(1);
   tmpData(:,2) = (tmpData(:,2)-muC(2))./sdC(2);
   tmpData(:,3) = (tmpData(:,3)-muC(3))./sdC(3);
   %%%%%%%%%%%%%%%%
   
   pvtC = tmpData(pvt,:);
   for m = 1:length(enumeratePath)
      tmp = enumeratePath{m};
      if isempty(intersect(tmp,pvt))
          dsT = sqrt((tmpData(tmp,1)-pvtC(1)).^2+(tmpData(tmp,2)-pvtC(2)).^2+(tmpData(tmp,3)-pvtC(3)).^2);
          if any(dsT<=delT)              
              comProx(k) = comProx(k)+1;
          end
      end      
   end  
end
comProx = comProx./noNeuronPaths;
clear dsT m k node2del tmpData muC sdC pvt

%% feature [EACH PATH]
%% tortuosity+bifurcation-angle+partition-asymmetry+concurrence+hierarchy+length+curvature+proximity 
tortuositY     = cell(noNeuronPaths,1);
bifurAnglePath = cell(noNeuronPaths,1);
caulescencE    = cell(noNeuronPaths,1);
concurrencE    = cell(noNeuronPaths,1);
hierarchY      = cell(noNeuronPaths,1);
seglengtH      = cell(noNeuronPaths,1);
competitioN    = cell(noNeuronPaths,1);


for k = 1:noNeuronPaths
   tmp = enumeratePath{k};  
   [B,ia,ib] = intersect(bifurList,tmp);
   bifurPerpath = [1;B;tmp(end)];
   tort = zeros(length(bifurPerpath)-1,1); 
   sgL = zeros(length(bifurPerpath)-1,1);
   %avLval = zeros(length(bifurPerpath)-1,1);
   
   for m=1:length(bifurPerpath)-1
      fI = find(tmp == bifurPerpath(m));
      fS = find(tmp == bifurPerpath(m+1));
      if fS == fI+1
          tort(m)=1;
          sgL(m) = norm(ipData(bifurPerpath(m),3:5)-ipData(bifurPerpath(m+1),3:5),2);
          %avLval = sgL(m);
      else
          segmenT = tmp(fI:fS);
          eucliD = norm(ipData(bifurPerpath(m),3:5)-ipData(bifurPerpath(m+1),3:5),2);
          coord = ipData(segmenT,3:5);
          coord = diff(coord,1);         
          dsT = sum(sqrt(coord(:,1).^2+coord(:,2).^2+coord(:,3).^2));
          tort(m) = dsT/eucliD;
          sgL(m) = dsT;
          %avLval(m) = sgL(m)/concurList(tmp(ib));
      end      
   end
   tortuositY{k} = tort;
   bifurAnglePath{k} = [bifurRoot;bifurAngle(ia)];
   caulescencE{k} = [bifurcauL;cauL(ia)];
   concurrencE{k} = [concurList(1);concurList(tmp(ib))];
   [~,IIX]   = sort([concurList(1);concurList(tmp(ib))],'descend');
   hierarchY{k} = IIX;
   seglengtH{k}   = sgL;
   competitioN{k} = [0;comProx(ia)];
end

clear k tort segmenT eucliD coord fS fI m dsT bifurPerpath ia ib sgL IIX
clear bifurAngle cauL concurList


%% Average path length computation
avL = 0;
for k = 1:noNeuronPaths
    avL = avL+sum(seglengtH{k});
end
avL = avL/noNeuronPaths;

%% order span calculation

%% pL calculation
pL = [];
for k = 1:noNeuronPaths
    pL = [pL;max(hierarchY{k})];
end
pL = max(pL);
%% feature assembly
opFeature = [tortuositY bifurAnglePath caulescencE concurrencE hierarchY seglengtH competitioN enumeratePath];


end