function [dsTNeuron,dstVecforPair,prL] = computeDSTInter(dataA,dataB)

fDim = 6;
dsTNeuron = 0;
noPath1 = size(dataA,1);
noPath2 = size(dataB,1);
dstCell = cell(noPath1,noPath2);           % Noting distance between a pair of paths

distanceMat = zeros(noPath1,noPath2);


for k = 1:noPath1
    tort1 = dataA{k,1};
    bifur1 = dataA{k,2};
    caul1 = dataA{k,3};
    conc1 = dataA{k,4};    
    bOrder1 = dataA{k,5} ;                  %branch order = hierarchy (more hierarchy = less important)
    bOrder1 = 1./(bOrder1+0.1);    
    bOrder1 = bOrder1./sum(bOrder1);    
    segL1 = dataA{k,6};
    compt1 = dataA{k,7};
    
    for m = 1:noPath2
        tort2 = dataB{m,1};
        bifur2 = dataB{m,2};
        caul2 = dataB{m,3};
        conc2 = dataB{m,4};
        bOrder2 = dataB{m,5};
        bOrder2 = 1./(bOrder2+0.1);
        bOrder2 = bOrder2./sum(bOrder2);        
        segL2 = dataB{m,6};
        compt2 = dataB{m,7};
        
        %% standard branch order       
       
        %tortuosity        
        maxTort = max([tort1;tort2])+0.0001;
        minTort = min([tort1;tort2]);
        tort11 = (tort1-minTort)/(maxTort-minTort);
        tort2 = (tort2-minTort)/(maxTort-minTort);
        if length(tort11) >= length(tort2)
            tort2 = [tort2;zeros(length(tort11)-length(tort2),1)] ;           
            dstTort = norm((tort1-tort2).*bOrder1,2);           
        elseif length(tort2) > length(tort1)           
            tort11 = [tort11;zeros(length(tort2)-length(tort11),1)];
            dstTort = norm((tort11-tort2).*bOrder2,2);             
        end
        
        
        %bifurAngle
        maxBifur = max([bifur1;bifur2])+0.0001;
        minBifur = min([bifur1;bifur2]);
        bifur11 = (bifur1-minBifur)/(maxBifur-minBifur);
        bifur2 = (bifur2-minBifur)/(maxBifur-minBifur);
        if length(bifur11) >= length(bifur2)
            bifur2 = [bifur2;zeros(length(bifur11)-length(bifur2),1)];
            dstBifur = norm((bifur11-bifur2).*bOrder1,2);            
        elseif length(bifur2) > length(bifur11)          
            bifur11 = [bifur11;zeros(length(bifur2)-length(bifur11),1)];
            dstBifur = norm((bifur11-bifur2).*bOrder2,2);
        end
        
        
        %caulescence
        maxCaul = max([caul1;caul2])+0.0001;
        minCaul = min([caul1;caul2]);
        caul11 = (caul1-minCaul)/(maxCaul-minCaul);
        caul2 = (caul2-minCaul)/(maxCaul-minCaul);
        if length(caul11) >= length(caul2)
            caul2 = [caul2;zeros(length(caul11)-length(caul2),1)];
            dstCaul = norm((caul11-caul2).*bOrder1,2);
        elseif length(caul2) > length(caul11)           
            caul11 = [caul11;zeros(length(caul2)-length(caul11),1)];
            dstCaul = norm((caul11-caul2).*bOrder2,2);
        end
        
        
        %concurrence
        maxConc = max([conc1;conc2])+0.0001;
        minConc = min([conc1;conc2]);
        conc11 = (conc1-minConc)/(maxConc-minConc);
        conc2 = (conc2-minConc)/(maxConc-minConc);
        if length(conc11) >= length(conc2)
            conc2 = [conc2;zeros(length(conc11)-length(conc2),1)];
            dstConc = norm((conc11-conc2).*bOrder1,2);
        elseif length(conc2) > length(conc11)           
            conc11 = [conc11;zeros(length(conc2)-length(conc11),1)];
            dstConc = norm((conc11-conc2).*bOrder2,2);
        end
        
        
        %segment length
        maxSegL = mean([segL1;segL2])+0.0001;
        minSegL = std([segL1;segL2]);
        segL11 = (segL1-minSegL)/(maxSegL-minSegL);
        segL2 = (segL2-minSegL)/(maxSegL-minSegL);
        if length(segL11) >= length(segL2)
            segL2 = [segL2;zeros(length(segL11)-length(segL2),1)];
            dstSegL = norm((segL11-segL2).*bOrder1,2)/sqrt(length(segL1));
        elseif length(segL2) > length(segL11)            
            segL11 = [segL11;zeros(length(segL2)-length(segL11),1)];
            dstSegL = norm((segL11-segL2).*bOrder2,2)/sqrt(length(segL2));
        end
                   
        
        %competition
        maxCompt = max([compt1;compt2])+0.0001;
        minCompt = min([compt1;compt2]);
        compt11 = (compt1-minCompt)/(maxCompt-minCompt);
        compt2 = (compt2-minCompt)/(maxCompt-minCompt);
        if length(compt11) >= length(compt2)
            compt2 = [compt2;zeros(length(compt11)-length(compt2),1)];
            dstCompt = norm((compt11-compt2).*bOrder1,2);
        elseif length(compt2) > length(compt11)           
            compt11 = [compt11;zeros(length(compt2)-length(compt11),1)];
            dstCompt = norm((compt11-compt2).*bOrder2,2);
        end
        
        
        %% distance  
        % User can try other non-linear functions 
        dstVec = [dstTort;dstBifur;dstCaul;dstConc;dstSegL;dstCompt];
        dstCell{k,m} = dstVec;        
        distanceMat(k,m) = 1/(1+exp(-[0.0304 0.0372 0.0417 0.1786 0.6489 0.0633]*dstVec));           
    end
end

% distanceMat
[dsTNeuron,prL] = distanceNeuronNew(distanceMat,dataA,dataB);

dstVecforPair = zeros(fDim,size(prL,1));
for nn = 1:size(prL,1)
   dstVecforPair(:,nn) = dstCell{prL(nn,1),prL(nn,2)}; 
end


end