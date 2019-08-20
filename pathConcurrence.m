function concurList = pathConcurrence(neuron1,concurList,pvt)

if pvt == 1
    return
else    
    pvt = neuron1(pvt,7);        
    concurList(pvt) = concurList(pvt)+1;  
    concurList = pathConcurrence(neuron1,concurList,pvt);
end    

end

