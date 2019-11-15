function  capacityOfSubsetMax=optimalSelected(Nr,Ns,Lr,SNR,H,antennaSubset)

    
    capacityOfSubsetMax=0;
    
    for k=1:nchoosek(Nr,Lr) 
        indexOfChannel=antennaSubset(k,:);
        H_sel=H(indexOfChannel,:);
        capacityOfSubset=log2(det(eye(Ns)+SNR/Ns*(H_sel'*H_sel))) ; 
        
        if(capacityOfSubset>capacityOfSubsetMax)
            capacityOfSubsetMax=capacityOfSubset;
        end
    end
