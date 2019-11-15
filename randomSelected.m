function  capacityOfSelected=randomSelected(Nr,Ns,Lr,SNR,H,fullAntenna);
if(Lr==Nr)
    capacityOfSelected=log2(det(eye(Ns)+SNR/Ns*(H'*H))) ;
else
    H_sel=[];
    for n=1:Lr   
        randomIndex=randint(1,1,[1,length( fullAntenna)]);
        H_sel=[H_sel;H(randomIndex,:)];
        fullAntenna(randomIndex)=[]; 
    end
    
    capacityOfSelected=log2(det(eye(Ns)+SNR/Ns*(H_sel'*H_sel))) ; 
end