function  capacityOfSelected=gorokohovSelected(Nr,Ns,Lr,SNR,H,fullAntenna);
if(Lr==Nr)
    capacityOfSelected=log2(det(eye(Ns)+SNR/Ns*(H'*H))) ;
else
    B=inv(eye(Ns,Ns)+SNR/Ns*(H'*H));
    for n=1:(Nr-Lr)  
        Alpha=[];
        for j=1:length(fullAntenna)  
            f=H(j,:);
            h=f';
            alpha=h'*B*h;
            Alpha=[Alpha alpha]; 
        end
        [minOfAlpha,index]=min(Alpha);  
        fullAntenna(index)=[]; 
        
        if (n<Nr-Lr)
            f=H(index,:);
            h=f';
            alpha=Alpha(index);
            a=B*h;
            B=B+a*a'/(Ns/SNR-alpha);
        end
    end
    H_sel=H(fullAntenna,:);
    
    capacityOfSelected=log2(det(eye(Ns)+SNR/Ns*(H_sel'*H_sel))) ; 
end