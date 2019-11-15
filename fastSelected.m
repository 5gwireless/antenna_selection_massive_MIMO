function  capacityOfSelected=fastSelected(Nr,Ns,Lr,SNR,H,fullAntenna);
if(Lr==Nr)
    capacityOfSelected=log2(det(eye(Ns)+SNR/Ns*(H'*H))) ;
else
    
    B=eye(Ns,Ns);
    Alpha=[];
    H_sel=[];%
    
    for j=1:Nr   %
        f=H(j,:);
        h=f';
        alpha=h'*h;
        Alpha=[Alpha alpha]; 
    end
    
    for n=1:Lr    
        [maxOfAlpha,index]=max(Alpha);  
        
        fullAntenna(index)=[]; 
        H_sel=[H_sel;H(index,:)];
        
        if (n<Lr)
            f=H(index,:);
            h=f';
            alpha=Alpha(index);
            a=(B*h)/sqrt((Ns/SNR)+alpha);
            B=B-a*a';                    
            
            Alpha(index)=[];
            
            for k=1:length( fullAntenna)
                Alpha(k)=Alpha(k)-(abs(a'*h))^2;   
            end
            
        end
        
    end
    capacityOfSelected=log2(det(eye(Ns)+SNR/Ns*(H_sel'*H_sel))) ; 
end