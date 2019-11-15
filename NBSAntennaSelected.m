function  capacityOfNBS=NBSAntennaSelected(Nr,Ns,Lr,SNR,H,fullAntenna)
 capacityOfNBS=0;
 selAntenna=fullAntenna;
 alpha=zeros(Nr,Ns);
 x=0;
 del=[];
 for k=1:Nr
        for l=1:Nr
            if k<=l
                alpha(k,l)=-1;  
            else
                hk=H(k,:);
                hl=H(l,:);
                x=abs(dot(hk,hl));
                alpha(k,l)=x;
            end
        end
    end
    for m=Nr:-1:Lr+1
        [p q]=find(max(max(alpha)));
        Xk=norm(H(p,:));
        Xl=norm(H(q,:));
        if Xk>=Xl
            alpha(q,:)=-1;
            alpha(:,q)=-1;
            del=[del,q];
        else
            alpha(p,:)=-1;
            alpha(:,p)=-1;
            del=[del,p];
        end
    end
    for n=1:length(del)
         x=del(n)
         selAntenna=[selAntenna(1:x-1),selAntenna(x+1:end)];
    end
    H_sel=H(selAntenna,:);
    capacityOfNBS=log2(det(eye(Ns)+SNR/Ns*(H_sel'*H_sel))) ; 

end

