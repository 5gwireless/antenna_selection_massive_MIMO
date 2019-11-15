%  Name:    Saad Mahboob
%  Date :   Dec 2013 
%  rev :    8
%  This code plots conventional MIMO antenna selection algorithms

close all;
clear all;
format long;

tic;
Ns=4; 
Nr=16;
Lr=4;
simulation=1000; %
 
m = 1;
v = 0.8;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));

capacityOfFastAver=[];
capacityOfRandomAver=[];
capacityOfOptimalAver=[];
capacityOfGorokohovAver=[];
capacityOfNBSAver=[];

for SNROfdB=0:20
    SNR= 10^(SNROfdB/10); 
    capacityOfFastSum=0;
    capacityOfRandomSum=0;
    capacityOfOptimalSum=0;
    capacityOfGorokohovSum=0;
    capacityOfNBSSum=0;
    
    antennaSubset=nchoosek([1:Nr],Lr);
    
    for sim=1:simulation
        
        D_shad = lognrnd(mu,sigma,1, Ns);%
%         D_b = diag(D_shad);%%
%         D_b = sqrt(D_b);%
        
        D_b =1
        H=sqrt(1/2)*(randn(Nr,Ns)+1j*randn(Nr,Ns));
        H = H* D_b;
        fullAntenna=[1:Nr];
        
        %% fastAntennaSelected
        capacityOfFastSelected=fastSelected(Nr,Ns,Lr,SNR,H,fullAntenna); 
        capacityOfFastSum=capacityOfFastSum+capacityOfFastSelected;
        
        
        %% randomAntennaSelected µÄ
        capacityOfRandomSelected=randomSelected(Nr,Ns,Lr,SNR,H,fullAntenna); 
        capacityOfRandomSum=capacityOfRandomSum+capacityOfRandomSelected;
        
        
        %% optimalAntennaSelected
        capacityOfOptimalSelected=optimalSelected(Nr,Ns,Lr,SNR,H,antennaSubset);
        capacityOfOptimalSum=capacityOfOptimalSum+capacityOfOptimalSelected;
        
        %% gorokohovAntennaSelected
        capacityOfGorokohovSelected=gorokohovSelected(Nr,Ns,Lr,SNR,H,fullAntenna);
        capacityOfGorokohovSum=capacityOfGorokohovSum+capacityOfGorokohovSelected;%
        
        %%NBSAntennnaSelected
        capacityOfNBSSelected=NBSAntennaSelected(Nr,Ns,Lr,SNR,H,fullAntenna);
        capacityOfNBSSum=capacityOfNBSSum+capacityOfNBSSelected;%
        
        
    end
    capacityOfFastAver=[capacityOfFastAver,capacityOfFastSum/simulation];
    capacityOfRandomAver=[capacityOfRandomAver,capacityOfRandomSum/simulation];
    capacityOfOptimalAver=[capacityOfOptimalAver,capacityOfOptimalSum/simulation];
    capacityOfGorokohovAver=[capacityOfGorokohovAver,capacityOfGorokohovSum/simulation];
    capacityOfNBSAver=[capacityOfNBSAver,capacityOfNBSSum/simulation];
    
    
end
X=[0:20];

%plot(X,capacityOfOptimalAver,'c',X,capacityOfFastAver,'k-o',X,capacityOfGorokohovAver,'b-+',X,capacityOfRandomAver,'r',X,capacityOfNBSAver,'b','LineWidth',4);
plot(X,capacityOfOptimalAver,'k.-',X,capacityOfFastAver,'k-o',X,capacityOfGorokohovAver,'k-+',X,capacityOfRandomAver,'k-.',X,capacityOfNBSAver,'k--','LineWidth',1);

legend('Optimal selected','Fast selected','Gorokohv selected','Random selected','NBS selected');
xlabel('SNR (dB)');
ylabel('Ergodic Capacity (bps/Hz)');
grid on;
toc;