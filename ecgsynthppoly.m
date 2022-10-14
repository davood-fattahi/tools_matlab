function ECG=ecgsynthppoly(L,HRmean,HRdev,pct,pmean,pdev,fs,noisdev)




%%
    
%%% scaling to [0 2*pi]
beta=range(pct,'all')./(2*pi);
gamma=pct(1);
pcphase=pct./beta+gamma;

for i=1:size(pct,1)
PolyOrder=size(pmean{i,1},2)-1;    
coefs=flip(pmean{i,1});
C=zeros(PolyOrder+1);
for j=0:PolyOrder
    for k=0:j
        C(j+1,k+1)=nchoosek(j,k)*(beta^k)*(gamma^(j-k))*coefs(j+1);
    end
end
coefs=sum(C,1);
P{i,1}=flip(coefs);
end

%%%%

%%
Ts=1./fs;
HR=HRmean+HRdev*randn(1);
ECG=zeros(3,L);
phase=zeros(L,1);
PP=deviate(P,pdev);
for i=2:L
    omega=(HR)*2*pi + noisdev(3).*randn(1);
    phase(i)=rem(phase(i-1)+omega*Ts,2*pi)+noisdev(1)*randn(1);
    if phase(i)< phase(i-1)/2
        HR=HRmean+HRdev*randn(1);
        PP=deviate(P,pdev);
    end
    
    k=(sum(((pcphase-phase(i))>=0),2)==1);
    p=mean(cell2mat(PP(k,1)),1);

    ECG(1,i)=phase(i);
    ECG(2,i)=polyval(p,phase(i))+(noisdev(2).*randn(1));
    ECG(3,i)=omega;
end

end


function PP=deviate(P,pdev)
for l=1:size(P,1)
    ppp=flip(P{l,1});
    Pdev=zeros(size(P{l,1}));
    for m=2%size(P{l,1},2)
        Pdev(m)=(pdev)^(m).*randn(1);
    end 
    PP{l,1}=P{l,1}+flip(Pdev);
end
end










