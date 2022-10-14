function afp =ecgfid_1beat(ecg,fs,hr,varargin)

% This function appriximately finds peaks and borders of the input ecg beat
% waves.   

% afp =Ecgfidapprox(ecg,fs,hr,sm)

% Inputs:
% ecg: one ecg beat,
% fs: sampleing frequency
% hr: approximated heart rate in Hz

% Outputs:
% afp: the vector of approximated fiducal points


% Davood Fattahi, 20/02/2022
% fattahi.d@gmail.com

ecg=ecg(:);
[~,Rp]=max(ecg);

Ps=Rp - floor(.4*(fs/hr)); Ps(Ps<1)=1;
Pe=Rp - floor(.1*(fs/hr));
[~,I]=max(abs(ecg(Ps:Pe))); Pp=I+Ps-1;

if Pp==Pe
    Pp=Pp-2;
elseif Pp==Ps
    Pp=Pp+2;
end

Qs=Pe;
Qe=Rp-floor(.025*(fs/hr)); 
[~,I]=min(ecg(Qs:Qe)); Qp=I+Qs-1;

if Qp==Qe
    Qp=Qp-2;
elseif Qp==Qs
    Qp=Qp+2;
end

Ss=Rp+floor(.04*(fs/hr)); 
Se=Rp+floor(.15*(fs/hr)); 
[~,I]=min(ecg(Ss:Se)); Sp=I+Ss-1;

if Sp==Se
    Sp=Sp-2;
elseif Sp==Ss
    Sp=Sp+2;
end

Rs=Qe;  Re=Ss;


Te=Rp+floor(.6*(fs/hr)); Te(Te>length(ecg))=length(ecg);
Ts=Rp+floor(.15*(fs/hr));
[~,I]=max(ecg(Ts:Te)); Tp=I+Ts-1;

if Tp==Te
    Tp=Tp-2;
elseif Tp==Ts
    Tp=Tp+2;
end

Us=Te;
Ue=size(ecg,1);    Ue(Ue>length(ecg))=length(ecg);
[~,I]=max(ecg(Us:Ue)); Up=I+Us-1;

if Up==Ue
    Up=Up-2;
elseif Up==Us
    Up=Up+2;
end

afp=[Ps Qs Rs Ss Ts Us ; Pp Qp Rp Sp Tp Up; Pe Qe Re Se Te Ue];
afp(afp<1)=1;
afp(afp>length(ecg))=length(ecg);


end

