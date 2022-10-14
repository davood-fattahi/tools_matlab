function [P, Q, R, S, T] = ecgWavesSoI_RRT(ecg, Rpeak, fs, RRf)
% 'RR interval & time based' version.
% 
% This function appriximately finds peaks and borders of the Q and T waves 
% in the input ecg.   


% Inputs:
% ecg: single channel ecg signal,
% R: R peak positions,
% fs: sampling frequency,
% RRf: R-peak referenced flag; if true, R peak is the reference of all the
% output points. if false, they are based on the samples index.

% Outputs:
% P, Q, R, S, T: Each of the outputs is a matrix, in which the first, second, and third 
% columns respectively contain onset, peak and offset of the P, Q, R, S, T waves. 


% Davood Fattahi, 5/5/2021
% fattahi.d@gmail.com



if length(Rpeak)==length(ecg)
    Rpeak = find(Rpeak);
end

ecg=ecg(:); Rpeak=Rpeak(:);
RR=Rpeak(2:end)-Rpeak(1:end-1); RR = [median(RR); RR; median(RR)];
P=nan(length(Rpeak),3); Q=nan(length(Rpeak),3);
R=nan(length(Rpeak),3); S=nan(length(Rpeak),3);
T=nan(length(Rpeak),3); 
% Uafp=nan(length(Rpeak),3);


for i= 1: length(RR)-1
    P(i,1)=Rpeak(i)-floor(.4*RR(i));
    P(i,3)=Rpeak(i)-floor(.1*RR(i));
    [~,I]=max(abs(ecg(P(i,1):P(i,3)))); P(i,2)=I+P(i,1)-1;
    
    Q(i,1)=Rpeak(i)-floor(.08*fs);    
    Q(i,3)=Rpeak(i)-floor(.02*fs); 
    [~,I]=min(ecg(Q(i,1):Q(i,3))); Q(i,2)=I+Q(i,1)-1;
    
    S(i,3)=Rpeak(i)+floor(.15*(RR(i+1)));
    S(i,1)=Rpeak(i)+floor(.04*(RR(i+1))); 
    [~,I]=min(ecg(S(i,1):S(i,3))); S(i,2)=I+S(i,1)-1;
    
    R(i,1)=Q(i,3);  R(i,3)=S(i,1); R(i,2) = Rpeak(i);

    T(i,3)=Rpeak(i)+floor(.6*RR(i+1));
    T(i,1)=Rpeak(i)+floor(.1*fs);
    [~,I]=max(ecg(T(i,1):T(i,3))); T(i,2)=I+T(i,1)-1;
end


P(P<1|P>length(ecg))=nan;
Q(Q<1|Q>length(ecg))=nan;
R(R<1|R>length(ecg))=nan;
S(S<1|S>length(ecg))=nan;
T(T<1|T>length(ecg))=nan;

% 
% fidp=[Q(:); R(:); S(:); T(:)]; fidp=fidp(~isnan(fidp));
% figure;
% plot(ecg); hold on; plot(fidp, ecg(fidp), 'r*');
% 


if RRf
    P=P-Rpeak; Q=Q-Rpeak; R=R-Rpeak; S=S-Rpeak; T=T-Rpeak;
end

end

