function [ecgB, rp, rps]=ecgBeats(ecg,fs,varargin)

% Syntax:
% [ecgB, rp, rps]=ecgmean(ecg,fs)
% [ecgB, rp, rps]=ecgmean(ecg,fs,peaks)
% [ecgB, rp, rps]=ecgmean(ecg,fs,peaks,rpr)
% [ecgB, rp, rps]=ecgmean(ecg,fs,[],rpr)


% Description:
% beat-wise ECG segmenting, from the single
% channel ecg signal (common morphology of lead II). First, if R-peaks are 
% not provided in inputs, they are detected using Pan-Tompkin method.
% Then the beats are segmented by aligning the R peaks, while each segment 
% contains 100*rpr percent of the previous RR interval and 100*(1-rpr) percent
% of the next one. No time warping is applied, and the segments are equalized
% in size using zeropadding.
%
%
% Inputs:
% ecg: 1-d vector of ecg record
% fs: sampling frequency
% peaks: vector of R-peaks positiopn, with tha same size of ecg, having
% true values at R-peak positions and false at other samples.
% rpr: R-peak position ratio in the segmented beats. 100*rpr percent 
% of the beats falls before R-peak and 100*(1-rpr) percent
% of them falls after R-peak. The default value is 0.4.
%
%
% Output:
% ecgB:ECG beats
% rp = R-peak position in the outputs.
% rps = R-peaks in the input ecg.
% 
% 
% Davood Fattahi, 12/08/2022
% update 14/10/2022: comment correction.
% fattahi.d@gmail.com



if nargin > 5 || (nargin==5) || (nargin>=4 && ~isscalar(varargin{2})) ||  ~isscalar(fs)
    error('wrong inputs');
end

if nargin==2 || (nargin>2 && isempty(varargin{1}))
    [~,peakspos,~]=pan_tompkin(ecg,fs,0);
    rps = false(size(ecg)); rps(peakspos)=true;
else
    rps = varargin{1};
    if length(rps)~=length(ecg)
        error('wrong size of peaks!');
    end
     peakspos = find(rps);
end

if nargin>=4 && isscalar(varargin{2})
    rpr=varargin{2};
else
    rpr=.4;
end



RRint=peakspos(2:end)-peakspos(1:end-1);
refRRint = ceil(prctile(RRint,90));
rp = ceil(rpr*refRRint);
RRint=[RRint(1) RRint RRint(end)];
ecgB=zeros(size(peakspos,2),refRRint);
for i=1:size(peakspos,2)
    ss=peakspos(i)-floor(rpr*RRint(i));
    if ss<1
        ss=1;
    end
    ee=peakspos(i)+floor((1-rpr)*RRint(i+1));
    if ee>size(ecg(:),1)
        ee=size(ecg(:),1);
    end
    s=rp - peakspos(i)+ss;
    e=rp - peakspos(i)+ee;
    if s<1
        ss=ss-s+1;
        s=1;
    end
    if e>refRRint
        ee=ee+refRRint-e;
        e=refRRint;
    end
    ecgB(i,s:e)=ecg(ss:ee);
end






