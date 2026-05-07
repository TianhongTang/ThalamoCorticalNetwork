function A=mtspec(TS,NW,K,pad,window,winstep);

% Multitaper Time-Frequency Spectrum
% function A=mtspec(TS,NW,K,pad,window,winstep);
% TS : input time series
% NW = time bandwidth parameter (e.g. 3 or 4)
% K = number of data tapers kept, usually 2*NW -1 (e.g. 5 or 7 for above)
% pad = padding for individual window. Usually, choose power
% of two greater than but closest to size of moving window.
% window = length of moving window
% winstep = number of of timeframes between successive windows


TS=TS(:)';
[E V]=dpss(window,NW,'calc');
[dum N]=size(TS);
A=zeros(round((N-window)/winstep), pad);
for j=1:((N-window)/winstep)
  eJ=zeros(1,pad);
  %TSM=([TS((j-1)*winstep+1: ...
		(j-1)*winstep+window) zeros(1, pad-window)])'*ones(1,K);
  TSM=TS((j-1)*winstep+[1:window])';
  J=fft(TSM(:,ones(1,K)).*(E(:,1:K)),pad)';
  %eJ=sum((abs(J)).^2);
  eJ=sum(J.*conj(J));
  A(j,:)=eJ/K;
end

A=A';
