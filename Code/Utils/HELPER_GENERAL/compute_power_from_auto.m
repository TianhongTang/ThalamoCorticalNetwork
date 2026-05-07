%******** compute power spectrum from autocorrelations ****
function tempinfo = compute_power_from_auto(auto,jackautos,jackrates,poissautos,renewautos,rate,NT,TT,TP)
%***** function info = compute_power_from_auto(auto,jackautos,jackrates,poissautos,renewautos,rate,NT,TT,TP)
%******
%** inputs:
%******** auto - autocorrelation computed (no normalization)
%******** jackautos - several copies, one trial left out
%******** jackrates - rates for each jacknife, must be exact
%******** poissautos - several copies, random permutes of Poisson
%******** renewautos - several copies, randoms of renewal process
%******** rate - mean firing rate, prob spike per unit time
%******** NT - number of trials that went into autocorrelation
%******** TT - total duration of each trial that went into autocor
%******** TP - time sampling of autocorrelations
%******
%*** outputs:
%
%******* tempinfo.plevel - Poisson predicted noise
%******* tempinfo.freqdiv - frequencies of returned power spectrum
%******* tempinfo.autopow - power spectrum
%******* tempinfo.autopow_s - jacknife estimate of power spectrum
%******* tempinfo.rautopow - renewal process spectrum
%******* tempinfo.rautopow_s - error bars from renewal process 
%******* tempinfo.urate - mean heigh of Poissons

   %T = 512;  % 256;   % fft in window of 256 ms, (T/TP) time points
   T = 256;
   
   %***** note, you might think basepow == mean(urate) ....
   %***** but alas, there is a slight bias that must be
   %***** corrected in which power is higher with fewer spikes
   [powauto,freqo,rat,tapers] = normauto_power(auto,T,TP,rate,[]);   
   %********************************************
   
   %********* first compute powers on Poisson samples
   urate = 0;
   pauto = [];
   for kk = 1:size(poissautos,1)
      [tempo,w,rate] = normauto_power(poissautos(kk,:),T,TP,rate,tapers); 
      pauto = [pauto ; tempo];
      urate = urate + rate;
   end
   urate = urate / size(poissautos,1);
   %*************************************************
  
   jauto = [];
   for kk = 1:size(jackautos,1)    
       jrate = jackrates(kk);
       [tempo,w] = normauto_power(jackautos(kk,:),T,TP,jrate,tapers); 
       jauto = [jauto ; tempo];    
   end
   powauto_s = sqrt( (size(jauto,1)-1) * var( jauto ) );   % Jacknife estimate
   %***************************
   
   %*********** error bounds from several renewal process samples
   rauto = [];
   for kk = 1:size(renewautos,1)
      [tempo,w] = normauto_power(renewautos(kk,:),T,TP,rate,tapers);
      rauto = [rauto ; tempo]; 
   end
   rpowauto = mean( rauto );
   rpowauto_s = std( rauto ) / sqrt( size(rauto,1) );    
   
   %************************************************
   tempinfo.pautopow = mean(pauto);
   tempinfo.pautopow_s = std(pauto) / sqrt( size(pauto,1) );
   %****************************
   tempinfo.urate = urate;
   tempinfo.freqdiv = freqo;
   tempinfo.autopow = powauto;
   tempinfo.autopow_s = powauto_s;
   tempinfo.rautopow = rpowauto;
   tempinfo.rautopow_s = rpowauto_s;
   
return;

%******** function to go from my normlaized autocorrelation to the
%******** corresponding power spectrum ****************
function [powo,freq,rate,tapers] = normauto_power(autoboy,N,TP,rate,tapers)
%*** function [powo,freq,rate,tapers] = normauto_power(auto,N,TP,rate,tapers)
%**** input:  autoboy - normalized autocorrelation (1 for poisson)
%****         N - duration in ms for window on power spectrum (suggest 256)
%****         TP - time sampling of autocorrelation, 0.5 ms
%****         rate - mean spike rate, prob of spike per unit time
%****         tapers - if desired, use tapering of autocorrelation
%**** output: powo - 1xF array of power values
%****         freq - 1xF array of frequency values
%****         rate - mean rate in autocorrelation
%****         tapers - return tapers used so they can be reused

%******** find closest power of 2 for FFT 
T = 2^(nextpow2(size(autoboy,2)) - 1);
T = min(floor(N/TP),T);  % use shorter window of N ms
   
%************** compute tapers to use if not already provided
if (isempty(tapers))
  if (1)  
      tap = hanning(2*T);
      tapers = [tap((T+1):(2*T))' tap(1:T)'];
      KW = size(tapers,1);
  else
      tapers = ones(1,(2*T));
      KW = 1;
  end
  %*******normalize for taper magnitude *********
  for kw = 1:KW  
    tapers(kw,:) = tapers(kw,:) / sum( abs( tapers(kw,:) ) );
  end
  %*********************
else
  KW = size(tapers,1);
end


FS = (1000*TP);  
w = FS * (0:(T-1)) / (TP*T);
 
cauto = autoboy(1);  % there would be a peak at delay t=0
                     % but we eliminate it here (set to value of t=1)
cauto = [cauto autoboy(1:T)];
cauto = [cauto autoboy((T-1):-1:1)];
cauto = cauto / rate;
cauto = cauto - 1;
%****************************
dauto = zeros(size(cauto));
dauto(1) = (1/rate)-1 - cauto(1);
                         
%************** this is where you could apply tapers ********
%********** now do the fft and get the amplitude per frequency
tempo = zeros(size(cauto));
for kw = 1:KW
   icauto = cauto .* tapers(kw,:);
   spowo = fft(icauto,size(icauto,2));
   spowo = 4.0 * spowo; %abs(spowo);    %power is (1/2) for fft, and (1/2) autocor averaging
                                % so multiply by 4 to return to A^2 value
   tempo = tempo + spowo;
end
tempo = tempo / KW;
%************************************************************

%******** compute the level due to perfect autocorrelation at orgin
dtempo = zeros(size(dauto));
for kw = 1:KW
   icauto = dauto .* tapers(kw,:);
   spowo = fft(icauto,size(icauto,2));
   spowo = 4.0 * spowo; %abs(spowo);    %power is (1/2) for fft, and (1/2) autocor averaging
                                % so multiply by 4 to return to A^2 value
   dtempo = dtempo + spowo;
end
dtempo = dtempo / KW;
%************************************

%********* fuse back dc offset, and then remove and replace by constant 1
%********* poisson will have exact constant of 1
ftempo = abs( tempo + dtempo );
ftempo = (ftempo - dtempo);
ftempo = ftempo + 1;
%***********************************


%******** the DC value at index 1 is not necessary, eliminate it 
powo = ftempo(2:end);
freq = w(2:end);

return;
