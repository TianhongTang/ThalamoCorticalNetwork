function nkl
    % HDR simulation with matrix setup
%HDRsimulation

   % you will need the code of Zoe to run this:
   % repeatedhistory.m
   % members.m
   % dcgethistory.m
[erDesign,efficency] = erDesignCounterBalanced(3,3,2,[],0);

% % %---------
% % %the following lines are the core of the efficiency calculation
% % 
% % if pluginCorrYes   %%%% THIS IS a MODEL for fMRI noise
% %   b=[0.406;0.8825]; % parameters from SPG fMRI data
% %   fittedACorr=autocorrFnct(b,1:erNtrials+(HDRdur-1)*rem(cutEndPad+1,2));
% %   Cninv=inv(toeplitz(fittedACorr));
% % else
% %   Cninv=eye(erNtrials+(HDRdur-1)*rem(cutEndPad+1,2));
% % end;
% % 
% % designEff = 1/trace(inv(eventMatrix'*Cninv*eventMatrix));

%--------------------------------------------------------------
function [erDesign,efficiency] = erDesignCounterBalanced(numConds, repFactor, ...
    balanceDepth, outFile, flagVerbose)
%
%  
%
%%  [erDesign,efficency] = erDesignCounterBalanced(3,3,2,[],0);
%% ---- 100: efficiency [0.003511 0.010618 83]: 0.008566
%% ---- efficiency: 0.008263 +/- 0.001213

Nsimulations = 10;
timebase = 1; % s
HDRdur = 30/timebase; % 30 sec: assumed # TRs to cover HDR fn

pluginCorrYes=1;
daleYes=1; % 1-usual Dale efficiency. 0-Fisher
% event matrices are defined based on vectors with the average
% energy removed
cutEndPad=0; %cutting reduces efficiency!

for isim=1:Nsimulations

erSeries = repeatedhistory(numConds, repFactor, ...
    balanceDepth, outFile, flagVerbose)';   % Zoe's counterbalanced design
erNtrials = length(erSeries);
erNConditions = numConds;

   % create EventMatrix
mEvent=zeros(erNConditions,erNtrials);
for k=1:erNConditions,
    ind = find(erSeries == k);
    mEvent(k,ind) = 1;
end;
if flagVerbose > 1
    figure(1)
    imagesc(mEvent);
    drawnow;
end

	% defining event matrix as convolution matrix
eventMatrix = makeEventMtrx(mEvent,HDRdur); 
if cutEndPad, eventMatrix=eventMatrix(1:erNtrials,:); end;
eventMatrix = eventMatrix-ones(size(eventMatrix,1),1)* ...
    sum(eventMatrix)/size(eventMatrix,1);
if flagVerbose > 1
    figure(3)
    imagesc(eventMatrix);
end

if pluginCorrYes
  b=[0.406;0.8825]; % parameters from SPG fMRI data
  fittedACorr=autocorrFnct(b,1:erNtrials+(HDRdur-1)*rem(cutEndPad+1,2));
  Cninv=inv(toeplitz(fittedACorr));
else
  Cninv=eye(erNtrials+(HDRdur-1)*rem(cutEndPad+1,2));
end;

if daleYes,
    designEff = 1/trace(inv(eventMatrix'*Cninv*eventMatrix));
else
    designEff = trace(eventMatrix'*Cninv*eventMatrix);
end;


if isim==1
    erDesign = zeros(erNtrials, Nsimulations);
    efficiency = zeros(Nsimulations,1);
end
efficiency(isim) = designEff;
designEffMax = max(efficiency);
designEffMin = min(efficiency(efficiency~=0));
fprintf('---- %d: efficiency [%f %f %d]: %f\n', isim,designEffMin(1),designEffMax(1),length(erSeries),designEff);
erDesign(:,isim) = erSeries';

end  %for: isim=1:Nsimulations

figure(4)
plot(efficiency)
ylim([designEffMin(1) designEffMax(1)]);
fprintf('---- efficiency: %f +/- %f\n', mean(efficiency), std(efficiency) );

return
%--------------------------------------------------------------
function HDRsimulation
%
%  simulation with matrix setup
%

Nsimulations = 10;

timebase = 1; % s
erStim = 2;   % s
erISI = 3;    % s
%%erTrialsEnd = [1 1 1 1 1 1];  % 18 s
HDRdur = 30/timebase; % 30 sec: assumed # TRs to cover HDR fn
erCond = [0 0.7 1];  % index 1=NULL EVENT, 2=50%, 3=100% contrast
f_model = 2;    % 1: Cohen model,  2: SPM model
impulseBroadeningFactor = 1;

pluginCorrYes=1;
daleYes=1; % 1-usual Dale efficiency. 0-Fisher
% event matrices are defined based on vectors with the average
% energy removed
cutEndPad=0; %cutting reduces efficiency!

for isim=1:Nsimulations
   % create paradigm
   % assuming range [1..erNConditions] defining Condition Values
erSeries = repeatedhistory(length(erCond),3,2,[],0)';   % Zoe's counterbalanced design
%erSeries = [1 2 1 3 2 3 2 2 1 2 2 2 3 3 2 1 1 3 3 3 1 3 1 1 1 2 3 1 2];
%erSeries = (erReadExpSeries([0 1 2])+1);
%erSeries = [erSeries erTrialsEnd];
erNtrials = length(erSeries);
erNConditions = length(erCond);

   % create HDR
erHDR = erImpulseFunction(timebase/impulseBroadeningFactor, f_model);
erHDRall = [erHDR(1:HDRdur)*erCond(1) erHDR(1:HDRdur)*erCond(2) erHDR(1:HDRdur)*erCond(3)]';

   % create EventMatrix
mEvent=zeros(erNConditions,erNtrials);
for k=1:erNConditions,
    ind = find(erSeries == k);
    mEvent(k,ind) = 1;
end;
figure(1)
imagesc(mEvent);
drawnow;

erStimulus = zeros(erNConditions,erNtrials*erISI/timebase);
for k=1:erNConditions,
    ind = find(erSeries == k);
    ind = (ind-1)*erISI/timebase+1;
    for iind=1:length(ind)
        erStimulus(k,ind(iind):ind(iind)+erStim/timebase-1) = 1;
    end
end;
%%figure(1)
%%imagesc(erStimulus);

	% defining event matrix as convolution matrix
eventMatrix = makeEventMtrx(mEvent,HDRdur); 
if cutEndPad, eventMatrix=eventMatrix(1:erNtrials,:); end;
eventMatrix = eventMatrix-ones(size(eventMatrix,1),1)* ...
    sum(eventMatrix)/size(eventMatrix,1);
figure(3)
imagesc(eventMatrix);

    % define StimulusConvolutionMatrix (SCM)
figure(4)
erSCM = makeEventMtrx(erStimulus,HDRdur); 
if cutEndPad, erSCM=erSCM(1:erNtrials*erISI/timebase,:); end;
plot(erSCM*erHDRall)

if pluginCorrYes
  b=[0.406;0.8825]; % parameters from SPG fMRI data
  fittedACorr=autocorrFnct(b,1:erNtrials+(HDRdur-1)*rem(cutEndPad+1,2));
  Cninv=inv(toeplitz(fittedACorr));
else
  Cninv=eye(erNtrials+(HDRdur-1)*rem(cutEndPad+1,2));
end;

%%Cninv1=eye(erNtrials*erISI/timebase+(HDRdur-1)*rem(cutEndPad+1,2));
if daleYes,
    designEff = 1/trace(inv(eventMatrix'*Cninv*eventMatrix))
else
    designEff = trace(eventMatrix'*Cninv*eventMatrix)
end;
	
if isim==1
    erDesign = zeros(erNtrials, Nsimulations);
    efficiency = zeros(Nsimulations,1);
end
erDesign(:,isim) = erSeries';
efficiency(isim) = designEff(1);

end  %for: isim=1:Nsimulations

plot(efficiency)

return
%--------------------------------------------------------------
function eventMatrix=makeEventMtrx(events,HDRdur)

% eventMatrix=makeEventMtrx(events,HDRdur)
% creates a length(events) X (HDRdur) convolution matrix 
% for each row of "events" and concatenates them

eventMatrix=[];
nVals=size(events,1);
for i=1:nVals,
   eventMatrix=[eventMatrix convmtx(events(i,:)',HDRdur)];
end;
%--------------------------------------------------------------
function impulseFun = erImpulseFunction(timebase, f_model)

f_verbose = 0;
%%f_model = 1;    % 1: Cohen model,  2: SPM model

deltat = timebase;
time = (0:deltat:30);   % sec

if (f_model == 1)
   % --- MODEL 2: Cohen 1997 -- 
   %     Gamma variate impulse model G(t, gama_a, gam_b, gam_c)
   gam_a = 0.008;
   gam_b = 8.60;
   gam_c = 0.547;
   if ( f_verbose )
      disp( sprintf('Gamma model: %f %f',gam_b, gam_c));
      disp( sprintf('Gamma model: time-to-peak [s]= %f', gam_b*gam_c) );
      disp( sprintf('Gamma model: FWHM [s]~ %f', 2.35 * sqrt(gam_b) * gam_c) );
   end
   impulseFun = funGamma([gam_a gam_b gam_c],time);
   norm = deltat/0.9973;
else
    impulseFun = funGammaSpm(time);
    norm = deltat/1.1036;
end

impulseFun = impulseFun*norm;

figure(3)
plot(time, impulseFun);
if f_verbose key; end

return
%--------------------------------------------------------------
function f = funGamma(a,x)
% f = funGauss(a,x):
% --- Cohen 1997 -- 
%     Gamma variate impulse model G(t, gama_a, gam_b, gam_c)
%   gam_a = 0.008;
%   gam_b = 8.60;
%   gam_c = 0.547;
%
%   JP Apr 2002

if length(a) < 3
   error('funGamma: length(a) < 3 !')
end

f_verbose = 0;

gam_a = a(1);
gam_b = a(2);
gam_c = a(3);

f = gam_a .* (x.^gam_b).* exp(-x/gam_c);
f = f/(sum(f)*(x(2)-x(1)));

if ( f_verbose )
    disp( sprintf('Gamma model: %f %f',gam_b, gam_c));
    disp( sprintf('Gamma model: time-to-peak [s]= %f', gam_b*gam_c) );
    disp( sprintf('Gamma model: FWHM [s]~ %f', 2.35 * sqrt(gam_b) * gam_c) );
    plot(time,f);
end
%--------------------------------------------------------------
function [hrf,p] = funGammaSpm(time,P)
%%% taken from:
%%% function [hrf,p] = spm_hrf(RT,P);
% returns a hemodynamic response function
% FORMAT [hrf,p] = spm_hrf(RT,[p]);
% RT   - scan repeat time
% p    - parameters of the response function (two gamma functions)
%
%							defaults
%							(seconds)
%	p(1) - delay of response (relative to onset)	   6
%	p(2) - delay of undershoot (relative to onset)    16
%	p(3) - dispersion of response			   1
%	p(4) - dispersion of undershoot			   1
%	p(5) - ratio of response to undershoot		   6
%	p(6) - onset (seconds)				   0
%	p(7) - length of kernel (seconds)		  32
%
% hrf  - hemodynamic response function
% p    - parameters of the response function
%_______________________________________________________________________
% @(#)spm_hrf.m	2.7 Karl Friston 99/05/17

% global parameter
%-----------------------------------------------------------------------
global fMRI_T; 
if isempty(fMRI_T), fMRI_T = 16; end;

% default parameters
%-----------------------------------------------------------------------
p     = [6 16 1 1 6 0 32];
if nargin > 1
      p(1:length(P)) = P;
end
p(7)  = time(length(time));

% modelled hemodynamic response function - {mixture of Gammas}
%-----------------------------------------------------------------------
%%dt    = RT/fMRI_T;
dt    = time(2)-time(1);
u     = [0:(p(7)/dt)] - p(6)/dt;
hrf   = spm_Gpdf(u,p(1)/p(3),dt/p(3)) - spm_Gpdf(u,p(2)/p(4),dt/p(4))/p(5);
%%%hrf   = hrf([0:(p(7)/RT)]*fMRI_T + 1);
%%%hrf   = hrf'/sum(hrf);
hrf   = hrf/(sum(hrf)*dt);

%--------------------------------------------------------------
