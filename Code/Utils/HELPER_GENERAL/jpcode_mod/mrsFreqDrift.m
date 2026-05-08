function freqArr = mrsFreqDrift(mrsDat, optin);
% %  freqArr = mrsFreqDrift(mrsDat, optin)
% %
% %  determine frequency drift by fitting phase changes in time domain
% %
% %  default options:
% %     opt(  
% %           'VERBOSE','1')
% %
% %  return: 
% %
% %  Functions called: 
% %  Tested for: 
% %
% %  Nov 2002 -  Josef Pfeuffer
% %
FCTNAME = 'mrsFreqDrift';

TR = mrsDat.acqp.PVM_RepetitionTime/1000/60;   % [min]
          
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.CUTOFF  = 0.10;
dopt.REFNUM  = 5;
dopt.PLOT    = 0;
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 1;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

if dopt.CUTOFF >0 & dopt.CUTOFF <= 1
    levelCutoff = dopt.CUTOFF;         %% relative to Maximum of Amplitude
    CutoffIndFix = 0;                  %% 0: DISABLE; manual setting of cutoff, 
else if dopt.CUTOFF > 1
        levelCutoff = 1;               %% relative to Maximum of Amplitude
        CutoffIndFix = dopt.CUTOFF;    %% 0: DISABLE; manual setting of cutoff, 
    else
        CutoffIndFix = 0;                  %% 0: DISABLE; manual setting of cutoff, 
    end
end
refnum = dopt.REFNUM;
f_plot = dopt.PLOT;                 %% [s] 0: NO PLOT; plot each linreg of data: pause(time)- 
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

   %% get info about target data
s_tdat = size(mrsDat.tdat);
nr = s_tdat(2);
refFidNum = min([refnum nr]);
freqArr = zeros(nr,2);

% refPhase = angle(mrsDat.tdat(:,refnum));
% for inr = 1:nr
%     deltaPhase = angle(mrsDat.tdat(:,inr)) - refPhase;
%     
% end

dat = mrsDat.tdat;
refFid = double(dat(:,refFidNum));
refFidAbs = abs(refFid);
refFidPhase = angle(refFid);   
refFidMax = max(refFidAbs);
CutoffInd = min( find(refFidAbs < (levelCutoff*refFidMax)) );
if isempty(CutoffInd) | (CutoffInd <=1)
    CutoffInd = length(refFid);
end
if CutoffIndFix > 0     % manual setting
    CutoffInd = min([CutoffIndFix length(refFid)]);
end
fprintf('CutoffInd = %d / %d\n', CutoffInd, length(refFid));
refFidCutoffInd = 1:CutoffInd;

datDeltaPhase = unwrap(angle(dat)); 
datDeltaPhase = datDeltaPhase(refFidCutoffInd,:);
refFidPhase = refFidPhase(refFidCutoffInd);

if (f_plot > 0)
    figure(4)
end
X = (1:length(datDeltaPhase(:,1)))'/(mrsDat.acqp.SW_h);
for i=1:nr
    Y = unwrap(refFidPhase - datDeltaPhase(:,i));
    % [ A0, A1, Yout ] = wlinearfit(X,Y,1);
    pcoeff= polyfit(X,Y,1);
    A1 = pcoeff(1);    % freq 
    A0 = pcoeff(2);  
    Yout = polyval([ A1, A0 ], X);
    freqArr(i,:) = [A0, A1/(2*pi)];    % phase in [rad] to Frequency shift
    if (f_plot > 0)
        fprintf('%d ',i)
        plot(X,Y)
        hold on
        plot(X,Yout)
        xlabel('acq time / ms')
        hold off
        pause(f_plot)
        plot(unwrap(refFidPhase - datDeltaPhase(:,i)))
    end 
end

if f_verbose
	figure(1)
	subplot(2,1,1)
	FX = (1:length(freqArr(:,1)))*TR;
	plot(FX,freqArr(:,1))
	ylabel('D phi0 [rad]')
	title(sprintf('Phase intercept (t=0)'))
	
	subplot(2,1,2)
	plot(FX,freqArr(:,2))
	xlabel('time [min]')
	ylabel('D frequency [Hz]')
	title(sprintf('Frequency fluctuation'))
end


%------------------------------------------------------------
