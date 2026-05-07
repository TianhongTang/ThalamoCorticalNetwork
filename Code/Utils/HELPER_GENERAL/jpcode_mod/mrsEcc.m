function mrsDat = mrsEcc(mrsDat, filenum, optin);
% %  mrsDat = mrsEcc(mrsDat, filenum, optin);
% %
% %  perform Eddy Current Corretion on MRS data
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
FCTNAME = 'mrsEcc';

FIGUREID = 4;
ECCREFNUM = 0;   % define #refFID, if 0: average over all FIDs
          
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.THRES   = 0.001;    %% amplitude related threshold for phase of reference FID
dopt.PLOT    = 1;
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

ampThreshold = dopt.THRES;
f_plot    = dopt.PLOT;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

eccDat = mrsRead(mrsDat.dir, filenum);

refnum = min([ECCREFNUM eccDat.info.nr]);
if refnum > 0 then
    eccFid = eccDat.tdat(:,refnum);
else
    eccFid = avg(eccDat.tdat,2);
end
eccFidTime  = (1:length(eccFid))'/eccDat.acqp.SW_h*1000;  % [ms]
eccFidAbs   = abs(eccFid);
eccFidPhase = angle(eccFid);
     %% perform thresholding of phase based on amplitude
[eccFidAbsMax maxInd] = max(eccFidAbs);
indThreshold = find(eccFidAbs < ampThreshold*eccFidAbsMax);
eccFidPhase(indThreshold) = 0;
eccFidPhase = unwrap(eccFidPhase);
%%if maxInd > 1
%%    eccFidPhase(1:maxInd-1) = 0;
%%end
    %% make sure NOT to phase correct anything beyond cutoff 
cutoffInd = indThreshold( find(indThreshold > maxInd) );
cutoffInd = cutoffInd(1);
eccFidPhase(cutoffInd:end) = 0;

if f_plot
	figure(FIGUREID)
	subplot(2,1,1)
	plot(eccFidTime,eccFidAbs)
	ylabel('abs(fid)')
	% hold on
	% for inr = 2:s_tdat(2)
	%     plot(eccFidTime,abs(eccDat.tdat(:,inr)))
	% end
	% hold off
	title(sprintf('ECC fid - cutoffIndex = %d (%.1f ms)', cutoffInd, eccFidTime(cutoffInd) ));
	
	subplot(2,1,2)
	plot(eccFidTime, eccFidPhase)
	xlabel('acq time / ms')
	ylabel('angle(fid)')
	% hold on
	% for inr = 2:s_tdat(2)
	%     plot(eccFidTime,unwrap(angle(eccDat.tdat(:,inr))))
	% end
	% hold off
    pause(0.1)
end

if mrsDat.flagEccDone > 0 
    return
end

   %%% perform ECC correction
if f_verbose
    fprintf('--- %s: ampThreshold = %f\n', FCTNAME, ampThreshold)
end
s_tdat = size(mrsDat.tdat);
td = s_tdat(1);
nr = s_tdat(2);
tdatEcc  = complex( zeros(td,nr) );
eccDat = complex( zeros(td,1) );
if length(eccFidPhase) < td
    eccDat(1:length(eccFidPhase)) = eccFidPhase;
else
    eccDat(:) = eccFidPhase(1:td);
end
for inr = 1:nr      %% 2D loop
    tdatEcc(:,inr) = mrsDat.tdat(:,inr) .* exp(-i*eccDat);
end

mrsDat.tdat   = tdatEcc;
mrsDat.eccDat = eccDat;   %% store ECC phase
mrsDat.flagEccDone = 1;

      
%------------------------------------------------------------
%------------------------------------------------------------
function  ecPlot(filenumarr)

P = [1 350];
ec1 = MRSecPlot('jptest.gI1',filenumarr(1));
ec2 = MRSecPlot('jptest.gI1',filenumarr(2));
ec3 = MRSecPlot('jptest.gI1',filenumarr(3));

figure(1)
plot( unwrap(angle( ec1.tdat(P(1):P(2)) )))
hold on
plot( unwrap(angle( ec2.tdat(P(1):P(2)) )))
plot( unwrap(angle( ec3.tdat(P(1):P(2)) )))
hold off

%------------------------------------------------------------
function  mrsDat = MRSecPlot(dir, filenum)

j = MRSrd(dir, filenum);

s_tdat = size(j.tdat);
time = (1:s_tdat(1))/j.acqp.SW_h*1000;  % [ms]

figure(3)
subplot(2,1,1)
plot(time,abs(j.tdat(:,1)))
ylabel('abs(fid)')
hold on
for inr = 2:s_tdat(2)
    plot(time,abs(j.tdat(:,inr)))
end
hold off
title('time domain data')

subplot(2,1,2)
plot(time,unwrap(angle(j.tdat(:,1))))
xlabel('acq time / ms')
ylabel('angle(fid)')
hold on
for inr = 2:s_tdat(2)
    plot(time,unwrap(angle(j.tdat(:,inr))))
end
hold off

mrsDat = j;
