function [spikeArr, spikeNum, img] = PVspikeDetect(dir, filenum, optin);
% %  [spikeArr, spikeNum, img] = PVspikeDetect(dir, filenum)
% %
% %  read ParaVision 2dseq file
% %  detects images affected by spike artefacts 
% %  writes a corrected '2dseq.despiked' file
% %
% %  default options:
% %     opt(  'THRES', 10,               % in Stdev's
% %           'NOISEROI', [1 10 1 10],
% %           'REPLACE', 1,
% %           'WRITE', 0,
% %           'VERBOSE', 1 )
% %
% %  return: 
% %
% %  Functions called: 
% %  Tested for: 
% %
% %  Aug 2003 -  Josef Pfeuffer
% %
FCTNAME = 'PVspikeDetect';
          
SAVENAMEAPP = '.despiked';

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.THRES    = 10;               % in Stdev's
dopt.NOISEROI = [1 10 1 10];
dopt.REPLACE  = 1;
dopt.WRITE    = 0;
dopt.VERBOSE  = 1;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

detectThreshold = dopt.THRES;
N               = dopt.NOISEROI;
f_replace       = dopt.REPLACE;
f_write         = dopt.WRITE;
f_verbose       = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

info = PVrd2dseq(dir,filenum,opt('GETINFO',1,'VERBOSE',0) );
img = PVrd2dseq(dir,filenum,opt('NR',[0 0],'VERBOSE',f_verbose) );
N(2) = min([N(2) info.nx]);
N(4) = min([N(4) info.ny]);
noise = img(N(1):N(2),N(3):N(4),:,:);

s_noise = [size(noise) 1 1]; 
noise = double( reshape(noise, s_noise(1)*s_noise(2),s_noise(3),s_noise(4)) );
noise = squeeze( mean(noise) ); 
s_noise = size(noise); 

NS = s_noise(1);
NR = s_noise(2);

    % determine spike threshold
noisearr = reshape(noise, NS*NR, 1);
medFiltLength = NR;
noisearrFiltered = medfilt1(noisearr, medFiltLength);
meanFiltered = mean(noisearrFiltered);
stdFiltered = std(noisearrFiltered);

detectVal = meanFiltered + detectThreshold*stdFiltered;
spikeArr = zeros(NS,NR);
spikeInd = find(noise > detectVal);
if ~isempty(spikeInd)
    spikeArr(spikeInd) = noise(spikeInd)/meanFiltered;
    spikeNum = length(spikeInd);
else
    spikeNum = 0;
end

if f_verbose
	figure(3)
	plot(reshape(noise,NS*NR,1));
	xlabelStr = sprintf('nSlices(%d)*nRepetitions(%d)  <%s> exp#%d ', NS, NR, dir, filenum);
	xlabel(xlabelStr)
	ylabel('mean intensity')
	titleStr = sprintf('%d  spike images found;  noise = %.2f +/- %.2f;  detectLevel > %.2f (%.1f STD)', ...
        spikeNum, meanFiltered, stdFiltered, detectVal, detectThreshold);
	title(titleStr)
	hold on
	plot([0 NS*NR],[detectVal detectVal],'r-')
	hold off
	
	figure(1)
	imagesc(spikeArr)
	xlabel('nRepetition')
	ylabel('nSlice')
end

     % simple version: replace with image nearby
if f_replace
    for isl=1:NS
        nrInd = find(spikeArr(isl,:) <= 0);
        if ~isempty(nrInd)
         for inr=1:NR
            if spikeArr(isl,inr) > 0
                j = find(nrInd < inr);
                if isempty(j)
                    j = find(nrInd >= inr);
                    if isempty(j)      % should not happen unless there are ONLY spikes
                        error(sprintf('%s: (%d/%d) NO images for replacement found', FCTNAME, isl, inr));
                    else
                        nrInd2 = min(nrInd(j));
                    end
                else
                    nrInd2 = max(nrInd(j));
                end
                img(:,:,isl,inr) = img(:,:,isl,nrInd2);
                fprintf('%s: (%d/%d) image replaced by (%d/%d)\n', FCTNAME, isl, inr, isl, nrInd2);
            end
         end
        end
    end
end

if f_write
    savename = strcat(info.file, SAVENAMEAPP);
    fid = fopen(savename, 'w', info.byteorder);
    if (f_verbose)
        fprintf('writing <%s>\n', savename);
    end
    fcount = fwrite(fid, img, info.precision); 
    fclose(fid);
    flen = info.nx*info.ny*info.nslices*info.nr;
    if (fcount < flen)
        error(sprintf('%s: only %d / %d written', FCTNAME, fcount, flen))
    end
end

%------------------------------------------------------------
