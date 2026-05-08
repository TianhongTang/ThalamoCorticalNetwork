function mrsDat = mrsConvdta(mrsDat, optin);
% %  mrsDat = mrsConvdta(mrsDat, optin)
% %
% %  convert BRUKER Avance data from digital to analog data format
% %
% %      sets mrsDat.flagConvdtaDone = 1(successful), 2(not done: analog data);
% %
% %% adapted from CJ: datConvdta = CSI1_Convdta(datRaw)
% %% converts digitally filtered Avance aquistion data to 'analog' data format
% %% full (single spectrum) version with additional tools can be found at D:/juchem/MATLAB/LCModel/convdta_test.m
% %% http://garbanzo.scripps.edu/nmrgrp/wisdom/dig.nmrfam.txt
% %% http://www.acdlabs.com/publish/offlproc.html
% %
% %  default options:
% %     opt(  
% %           'VERBOSE','1')
% %
% %  return: mrsDat.flagConvdtaDone = 1(successful), 2(not done: analog data);
% %
% %  Functions called: AvanceConvList()
% %  Tested for: 
% %
% %  Nov 2002 -  Josef Pfeuffer
% %
FCTNAME = 'mrsConvdta';

FUDGEZEROPHASE  = 90;            % 224: determined by juchem/MATLAB/general/SaFOrderDetermination.m
f_Debug         = 0;
f_zeroEnd       = 0.05;          % cut end of FID
          
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.SHIFT   = 0;    % manual shift
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 1;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

shiftManual = dopt.SHIFT;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%
      
if mrsDat.flagConvdtaDone > 0 
    return
end
if ~strcmp(mrsDat.acqp.DIGMOD, 'digital_mode')
    fprintf('%s: NO convdta performed - DIGMOD<%s>\n', FCTNAME, mrsDat.acqp.DIGMOD);
    mrsDat.flagConvdtaDone = 2;
    return
end

    %----- determination of PhaseFac & datShift
if shiftManual <= 0 
    ConvDat = AvanceConvList(mrsDat.acqp.DECIM, mrsDat.acqp.DSPFVS);
    if f_verbose
        fprintf('--- %s: shift = %.2f\n', FCTNAME, ConvDat)
    end
else
    ConvDat = shiftManual;
    fprintf('--- %s: manual shift = %.2f\n',FCTNAME,shiftManual);
end
datShift = fix(ConvDat);
PhaseFac = mod(ConvDat,1);
FirstPhase  = PhaseFac*360;      % 48/12 -> 190deg (techn. report), or juchem/MATLAB/general/SaFOrderDetermination.m  
                                 % 64/12 -> 139deg        % Bruker notation for phasing is used
s_tdat = size(mrsDat.tdat);
td = s_tdat(1) - datShift;
nr = s_tdat(2);
tdatConvdat = complex( zeros(td,nr) );
if FirstPhase ~= 0
    FirstPhaseArr = (0:FirstPhase/(td-1):FirstPhase)'; 
end

for inr = 1:nr      %% 2D loop
    
    datRaw = mrsDat.tdat(:,inr);
    
    %------- circular shift/data rotation --------------   %% left shift: cut off
    dat = datRaw(datShift+1:end);
    %------- zero/1st order phase correction -------
    if FirstPhase > 0 | FUDGEZEROPHASE ~= 0
        datFt = fftshift(fft(dat,[],1),1); 
        if FirstPhase > 0
            datFt = datFt .* exp(i*pi/180*(FUDGEZEROPHASE + FirstPhaseArr));
        else
            datFt = datFt .* exp(i*pi/180*(FUDGEZEROPHASE));
        end
        dat = ifft(ifftshift(datFt,1),[],1);
    end
    
    tdatConvdat(:,inr) = dat;

    if(f_Debug)
        c_ChFig = 0;
        figh = figure(5+c_ChFig);
        set(figh,'Name',FCTNAME);
        subplot(2,2,1)
        plot(real(fftshift(fft(datRaw))));
        % axis([1 TD2 -10e5 10e5])
        title('real(fft(datRaw))')
        subplot(2,2,2)
        plot(imag(fftshift(fft(datRaw))));
        % axis([1 TD2 -10e5 10e5])
        title('imag(fft(datRaw))')
        subplot(2,2,3)
        plot(real(fftshift(fft(tdatConvdat(:,inr)))));
        % axis([1 TD2 -10e5 10e5])
        title('real(fft(datConvdta))')
        subplot(2,2,4)
        plot(imag(fftshift(fft(tdatConvdat(:,inr)))));
        % axis([1 TD2 -10e5 10e5])
        title('imag(fft(datConvdta))')     
        key(inr)
    end
end

if f_zeroEnd > 0 
    s_tdat = size(tdatConvdat);
    tdatConvdat = tdatConvdat(1:fix(s_tdat(1)*(1-f_zeroEnd)),:);
end

mrsDat.tdat = tdatConvdat;
mrsDat.flagConvdtaDone = 1;

return

%------------------------------------------------------------
function ConvDat = AvanceConvList(Decim, Pspfvs);
%
%
%% function that uses DECIM and DSPFVS (acqus file, part of global acqp parameter)
%% to give back the corresponding value for data rotation and 1st order phasing
%% -> technical report of Westler & Abildgaard, 1996
%% http://garbanzo.scripps.edu/nmrgrp/wisdom/dig.nmrfam.txt
%% http://www.acdlabs.com/publish/offlproc.html
%% 10-2002, Christoph Juchem
%
%

Dspfvs10Vec = [44.75 33.5 66.625 59.0833 68.5625 60.375 69.5313 61.0208 70.0156 61.3438 70.2578 61.5052 70.3789...
               61.5859 70.4395 61.6263 70.4697 61.6465 70.4849 61.6566 70.4924];
Dspfvs11Vec = [46 36.5 48 50.1667 53.25 69.5 72.25 70.1667 72.75 70.5 73 70.6667 72.5 71.3333 72.25 71.6667 72.125 71.8333 72.0625 71.9167 72.0313];
Dspfvs12Vec = [46.311 36.53 47.87 50.229 53.289 69.551 71.6 70.184 72.138 70.528 72.348 70.7 72.524 0 0 0 0 0 0 0 0];

DatMat = zeros(3,21);
DatMat = [Dspfvs10Vec(1,:)',Dspfvs12Vec(1,:)',Dspfvs12Vec(1,:)'];

switch Decim
    case 2, DecimInd=1;
    case 3, DecimInd=2;
    case 4, DecimInd=3;
    case 6, DecimInd=4;
    case 8, DecimInd=5;
    case 12, DecimInd=6;
    case 16, DecimInd=7;
    case 24, DecimInd=8;
    case 32, DecimInd=9;
    case 48, DecimInd=10;
    case 64, DecimInd=11;
    case 96, DecimInd=12;
    case 128, DecimInd=13;
    case 192, DecimInd=14;
    case 256, DecimInd=15;
    case 384, DecimInd=16;
    case 512, DecimInd=17;
    case 768, DecimInd=18;
    case 1024, DecimInd=19;
    case 1536, DecimInd=20;
    case 2048, DecimInd=21;
    otherwise, error(sprintf('DECIM does not have a legal value\n'));
end

switch Pspfvs
    case 10, PspfvsInd=1;
    case 11, PspfvsInd=2;
    case 12, PspfvsInd=3;
    otherwise, error(sprintf('PSPFVS does not have a legal value\n'));
end

if(PspfvsInd==3 & DecimInd>13)
    error(sprintf('combination of DECIM and PSPFVS is not possible\n'));
end

ConvDat = DatMat(DecimInd,PspfvsInd);
%------------------------------------------------------------
