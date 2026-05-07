function SesSystemTest(SESSION, optin);
% %  SesSystemTest(SESSION);
% %  get imaging parameters directly from the acqp/reco
% %  create data array/files
% %
% %  default options:
% %     opt(  
% %           'VERBOSE','1')
% %
% %  return: global t1 t2 t2s
% %
% %  Functions called: 
% %  Tested for: 
% %
% %  Apr 2003 -  Josef Pfeuffer
% %
FCTNAME = 'SesSystemTest';

global STDPATH acqp reco
global noise noisefig spike

savename = 'SesSystemTest.mat';

%%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.READMAT = 0;
dopt.WRITE   = 0;
dopt.VERBOSE = 0;

      % --- arg handling
nargVars = 1;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

f_readmat = dopt.READMAT;
f_write   = dopt.WRITE;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

Ses = goto(SESSION);
%%infogen = getpvpars(SESSION,'t1',1)

if f_readmat
    fprintf('\n<%s> loading\n', savename)
    load(savename);
    f_t1 = 0;
    f_t2 = 0;
    f_t2s = 0;
    f_write = 0;
else
    noiseind = Ses.grp.noise.exps;
    if noiseind(1) <= 0
        f_noise = 0;
    else
        f_noise = 1;
    end
    
    noisefigind = Ses.grp.noisefig.exps;
    if noisefigind(1) <= 0
        f_noisefig = 0;
    else
        f_noisefig = 1;
    end

    spikeind = Ses.grp.spike.exps;
    if spikeind(1) <= 0
        f_spike = 0;
    else
        f_spike = 1;
    end
%     t1 = [];       %% loading new study: clear all data
%     t2 = [];
%     t2s = [];
%     t1map = [];    %% loading new study: clear all maps
%     t2map = [];
%     t2smap = [];
end

    %%%%% read noise data
if (f_noise == 1)
	for ind=1:length(noiseind)
		ExpNo = noiseind(ind);
		dataPath.stdpath	= Ses.DataMri;
		dataPath.dir		= Ses.dirname;
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);
		STDPATH.pv = dataPath.stdpath;
		
		if (ind == 1)
            info = PVrdFid(dataPath.dir, dataPath.filenum, opt('GETINFO',1,'VERBOSE',f_verbose) );
            noise.dat = zeros(info.nx/2, info.ny, length(noiseind), info.nslices, info.nr);
            noise.stat = zeros(3, length(noiseind));  % in percent
            noise.Ses = Ses;
            noise.acqp = acqp;
            noise.reco = reco;
            normfac = 2^16;
        end
        tdat = PVrdFid(dataPath.dir, dataPath.filenum, opt('VERBOSE',f_verbose) );
		noise.dat(:,:,ind,:,:) = real(double(tdat(:,1)));
        noise.stat(1,ind) = median(real(double(tdat(:,1)))) / normfac * 100;
        noise.stat(2,ind) = std(real(double(tdat(:,1)))) / normfac * 100;
        noise.stat(3,ind) = (max(real(double(tdat(:,1)))) - min(real(double(tdat(:,1))))) / (2*normfac) * 100;
        
        figure(3)
        plot(real(double(tdat(:,1))),'k')
        ylim([-normfac normfac])
        title( sprintf('median = %f %%,  std = %f %%,  p-p = %f %%', ...
             noise.stat(1,ind), noise.stat(2,ind), noise.stat(3,ind)) )
        xlabel( sprintf('<%s/%d> ExpNo<%d>', dataPath.dir, dataPath.filenum, ExpNo ))
        fprintf('%s %s, RG=%d, SP=<%.1f,%.1f>, <%.4f,%.2f,%.2f>%%\n', acqp.PULPROG, acqp.DIGMOD, acqp.RG, acqp.SP(1), acqp.SP(2),...
            noise.stat(1,ind), noise.stat(2,ind), noise.stat(3,ind))
        key
	end
end  % f_t1

    %%%%% read noise figure data
if (f_noisefig == 1)
	for ind=1:length(noisefigind)
		ExpNo = noisefigind(ind);
		dataPath.stdpath	= Ses.DataMri;
		dataPath.dir		= Ses.dirname;
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);
		STDPATH.pv = dataPath.stdpath;
		
		if (ind == 1)
            info = PVrdFid(dataPath.dir, dataPath.filenum, opt('GETINFO',1,'VERBOSE',f_verbose) );
            %%noisefig.dat = zeros(info.nx/2*info.ny, 1, length(noisefigind), info.nslices, info.nr);
            noisefig.dat = zeros(info.nx/2*info.ny, 1, length(noisefigind), 1, 1);
            noisefig.stat = zeros(4, length(noisefigind));  % in percent
            noisefig.NF  = [0 0];
            noisefig.Ses = Ses;
            noisefig.acqp = acqp;
            noisefig.reco = reco;
            normfac = 2^16;       % PV scales 32 bit data to [-2^16 .. +2^16]
        end
        tdat = PVrdFid(dataPath.dir, dataPath.filenum, opt('VERBOSE',f_verbose) );
        s_tdat = size(tdat);
        tdat = reshape(tdat(:,:,1,1), s_tdat(1)*s_tdat(2), 1);
		noisefig.dat(:,:,ind,:,:) = double(tdat(:,:,1,1));
        noisefig.stat(1,ind) = median(real(double(tdat(:,1)))) / normfac * 100;
        noisefig.stat(2,ind) = std(real(double(tdat(:,1)))) / normfac * 100;
        noisefig.stat(3,ind) = (max(real(double(tdat(:,1)))) - min(real(double(tdat(:,1))))) / (2*normfac) * 100;
        noisefig.stat(4,ind) = sqrt( sum(abs(double(tdat(:,1))).^2) / length(tdat(:,1))) / normfac * 100;  %RMS
        
        figure(3)
        plot(real(double(tdat(:,1))),'k')
        ylim([-normfac normfac])
        title( sprintf('median = %f %%,  std = %f %%,  p-p = %f %%, RMS = %.2f %%', ...
             noisefig.stat(1,ind), noisefig.stat(2,ind), noisefig.stat(3,ind), noisefig.stat(4,ind)) )
        xlabel( sprintf('<%s/%d> ExpNo<%d>', dataPath.dir, dataPath.filenum, ExpNo ))
        fprintf('%s %s, RG=%d, SP=<%.1f,%.1f>, <%.4f,%.2f,%.2f>%% <%.2f>%%\n', acqp.PULPROG, acqp.DIGMOD, acqp.RG, acqp.SP(1), acqp.SP(2),...
            noisefig.stat(1,ind), noisefig.stat(2,ind), noisefig.stat(3,ind), noisefig.stat(4,ind))
        if ind < length(noisefigind)
            key
        end
	end
    noisefig1 = CalculateNoiseFigure(noisefig.dat(:,:,1), noisefig.dat(:,:,2), Ses.grp.noisefig.attCal);
end  % f_t1

    %%%%% read noise data
if (f_spike == 1)
	for ind=1:length(spikeind)
		ExpNo = spikeind(ind);
		dataPath.stdpath	= Ses.DataMri;
		dataPath.dir		= Ses.dirname;
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);
		STDPATH.pv = dataPath.stdpath;
		
		if (ind == 1)
            info = PVrdFid(dataPath.dir, dataPath.filenum, opt('GETINFO',1,'VERBOSE',f_verbose) );
            spike.dat = zeros(info.nx/2*info.ny*info.nslices*info.nr, length(spikeind) );
            spike.stat = zeros(3, length(spikeind));  % in percent
            spike.Ses = Ses;
            spike.acqp = acqp;
            spike.reco = reco;
            normfac = (2^32/2-1);
        end
        tdat = PVrdFid(dataPath.dir, dataPath.filenum, opt('VERBOSE',f_verbose) );
        tdat = reshape(tdat, prod(size(tdat)), 1);
        spike.dat(:,ind) = real(double(tdat(:,1)));
        spike.stat(1,ind) = median(real(double(tdat(:,1)))) / normfac * 100;
        spike.stat(2,ind) = std(real(double(tdat(:,1)))) / normfac * 100;
        spike.stat(3,ind) = (max(real(double(tdat(:,1)))) - min(real(double(tdat(:,1))))) / 2^32 * 100;
        plot(real(double(tdat(:,1))),'k')
        ylim([-2^32/2 2^32/2-1])
        title( sprintf('median = %f %%,  std = %f %%,  p-p = %f %%', ...
             spike.stat(1,ind), spike.stat(2,ind), spike.stat(3,ind)) )
        xlabel( sprintf('<%s/%d> ExpNo<%d>', dataPath.dir, dataPath.filenum, ExpNo ))
        fprintf('%s %s, RG=%d, SP=<%.1f,%.1f>, <%.4f,%.2f,%.2f>%%\n', acqp.PULPROG, acqp.DIGMOD, acqp.RG, acqp.SP(1), acqp.SP(2),...
            spike.stat(1,ind), spike.stat(2,ind), spike.stat(3,ind))
        key
	end
end  % f_t1

if f_write
    key('SAVE ? ')
    save(savename,'t1','t2','t2s','t1map','t2map','t2smap');
    fprintf('\n<%s> saved\n', savename)
    
    %load(savename);
    %%STIMwrSdt( t1map,sprintf('%s%s/%d_t1map', t1.Ses.DataMri,t1.Ses.dirname, t1.Ses.expp(t1.Ses.grp.t1.exps(1)).scanreco(1)), opt('ROTATE',0));

    if ~isempty(t1)
        s_dat = size(t1.dat);
        for inr=1:s_dat(3)
            t1dat(:,:,:,inr) = t1.dat(:,:,inr,:);
        end
        STIMwrSdt( t1dat, sprintf('%s%s/%d_t1dat', t1.Ses.DataMri,t1.Ses.dirname,t1.Ses.expp(t1.Ses.grp.t1.exps(1)).scanreco(1)), opt('ROTATE',0));
    end
end


%------------------------------------------------------------
function noisefig = CalculateNoiseFigure(fid1, fid2, excessNoise)

f_verbose = 1;
f_temperature = 0;

fid1 = fid1(:,1,1,1);  % take only 1D data
fid2 = fid2(:,1,1,1);

    %/*------------------------------------------------------------*
    %* - Separate FIDs in memory into real & imaginary components
    %*------------------------------------------------------------*/
fid1r = real(fid1);
fid1i = imag(fid1);
fid2r = real(fid2);
fid2i = imag(fid2);

    %/*------------------------------------------------------------*
    % * - Remove DC offsets from each channel for both FIDs
    % *------------------------------------------------------------*/
fid1r = fid1r - mean(fid1r);
fid1i = fid1r - mean(fid1i);
fid2r = fid2r - mean(fid2r);
fid2i = fid2r - mean(fid2i);

% %     /*------------------------------------------------------------*
% %      * - Determine average value of each channel of both FIDs
% %      *------------------------------------------------------------*/
avg1 = mean(fid1r + fid1i);
avg2 = mean(fid2r + fid2i);

% %     /*------------------------------------------------------------*
% %      * - Determine variance of each channel
% %      *------------------------------------------------------------*/
var1 = sum( (fid1r - avg1).^2 + (fid1i - avg1).^2 );
var2 = sum( (fid2r - avg2).^2 + (fid2i - avg2).^2 );

NoiseLevel1 = sqrt(var1 / (length(fid1r) + length(fid1i)) );
NoiseLevel2 = sqrt(var2 / (length(fid2r) + length(fid2i)) );

% % /*------------------------------------------------------------*
% %  * - Calculate Variance and Noise Figure from current values
% %  *------------------------------------------------------------*/

Var = NoiseLevel1^2 / NoiseLevel2^2;
Temp = 300 / 77;

if (NoiseLevel1 < NoiseLevel2) 
    Var = 1/Var;
    Temp = 1/Temp;
end
if f_temperature
    excessNoise = - 10.0*log10( 1 - Temp );
end

noiseRatio = - 10.0*log10( Var - 1);
NoiseFigure = excessNoise + noiseRatio;

if f_verbose
    fprintf( 'NoiseLevel (dataset #1) = %.2f\n', NoiseLevel1 );
    fprintf( 'NoiseLevel (dataset #2) = %.2f\n', NoiseLevel2 );
    fprintf( 'Ratio of Variances      = %f\n', Var );
    fprintf( 'Noise ratio             = %f dB\n', noiseRatio );
    if f_temperature
        fprintf( 'Ratio of Temperatures   = %f\n', Temp );
    end
    fprintf( 'Noise figure (ENR=%.2f dB)  = %f dB\n', excessNoise, NoiseFigure );
end

noisefig.NF = NoiseFigure;

%------------------------------------------------------------
% % /*:=INFO=:*******************************************************
% %  * AU: calc_noisefig
% %  *
% %  * Description :
% %  *   This AU program is used to calculate the noise figure from 
% %  *   either two individual EXPNOs or from the observed rms noise
% %  *   of two GSPs at different temperatures.
% %  *
% %  * Usage:
% %  *   Load the ONEPULSE protocol, set IMND_spect_acq_size equal to
% %  *   at least 8K.  Set the RG to nearly 100%.
% %  *   Connect a BNC cable to the output of the preamplifier to be
% %  *   measured.  The other end of the BNC should have a 50 ohm
% %  *   resister soldered between the cable and the ground.
% %  *   Record data in one of two manners:
% %  *   (1) Immerse resister in a cold solution of known temperature
% %  *       and record the area of the "FID".  Repeat for solution
% %  *       at warmer temperature.  Execute this AU program and enter
% %  *       these areas when asked for the 'rms Noise' for a dataset.
% %  *   (2) Immerse resister in a cold solution of known temperature
% %  *       and acquire an FID using GOP.  Clone the scan and repeat 
% %  *       for solution at warmer temperature.  Select initial EXPNO 
% %  *       and execute this AU program.  Enter the EXPNOs and 
% %  *       temperatures for each dataset.
% %  *   The Noise Figure will be output upon the conclusion of the AU.
% %  *   Typically, we measure at liquid N2 (77 K) and then at room 
% %  *   temperature (300 K).
% %  *
% %  * The original code for this AU was borrowed from subtract_fids
% %  *
% %  * Written by RER @ BII - Tue Aug  5 09:47:41 EDT 1997
% %  * Date of last modification: 
% %  *
% %  * --------------------------------------------------------------
% %  * Interface for AUTOMATIC EXECUTION MODE:
% %  *
% %  * [n_fig expno_1 expno_2 T1 T2 ]
% %  * [n_fig Noise1  Noise2  T1 T2 ]
% %  *
% %  * explanation:
% %  *
% %  * n_fig        f = calculate noise figure from FIDS of two datasets
% %  *              r = calculate noise figure from recorded RMS noise 
% %  *
% %  * if (n_fig == f):
% %  * expno_1      EXPNO of the first  dataset
% %  * expno_2      EXPNO of the second dataset
% %  *
% %  * if (n_fig == r):
% %  * Noise1       rms noise level of the first dataset
% %  * Noise2       rms noise level of the second dataset
% %  *
% %  * T1           Temperature of the first dataset
% %  * T2           Temperature of the second dataset
% %  * 
% %  *::=info=:******************************************************/
% %     /*------------------------------------------------------------*
% %      * - Separate FIDs in memory into real & imaginary components
% %      *------------------------------------------------------------*/
% % 
% %     i = j = 0;
% %     do
% %     {
% % 	fid1r[j] = fid1[i];
% % 	fid1i[j] = fid1[i+1];
% % 	fid2r[j] = fid2[i];
% % 	fid2i[j] = fid2[i+1];
% % 	i += 2;
% % 	j += 1;
% %     } 
% %     while (i < ACQ_size_1);
% % 
% %     /*------------------------------------------------------------*
% %      * - Remove DC offsets from each channel for both FIDs
% %      *------------------------------------------------------------*/
% % 
% %     sum1r = sum1i = sum2r = sum2i = 0.0;
% %     for (i=0; i<ACQ_size_1/2; i++)
% %     {
% % 	sum1r += (double)fid1r[i];
% % 	sum1i += (double)fid1i[i];
% % 	sum2r += (double)fid2r[i];
% % 	sum2i += (double)fid2i[i];
% %     }
% % 
% %     DC1r = sum1r / (double)(ACQ_size_1/2);
% %     DC1i = sum1i / (double)(ACQ_size_1/2);
% %     DC2r = sum2r / (double)(ACQ_size_1/2);
% %     DC2i = sum2i / (double)(ACQ_size_1/2);
% % 
% %     for (i=0; i<ACQ_size_1/2; i++)
% %     {
% % 	fid_1r[i] = (int)((double)fid1r[i] - DC1r);
% % 	fid_1i[i] = (int)((double)fid1i[i] - DC1i);
% % 	fid_2r[i] = (int)((double)fid2r[i] - DC2r);
% % 	fid_2i[i] = (int)((double)fid2i[i] - DC2i);
% %     }
% % 
% %     sum1 = sum2 = 0.0;
% %     for (i=0; i<ACQ_size_1/2; i++)
% %     {
% % 	sum1 = sum1 + (double)fid_1r[i] + (double)fid_1i[i];
% % 	sum2 = sum2 + (double)fid_2r[i] + (double)fid_2i[i];
% %     }
% % 
% %     /*------------------------------------------------------------*
% %      * - Determine average value of each channel of both FIDs
% %      *------------------------------------------------------------*/
% % 
% %     avg1 = sum1 / (double)ACQ_size_1;
% %     avg2 = sum2 / (double)ACQ_size_2;
% % 
% %     /*------------------------------------------------------------*
% %      * - Determine variance of each channel
% %      *------------------------------------------------------------*/
% %     
% %     var1 = var2 = 0.0;
% %     for (i=0; i<ACQ_size_1/2; i++)
% %     {
% % 	var1 += ((double)fid_1r[i]-avg1)*((double)fid_1r[i]-avg1);
% % 	var1 += ((double)fid_1i[i]-avg1)*((double)fid_1i[i]-avg1);
% % 	var2 += ((double)fid_2r[i]-avg2)*((double)fid_2r[i]-avg2);
% % 	var2 += ((double)fid_2i[i]-avg2)*((double)fid_2i[i]-avg2);
% %     }
% %     
% %     NoiseLevel1 = sqrt(var1 / (double)(ACQ_size_1-1));
% %     NoiseLevel2 = sqrt(var2 / (double)(ACQ_size_2-1));
% % }
% % 
% % /*------------------------------------------------------------*
% %  * - Calculate Variance and Noise Figure from current values
% %  *------------------------------------------------------------*/
% % 
% % if (temperature_1 < temperature_2) 
% % {
% %     Var = (NoiseLevel1 * NoiseLevel1) / (NoiseLevel2 * NoiseLevel2);
% %     Temp = (temperature_1 / temperature_2);
% % } 
% % else 
% % {
% %     Var = (NoiseLevel2 * NoiseLevel2) / (NoiseLevel1 * NoiseLevel1);
% %     Temp = (temperature_2 / temperature_1);
% % }
% % 
% % #if DEBUG_ON
% % printf( "NoiseLevel (dataset #1) = %.2f\n", NoiseLevel1 );
% % printf( "NoiseLevel (dataset #2) = %.2f\n", NoiseLevel2 );
% % printf( "Ratio of Variances = %.2f\n", Var );
% % printf( "Temperature (dataset #1) = %.2f degrees K\n", 
% % 	temperature_1);
% % printf( "Temperature (dataset #2) = %.2f degrees K\n", 
% % 	temperature_2);
% % printf( "Ratio of Temperatures = %.2f\n", Temp );
% % #endif
% % 
% % NoiseFigure = -10.0*log10((1.0-Var)/(1.0-Temp));
% % 
% % #if DEBUG_ON
% % printf( "calc_noisefig: Noise Figure = %.5f\n", NoiseFigure );
% % #endif
