function dat = dorkTrace( dir, filenum, optin )
% %  dat = dorkTrace( dir, filenum, optin )
% %
% %  Analysis of Dynamic Off-Resonance Effects in K-space (DORK)
% %      J. Pfeuffer, P.-F. van de Moortele, K. Ugurbil,
% %          X. Hu, G. H. Glover, Magn.Res.Med 2001
% %
% %  default options:
% %     opt(  'NS',[0 0],       % slice range
% %           'NR',[0 0],       % range in time series
% %           'REFNUM',1,       % reference FID #
% %           'FIDFILE','fid',  % fid/fid.orig/ser
% %           'DETREND',0,
% %           'FUDGE' = [0 0 0 0];
% %           'VERBOSE','1')
% %
% %  return: dat.deltaAmp  (nslices, nsegments, nrepetitions);
% %          dat.deltaFreq (nslices, nsegments, nrepetitions);
% %          dat.acqp      (ACQP struct)
% %          global acqp/reco parameters
% %
% %  Tested for: onepulse, EPI
% %
% %  Sep 2001 -  Josef Pfeuffer
% %  Oct 2002 - JP: opt DETREND,SLICE,FUDGE
% %
FCTNAME = 'dorkTrace';

global acqp navFidDat navDat

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.NS = [0 0];
dopt.NR = [0 0]; 
dopt.REFNUM = 1;
dopt.FIDFILE = 'fid';
dopt.DETREND = 0;
dopt.FUDGE = [0 0 0 0];
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

nsrange     = dopt.NS;
nrrange     = dopt.NR; 
refNumFid   = max(max(dopt.REFNUM,1),nrrange(1));
flagDetrend = dopt.DETREND;
f_verbose   = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

info = PVrdFid(dir, filenum, opt('FIDFILE',dopt.FIDFILE,'GETINFO',1,'VERBOSE',f_verbose));
   % read a SINGLE reference FID and ACQP parameters: => int32 type, 4D
tdat = PVrdFid(dir, filenum, opt('NR',[refNumFid refNumFid],'FIDFILE',dopt.FIDFILE, ...
    'FUDGE',dopt.FUDGE,'VERBOSE',f_verbose) );  

               
if (acqp.EPI_TE_eff > 0)
   te = acqp.EPI_TE_eff;
else
   goo = input('echo time TE/ms [50] ');
   if isempty(te) te = 50; else  te = goo(1); end
end
te = te*0.001;
sw = acqp.SW_h;
fprintf('EPI: te = %.2f ms, sw = %.3f kHz\n', te*1e3, sw*1e-3);
if strcmpi(acqp.EPI_segmentation_mode, 'No_Segments')
   nseg = 1;
else
   nseg = acqp.IMND_numsegments;
end
ref_pos= 1;       % position of non-PE image
nnav= 10;         % only used for DTYPE = 2
navecho = 0;      % navecho lines at beginning
DTYPE = 1;        % 1: EPI, 2: SPIRAL
NIREF = max([1 ref_pos-5]);	 % position of phase reference
nslices = info.nslices;
nr = info.nr;
if (nsrange(1) ~= 0)
    ns1 = max( [nsrange(1) 1] );
    ns2 = min( [nsrange(2) nslices] );
    if (ns2 < ns1) ns2 = nslices; end
else
    ns1 = 1;
    ns2 = nslices;
end
if (nrrange(1) ~= 0)
    nr1 = max( [nrrange(1) 1] );
    nr2 = min( [nrrange(2) nr] );
    if (nr2 < nr1) nr2 = nr; end
 else
    nr1 = 1;
    nr2 = nr;
end
nslices_12 = ns2 - ns1 + 1;
nr_12 = nr2 - nr1 + 1;

if f_verbose
figure(2);
colormap(gray);
subplot(1,1,1);
plotdat = double(tdat(:,:,1,1));
imagesc(abs(plotdat));
title('k-space ABS()')
%%axis image;
end

s_tdat = size(tdat);
if (DTYPE == 1)
   dim1 = s_tdat(1);
   dim2 = s_tdat(2)/nseg;
   dim3 = nseg;
   dim4 = prod(s_tdat)/(dim1*dim2*dim3);
   fprintf('max = %8.2f,  min = %8.2f \n', ...
            maxall(abs(double(tdat))), minall(abs(double(tdat))) );
   tdat = reshape(tdat, dim1, dim2, dim3, dim4);
else 
   dim1 = s_tdat(1);
   dim2 = 1;
   dim3 = nseg;
   dim4 = prod(s_tdat)/(dim1*dim2*dim3);
   fprintf('max = %8.2f,  min = %8.2f\n', maxall(abs(tdat)), minall(abs(tdat)));
   tdat = reshape(tdat, dim1, dim2, dim3, dim4);
end
nx=dim1;
ny=dim2;
%% dim3=nseg;
ni=dim4;

%%key;      % display kspace

% ----- determine timing
if (DTYPE == 1)
   dwelltime = (1:nx)'./sw;      % timing of single readout line
                                 % here assumed: same DW, linear sampling
   tline = dwelltime( fix(sum(size(dwelltime))/2) );  %.32
   tacq = zeros(nx, ny);
   for ipe=1:ny
      tacq(:,ipe) = ( dwelltime + (ipe-ny/2-navecho+0.5)*tline );
   end
   %%%tacq = (reshape(1:nx*ny, nx, ny)) / sw;
else 
   tacq = (1:nx- nnav )/sw;
end

% ----- extract navigator data (nnav,ni,navmulti)
if (DTYPE == 1)
   temulti = zeros(nseg);
   navdata = complex( zeros(nslices_12, nseg, nr_12) );  
   for islice=ns1:ns2
   for iseg=1:nseg
      imgtmp1 = abs( double(tdat(:,:,iseg,islice)) );
      imgtmp1(:,1) = 0;		% don't look at navigator in first line
      
      maxind = find( imgtmp1 >= maxall(imgtmp1));
      if (~isempty(maxind))
         k0ind = max(maxind);
      else 
         error('max NOT found');
      end
      nx_nav = mod(k0ind, nx);
      ny_nav0 = fix(k0ind / nx)+1;
      ny_nav = ny_nav0 + ny*(iseg-1);
      if f_verbose
          fprintf('kspace max at (%d,%d) %d   (seg#%d  sl#%d)\n',nx_nav, ny_nav, k0ind, iseg, islice);
      end
      temulti(iseg) = tacq(nx_nav, ny_nav0) + te;      
      navdata(islice,iseg,:) = squeeze( PVrdFid(dir, filenum, ...
         opt('NS',[islice islice],'NR',[nr1 nr2],'XY',[nx_nav ny_nav],'FIDFILE',dopt.FIDFILE,...
         'FUDGE',dopt.FUDGE,'VERBOSE',0) )); 
   end
   end
else    
   stop   %% spiral NOT YET TESTED/IMPLEMENTED
   navmulti = 1;
   temulti = te;
   navdata = reshape( double(tdat(:,:,:,1:nnav)), nseg*ni, navmulti, nnav );
   % all segments equally treated
end

     % create dat struct (RETURN struct)
dat.deltaAmp = zeros(nslices_12, nseg, nr_12);
dat.deltaFreq = zeros(nslices_12, nseg, nr_12);
dat.acqp = acqp;
%%dat.PSDfreq = zeros(ceil(nr_12/2)+1);
%%dat.PSDamp = zeros(nslices_12, nseg, ceil(nr_12/2)+1);
%%dat.PSDoffres = zeros(nslices_12, nseg, ceil(nr_12/2)+1);


for islice=ns1:ns2
expLabel = sprintf('<%s> exp#%d  sl#%d  TR=%.3f s', dir, filenum, islice, acqp.IMND_rep_time);
for iseg=1:nseg
   navtmp = squeeze(navdata(islice,iseg,:));
   TEnav = temulti(iseg);
   Snav = navtmp;
   %%Sref = avg( Snav );
   Sref = Snav(refNumFid);       %% THIS is the reference point to compare different SLICE/SEGMENTS
   navcor = ( Sref/abs(Sref) ).*( conj(Snav)./abs(Snav) );
   df = unwrap( angle(navcor) )/TEnav/(2*pi) ;
   fprintf('%s: TEnav = %.2f ms (seg#%d  sl#%d)\n', FCTNAME, TEnav*1e3, iseg, islice);

   if f_verbose
       figure(3);
       subplot(1,1,1)
       [pmax j pmin j]= maxall([abs(navtmp),abs(real(navtmp))]);
       plot( abs(real(navtmp) ))
       xlim([0 min([nr_12, 30])]);
       ylim([pmin pmax]);
       hold on
       plot( abs(navtmp) )
       hold off
   end
   
   if (nr_12 > 1)
      if flagDetrend
          std_a = std(detrend(abs(Snav)));
          mean_a = mean(abs(Snav));
          std_ph = std(detrend(angle(Snav)));    % rad
          mean_ph = mean(angle(Snav));  % rad
          std_f = std(detrend(df));
          mean_f = mean(df);
          titlestr1 = sprintf('SNAV: mean=%f;  stdev=%f;  stdev%%(detrended)=%f %%', mean_a, std_a, std_a/mean_a*100);
          titlestr2 = sprintf('NAVIGATOR(detrended): std=%f Hz;   PHASE: %f (%f) deg', std_f, mean_ph*180/pi, std_ph*180/pi);
      else
          std_a = std(abs(Snav));
          mean_a = mean(abs(Snav));
          std_ph = std(angle(Snav));    % rad
          mean_ph = mean(angle(Snav));  % rad
          std_f = std(df);
          mean_f = mean(df);          
          titlestr1 = sprintf('SNAV: mean=%f;  stdev=%f;  stdev%%=%f %%', mean_a, std_a, std_a/mean_a*100);
          titlestr2 = sprintf('NAVIGATOR: std=%f Hz;   PHASE: %f (%f) deg', std_f, mean_ph*180/pi, std_ph*180/pi);
      end
      
      if f_verbose
      figure(1);
      subplot(2,1,1)
      plot( (abs(Snav)/mean_a-1)*100 )
      %ylim([-4,4])
      ylabel('delta amplitude [%]')
      title(titlestr1)
      
      subplot(2,1,2)
      plot( df )
      %ylim([-6,6])
      ylabel('delta offres [Hz]')
      title(titlestr2)
      xlabel(expLabel);
      end
   end

     % return structure: dat (continued BELOW!)
deltaAmp = (abs(Snav)/mean_a-1)*100;
dat.deltaAmp(islice,iseg,:) = deltaAmp;
dat.deltaFreq(islice,iseg,:) = df;

    % Power Spectral Density
rep = acqp.IMND_rep_time*nseg;
repFreq = 1/rep;
repFreqArr = repFreq*(0:ceil(nr/2))/nr;
AmpPSD = fft(detrend(deltaAmp),nr);
FreqPSD = fft(detrend(df),nr);
AmpPSD=abs(AmpPSD.*conj(AmpPSD)/nr);
FreqPSD=abs(FreqPSD.*conj(FreqPSD)/nr);

     % return structure: dat 
%%dat.PSDfreq = repFreqArr;
%%dat.PSDamp(islice,iseg,:) = AmpPSD(1:length(repFreqArr));
%%dat.PSDoffres(islice,iseg,:) = FreqPSD(1:length(repFreqArr));

if f_verbose
figure(3);
subplot(2,1,1)
plot(repFreqArr,AmpPSD(1:(length(repFreqArr))))
ylim([0,1])
title('FFT for Amplitude and k-center Offresonance')
ylabel('PSD amplitude')
xlabel('Frequency [Hz]')

subplot(2,1,2)
plot(repFreqArr,FreqPSD(1:(length(repFreqArr))))
ylim([0,1])
ylabel('PSD offres');
xlabel(expLabel);
end

   if (nseg*(ns2-ns1+1) > 1 & islice*iseg < nseg*ns2 & f_verbose > 0) key; end
end    %for iseg=1:nseg
end    %for islice=ns1:ns2

