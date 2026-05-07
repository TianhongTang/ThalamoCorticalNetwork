function DATA = dorkTrace2( dir, filenum, optin )
% %  DATA = dorkTrace2( dir, filenum, optin )
% %
% %  Analysis of Dynamic Off-Resonance Effects in K-space (DORK)
% %      J. Pfeuffer, P.-F. van de Moortele, K. Ugurbil,
% %          X. Hu, G. H. Glover, Magn.Res.Med 2001
% %
% %  default options:
% %     opt(  'NS',[0 0],       % slice range
% %           'NR',[0 0],       % range in time series
% %           'VERBOSE','1')
% %
% %  return: series of center k-space data 
% %               (in PRECISION type: NOT DOUBLE !!! -> faster)
% %          global acqp/reco parameters
% %
% %  Tested for: onepulse, EPI
% %
% %  Sep 2001 -  Josef Pfeuffer
% %  BP: - returning structure DATA now
% %      - plotting PSD in figure 3
% %
FCTNAME = 'dorkTrace2';

global acqp

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.NS = [0 0];
dopt.NR = [0 0]; 
dopt.DETREND = 0;
dopt.SLICE = 0;
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
flagDetrend = dopt.DETREND;
slicenum    = dopt.SLICE;
f_verbose   = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

   % read a single FID reference image + ACQP parameter: => int32 type, 4D
if nrrange(1) == 0
   tdat = PVrdFid2(dir, filenum, opt('NR',[1 1],'FUDGE',dopt.FUDGE,'VERBOSE',f_verbose) );      
else
   tdat = PVrdFid2(dir, filenum, opt('NR',[nrrange(1) nrrange(1)],'FUDGE',dopt.FUDGE,'VERBOSE',f_verbose) );  
end
               
if (acqp.EPI_TE_eff > 0)
   te = acqp.EPI_TE_eff;
else
   goo = input('echo time TE/ms [50] ');
   if isempty(te) te = 50; else  te = goo(1); end
end
te = te*0.001
sw = acqp.SW_h
if strcmpi(acqp.EPI_segmentation_mode, 'No_Segments')
   nseg = 1;
else
   nseg = acqp.IMND_numsegments
end
ref_pos= 1;       % position of non-PE image
nnav= 10;         % only used for DTYPE = 2
navecho = 0;      % navecho lines at beginning
DTYPE = 1;        % 1: EPI, 2: SPIRAL
NIREF = max([1 ref_pos-5]);	 % position of phase reference
nslices = acqp.NSLICES;
nr = acqp.NR;
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

s_tdat = size(tdat)
if (DTYPE == 1)
   dim1 = s_tdat(1);
   dim2 = s_tdat(2)/nseg	;
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

%%input('===> ');      % display kspace

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
% --- slice handling here omitted: first slice taken
slicenum = min([nslices slicenum]);
if (slicenum < 1)
   islice = round(nslices/2);
else
   islice = slicenum;
end
if (DTYPE == 1)
   temulti = zeros(nseg);
   navdata = complex( zeros(nslices_12, nr_12, nseg) );  
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
      fprintf('kspace max at (%d,%d) %d\n',nx_nav, ny_nav, k0ind);
      temulti(iseg) = tacq(nx_nav, ny_nav0) + te;      
      navdata(:,:, iseg) = squeeze( PVrdFid2(dir, filenum, ...
         opt('NS',[ns1 ns2],'NR',[nr1 nr2],'XY',[nx_nav ny_nav],'FUDGE',dopt.FUDGE,'VERBOSE',f_verbose) )); 
   end
else    
   stop
   navmulti = 1;
   temulti = te;
   navdata = reshape( double(tdat(:,:,:,1:nnav)), nseg*ni, navmulti, nnav );
   % all segments equally treated
end
expLabel = sprintf('<%s> exp#%d  sl#%d  TR=%.1f s', dir, filenum, islice, acqp.IMND_rep_time);

for imulti=1:nseg
   navtmp = squeeze(navdata(islice,:,imulti));
   TEnav = temulti(imulti)
   Snav = navtmp;
   Sref = avg( Snav );
   navcor = ( Sref/abs(Sref) ).*( conj(Snav)./abs(Snav) );
   df = unwrap(angle(navcor))/TEnav/(2*pi);

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
      %ylim([-1,1])
      ylabel('D amplitude [%]')
      title(titlestr1)
      
      subplot(2,1,2)
      plot( df )
      %ylim([-1,1])
      ylabel('D offres [Hz]')
      title(titlestr2)
      xlabel(expLabel);
      end
   end
   if (nseg > 1 & imulti < nseg & f_verbose > 0) input('===> '); end
end
dat = navdata;

% BP: instead of dat, return structure DATA 
DATA.deltaAmp=(abs(Snav)/mean_a-1)*100;
DATA.deltaFreq=(df);
DATA.acqp=acqp;

% BP: Display Power Spectral Density
rep=acqp.IMND_rep_time*nseg;
repFreq=1/rep;
nr=acqp.NR;
%nr=nr_12;

Amp=fft(detrend(DATA.deltaAmp,nr));
Ores=fft(detrend(DATA.deltaFreq,nr));

psdOres=abs(Ores.*conj(Ores)/nr);
psdAmp=abs(Amp.*conj(Amp)/nr);
f=repFreq*(0:ceil(nr/2))/nr;
whos f
whos psdOres

if f_verbose
figure(3);
subplot(2,1,1)
plot(f,psdAmp(1:(length(f))))
ylim([0,1])
title('FFT for Amplitude and k-center Offresonance')
ylabel('PSD Amplitude')
xlabel('Frequency [Hz]')

subplot(2,1,2)
plot(f,psdOres(1:(length(f))))
ylim([0,1])
ylabel('PSD Offresonance');
xlabel(expLabel);
end

DATA.f = f;
DATA.psdAmp = psdAmp;
DATA.psdOres = psdOres;

