function [spikeArr, spikeNum, img] = PVspikeDetectFid(dir, filenum, optin);
% %  [spikeArr, spikeNum, img] = PVspikeDetectFid(dir, filenum)
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
% %------------------------------
% %spike  artifacts removing    2-8-01
% % The spike in the k-space can add the strong noise to the image
% %The spike is detected via the comparasion with the std of the neighbor points.
% %data_in is five dimention raw data 
% 
% %This is the main program
% % the support subroutine programs are
% %		rd_raw_fid.m
% %		rdFid.m
% %		spike_flt_1d_std_66_0.m
% %		spike_flt_1d_std_66_2
% %		spike_flt_1d_std_intp.m
% %		wr_raw_fid
% %The current parameters are not including the percentage change for the mean. Here it is 0.05 
% % defined in 	spike_flt_1d_std_66_2.m
% %modified in Nov. 20, 2001
% %------------------------------
% %  Functions called: 
% %  Tested for: 
% %
% %  Aug 2003 -  Josef Pfeuffer
% %
FCTNAME = 'PVspikeDetectFid';
          
global STDPATH data_in info acqp

SAVENAMEAPP = '.despiked';

STDPATH.pv = '//wks5/josef/mridata_wks5/' 
dir = 'spikeJ02.lw1';
filenum = 15;
% STDPATH.pv = '//wks5/josef/mridata_wks4/' 
% dir = 'jptest.lw1';
% filenum = 7;
fidfile = 'fid';
f_read = 0;
f_despike = 1;

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.THRES    = 10;               % in Stdev's
dopt.NOISEROI = [1 10 1 10];
dopt.REPLACE  = 1;
dopt.WRITE    = 0;
dopt.VERBOSE  = 1;

      % --- arg handling
nargVars = 0;
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

  %the obtained the values
%   main_dir 	= char(answer(1)) 
%   thres = str2num(char(answer(2)))
%   out_filename = char(answer(3))
  

%load the raw data
%[fdata,ref_scan_loc] = rd_raw_fid(main_dir);

if f_read
    info = PVrdFid(dir, filenum, opt('FIDFILE',fidfile,'GETINFO',1,'VERBOSE',0));
    data_in = PVrdFid(dir, filenum, opt('FIDFILE',fidfile,'VERBOSE',f_verbose));
%%fid(x,y,k,n,t) = k-space point corresponding to readout step x, line y, slice k, segment n and time t
     data_in = data_in(65:192,1:128,:,:);
     info.nx=256;
     info.ny=128;
    
    nseg       = acqp.IMND_numsegments;			
    if strcmp(acqp.EPI_segmentation_mode,'No_Segments')    % glitch for EPI
        nseg   = 1;			
    end
    data_in = reshape(data_in,info.nx/2,info.ny/nseg,nseg,info.nslices,info.nr);
else
    
end
ref_scan = 0;
thres = 2.0;        %thres = 3 --> 99.7%
					%thres = 2.58 --> 99.0%
					%thres = 2.0  --> 95.0%

ss = size(data_in)
if length(ss) < 5
   disp('wrong data ')
   return
 end
    
%locate the spike position   and correct the spike with interpolation

amp_tmp = abs(double(data_in));
data_out = data_in;

if f_despike
	%obtain one standard deviation from he high k-space
	tmp1 = amp_tmp(ss(1)-2,ss(2)-2,1,1,:);
	[tmp_out, tmp_flg] = spike_flt_1d_std_66_0(tmp1,thres,ref_scan);
	std_val = std(tmp_out)
	
	for lp4=1:ss(4)
       for lp3 = 1:ss(3)
          fprintf('nseg<%d> slice<%d>\n',lp3,lp4);
          for lp2 = 1:ss(2)
             for lp1 = 1:ss(1)
                tmp1 = squeeze(amp_tmp(lp1,lp2,lp3,lp4,:));
                tmp_out = spike_flt_1d_std_66_2(tmp1,thres,ref_scan,std_val);
                
                tmp2 = squeeze(data_in(lp1,lp2,lp3,lp4,:));
                tmp = spike_flt_1d_std_intp(tmp2,tmp_out);
                data_out(lp1,lp2,lp3,lp4,:) = tmp;
  figure(2)
  plot(tmp_out)
  figure(1)
  plot(abs(double(tmp2)))
  figure(3)
  plot(abs(double(tmp)))
             end
          end
       end
	end
end

if f_write
    %%wr_raw_fid(data_out, main_dir, out_filename);
end

stop      
%------------------------------------------------------------
function [data_out, data_flg] = spike_flt_1d_std_66_0(data_in,thres,ref_scan)

%data_in is amplitude of one dimension

data_out = squeeze( data_in );
ss = size(data_out);
data_flg = zeros(ss);
tt_flt = medfilt1(data_out,10);

for lp5=10:ss(1)
    std_val = std(data_out(lp5-9:lp5-1));
    std_mean = mean(data_out(lp5-9:lp5-1));	%added  3-5-01
   
    if  abs(data_out(lp5)-tt_flt(lp5)) > thres * std_val
        %if lp5 ~= ref_scan & abs(data_out(lp5)-std_mean) > 0.1*std_mean+thres * std_val    %changed 3-5-01
        
        data_out(lp5) = mean(data_out(lp5-2:lp5-1));
        data_flg(lp5) = 1; 
    end
end
%------------------------------------------------------------
function data_out = spike_flt_1d_std_66_2(data_in,thres,ref_scan,std_val)

%data_in is amplitude of one dimention

ss = size(data_in);
data_out = data_in;
data_flg = zeros(ss);
tt_flt = medfilt1(data_in,10);

for lp5=10:ss(1)
    std_mean = mean(data_out(lp5-9:lp5-1));	%added  3-5-01
    if lp5 ~= ref_scan & abs(data_out(lp5)-std_mean) > 0.05*std_mean+thres * std_val    %changed 3-5-01
        data_out(lp5) = mean(data_out(lp5-2:lp5-1));
        data_flg(lp5) = 1; 
    end
end
data_out = data_flg;
%------------------------------------------------------------
function data_out = spike_flt_1d_std_intp(data_in,flag)

%data_in is original data of one dimention
%

ss = size(data_in);

data_out = data_in;

for lp5=4:ss(1)
   if flag(lp5)
      tmp_r = mean(real(data_out(lp5-2:lp5-1)));
      tmp_i = mean(imag(data_out(lp5-2:lp5-1)));
     data_out(lp5) = tmp_r + i*tmp_i;
    end
 
 end
%------------------------------------------------------------
%------------------------------------------------------------

%------------------------------------------------------------
function PVspikeDetectDISABLED

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
