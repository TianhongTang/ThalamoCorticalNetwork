function saccs = saccadeOKNAnalysis1D (em_list,samp_time,resamp_time)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  FUNCTION:
%    saccadeAnalysis1D(emlist, samp_time, resamp_time) 
%
%  PURPOSE: 
%    This is a generic eye movement analysis function that can be used
%    to identify most saccades.  It is relatively robust to fixational
%    movements, as well as nonzero intersaccadic eye velocities (e.g. OKN).
%    This program does checks a single eye trace (say, vertical) and 
%    extracts saccades in several steps:  
% 		1) Velocity thresholding
%		2) Check for minimum amplitude,interval
%               3) Finds characteristics of fixation periods
%               4) Returns structure with relevant data
%
%   Since it is not possible to make a package that is flexible enough to
%   handle every eye movement case, there are some changeable parameters
%   (capital letters, below) that can be used to 'tune' the performance of 
%   saccade detection, etc.  This will probably be necessary in the beginning.
%   I recommend setting the SHOW_FIGURES to 1 so one can see some relevant
%   distributions, etc.
%
%
%   D.Leopold  7-APR-00
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
   error('USAGE: saccadeAnalysis1D(em [samp_time] [resamp_time])');
end
if nargin < 2
  samp_time = 5;     % milliseconds	
end
if nargin < 3
  resamp_time = 1;   % milliseconds	
end

saccs = [];                     % initialize the return structure

VEL_THRESHOLD          = 40.0;   % Thresh velocity for "candidate" saccades
FORM_MS_HWIN	       = 60.0;  % how far from saccade begin/end to consider
SLOPE_MS_HWIN	       = 120.0; % how far from saccade begin/end to consider
BLINK_HWIN	       = 120.0; % how far from saccade begin/end to consider
N_MONOTONIC	       = 30;     % number of continuous ms to identify saccade
VELCONVSTD	       = 3.0;   % milliseconds velocity gaussian smoothing
MINIMUM_AMPLITUDE      = 1;  % smallest detectable amplitude for the algorithm (used to be 0.2)
MINIMUM_INTERVAL       = 100.0;  % smallest detectable interval
%MINIMUM_RESIDUAL       = 0.5;   % residual for the linear regression for fixation
MINIMUM_RESIDUAL       = 4.5;   % residual for the linear regression for fixation
MINIMUM_SLOPE_DIFF     = 0.002; % diff presacc slope and putative 'saccade' slope
CORRECTIVE_SACC_TIME   = 50;    % time in which to consider a saccade 'corrective'
BLINK_VEL_THRESH       = 500;	% deg/sec, velocity threshold for detecting blinks
BLINK_DISP_THRESH      = 100.0;	% deg/sec, displacement threshold for detecting blinks
SHOW_FIGURES	       = 1;     % display some raw traces, distributions, etc
EARLIEST_SACCADE       = 100;   % forget stuff in the very beginning of the trial
MIN_AMP_DUR_RATIO      = 0.005; % do not accept low amplitude saccades with long durations...


% Fit a spline to the eye movement data
%
em = spline(samp_time*[1:length(em_list)], em_list,...
	 resamp_time*[1:(samp_time/resamp_time)*length(em_list)]);
% Figure out the instantaneous velocities, convert to deg/sec
%
emvel = diff(em)*1000.;

% Velocity threshold -- this is a free parameter
k = normpdf([-5:5],0,VELCONVSTD);
absvel =  doConv(abs(emvel),k);
bigvels = absvel >= VEL_THRESHOLD;


% at this point, bigvels is a list of lists with ones 
% indicating that a saccade candidate is present
%
% Now we figure out the beginning and ending times for the 
% candidate saccades...these give us indices
%
diffvels = diff(bigvels);
upindxs  = find(diffvels == 1);
downindxs  = find(diffvels == -1);


% make sure there is no frameshift error here
%
if length(upindxs) == 0 | length(downindxs) == 0
   return
end
firstup   = upindxs(1);
firstdown = downindxs(1);
offset    = (firstup > firstdown);
nups	  = length(upindxs);
ndowns    = length(downindxs)-offset;
numsacc   = min(nups,ndowns);
if numsacc == 0
   return
end


% Do roughly the same thing for blinks.
realbigdisp    = abs(em) > BLINK_DISP_THRESH;
if sum(realbigdisp)
  blinkdiffvels  = diff(realbigdisp);
  blinkupindxs   = find(blinkdiffvels == 1);
  blinkdownindxs = find(blinkdiffvels == -1);
  if length(blinkupindxs) & length(blinkdownindxs)
     blinkfirstup   = blinkupindxs(1);
     blinkfirstdown = blinkdownindxs(1);
     blinkoffset    = (blinkfirstup > blinkfirstdown);
     blinknups      = length(blinkupindxs);
     blinkndowns    = length(blinkdownindxs)-blinkoffset;
     numblinks      = min(blinknups,blinkndowns);
     whichblinkindxs = [1:numblinks];
     blinkups        = blinkupindxs(whichblinkindxs);
     blinkdowns      = blinkdownindxs(whichblinkindxs+blinkoffset);
  else 
     numblinks 	     = 0;
  end  
else 
  numblinks 	 = 0;
end
OK_blinkindx  = [];
blinkamps      = [];
for i=1:numblinks
  tmpblinkamps = abs(em(blinkups(i):blinkdowns(i)));
  OK_blinkindx = [OK_blinkindx maxIndex(tmpblinkamps)+blinkups(i)];
  blinkamps    = [blinkamps max(tmpblinkamps)]; 	
end	


% again, make sure you only take the valid saccades by making 
% this selection (choosing) list
%
whichindxs = [1:numsacc];
ups        = upindxs(whichindxs);
downs      = downindxs(whichindxs+offset);


OK_indx  = [];
peakvels = [];
for i=1:length(ups)
  tmpvel   = absvel(ups(i):downs(i));
  OK_indx  = [OK_indx maxIndex(tmpvel)+ups(i)];
  peakvels  = [peakvels max(tmpvel)]; 	
end	



% based on the peak velocity, look before and after saccade 
% by a certain time 

hdur                  = round(10 + peakvels*0.05);  % msec
befores               = ups-hdur;
befores(befores < 1)  = 1;
beforepos             = em(befores);
afters                = downs+hdur;
maxind                = length(em);
afters(afters>maxind) = maxind;
afterpos              = em(afters);
amps 		      = abs(beforepos-afterpos);

% now have info for all the saccade candidates;
%
tmpsaccs.em            	= em;
tmpsaccs.emvel         	= emvel;

tmpsaccs.C_amps         = amps;
tmpsaccs.C_times	= OK_indx;
tmpsaccs.C_start        = befores;
tmpsaccs.C_stop	      	= afters;
tmpsaccs.C_intervals    = [diff(tmpsaccs.C_times) 100000];

% figure out some sacc stats, and save waveforms

formlen = 2*FORM_MS_HWIN+1; 
allforms = zeros(formlen,length(tmpsaccs.C_times));

for i = 1:length(tmpsaccs.C_times)
   start       = round(max(1,tmpsaccs.C_times(i)-FORM_MS_HWIN));	
   stop        = round(min(maxind,tmpsaccs.C_times(i)+FORM_MS_HWIN));	
   form	       = tmpsaccs.em([start:stop]);	
   if (length(form) < formlen)
      form = [form zeros(1,formlen-length(form))];
   end
   form = form-form(1);	
   allforms(:,i) = form';
end
tmpsaccs.C_peakvels   = tmpsaccs.emvel(tmpsaccs.C_times);
tmpsaccs.C_sign       = sign(tmpsaccs.C_peakvels);
tmpsaccs.C_forms      = allforms;    	

%
% use monotonicity to find real start and stop points, as well as eliminate any saccades that 
% look a little 'fishy';
%
k = ones(1,N_MONOTONIC)/N_MONOTONIC;
for i=1:length(tmpsaccs.C_times)
   % this first part is a monotonicity criterion...convolution with a 'plateau' 	
   signform 		= sign(diff(tmpsaccs.C_forms(:,i)))';
   mono 		= doConv(signform,k);
   mono			= round(mono*100)/100;
   offs 		= find(abs(mono)>=1.0);
   if length(offs)
      startoff		= round(offs(1)-N_MONOTONIC/2);     	 % take the first one
      stopoff		= offs(length(offs))+N_MONOTONIC/2;
      tmpsaccs.C_start(i)  = round(tmpsaccs.C_times(i)-FORM_MS_HWIN+startoff);
      tmpsaccs.C_stop(i)   = round(tmpsaccs.C_times(i)-FORM_MS_HWIN+stopoff);
      S_selmono(i)      = 1;	
   else
      tmpsaccs.C_start(i)  = 1;   		% fake values that will be selected out...
      tmpsaccs.C_stop(i)   = 2;
      S_selmono(i)      = 0;	
   end
   if tmpsaccs.C_start(i) < 1
      tmpsaccs.C_start(i)  = 1;
      tmpsaccs.C_stop(i)   = 2;
      S_selmono(i)         = 0;	
   end	
   if tmpsaccs.C_stop(i) > length(tmpsaccs.em) 
      tmpsaccs.C_start(i)  = 1;
      tmpsaccs.C_stop(i)   = 2;
      S_selmono(i)        = 0;	
   end	
end
tmpsaccs.C_durs = tmpsaccs.C_stop-tmpsaccs.C_start;

%
% Extract candidate pre- and post- saccadic fixation periods to help figure out if 
% a saccade has occurred
%
listlen = SLOPE_MS_HWIN-25+1;
alllist0 = zeros(listlen,length(tmpsaccs.C_times));
alllist1 = zeros(listlen,length(tmpsaccs.C_times));

for i = 1:length(tmpsaccs.C_times);
   if tmpsaccs.C_times(i) <= SLOPE_MS_HWIN
      C_prefixsel(i) = 0;
   else 
      C_prefixsel(i) = 1;
      beg0       = round(tmpsaccs.C_times(i)-SLOPE_MS_HWIN);	
      end0       = round(tmpsaccs.C_times(i)-25);
      beg1       = round(max(1,tmpsaccs.C_times(i)+25));
      end1       = round(min(maxind,tmpsaccs.C_times(i)+SLOPE_MS_HWIN));	
      list0      = tmpsaccs.em([beg0:end0]);
      if (length(list0) < listlen)
         list0 = [list0 zeros(1,listlen-length(list0))];
      end
      list0 = list0-list0(1);	
      alllist0(:,i) = list0';	
      list1      = tmpsaccs.em([beg1:end1]);
      if (length(list1) < listlen)
         list1 = [list1 zeros(1,listlen-length(list1))];
      end
      list1 = list1-list1(1);	
      alllist1(:,i) = list1';	
   end
end
tmpsaccs.C_presaccfix  = alllist0;
tmpsaccs.C_postsaccfix = alllist1;


%
% Fit a linear regression to each of the pre- and post- periods.  Find
% residuals.
% 
t = 1:length(tmpsaccs.C_presaccfix(:,1));
allslope0 = zeros(1,length(tmpsaccs.C_times));
allslope1 = allslope0;
allresid0 = allslope0;
allresid1 = allslope1;
for i = 1:length(tmpsaccs.C_times)
  [p0,s0] = polyfit(t,tmpsaccs.C_presaccfix(:,i)',1);   
  allslope0(i) = p0(1);
  allresid0(i) = s0.normr;
  [p1,s1] = polyfit(t,tmpsaccs.C_postsaccfix(:,i)',1);   
  allslope1(i) = p1(1);
  allresid1(i) = s1.normr;
end
tmpsaccs.C_presaccfixslope  = allslope0;   
tmpsaccs.C_postsaccfixslope = allslope1;   
tmpsaccs.C_presaccresid     = allresid0;   
tmpsaccs.C_postsaccresid    = allresid1;   

%
% Figure out if the slope of the line connecting the before and after points
% is significantly different from the before and after slopes.

for i = 1:length(tmpsaccs.C_times)
   startindx   = tmpsaccs.C_start(i);
   stopindx    = tmpsaccs.C_stop(i);
   start       = tmpsaccs.em(startindx);
   stop        = tmpsaccs.em(stopindx);
   allslope(i) = (stop-start)/(stopindx-startindx);
end
tmpsaccs.C_transaccslope    = allslope;


%
%  Finally determine the valid saccade candidates
% 

S_selamp  	         = tmpsaccs.C_amps > MINIMUM_AMPLITUDE;
S_selearly 	         = tmpsaccs.C_times > EARLIEST_SACCADE;
S_stable                 = tmpsaccs.C_presaccresid < MINIMUM_RESIDUAL...
			   & tmpsaccs.C_postsaccresid < MINIMUM_RESIDUAL;
S_slopesel	         = abs(tmpsaccs.C_transaccslope-...
			   tmpsaccs.C_presaccfixslope) > MINIMUM_SLOPE_DIFF;
S_ampdursel		 = (tmpsaccs.C_amps./(tmpsaccs.C_durs-N_MONOTONIC+1)) > MIN_AMP_DUR_RATIO;
N_selamp = sum(S_selamp);
N_selearly = sum(S_selearly);
N_selstable = sum(S_stable);
N_selslope = sum(S_slopesel);
N_ampdursel = sum(S_ampdursel);


%S_sel		         = find(S_selearly & S_selamp & S_selmono & S_stable);
%S_sel		         = find(S_selearly & S_selamp & S_selmono & S_stable & S_slopesel...
%	 & S_ampdursel);
S_sel		         = find(S_selearly & S_selamp & S_selmono & S_ampdursel);


saccs.S_amps		 = tmpsaccs.C_amps(S_sel);
saccs.S_times	      	 = tmpsaccs.C_times(S_sel);
saccs.S_peakvels         = tmpsaccs.C_peakvels(S_sel);
saccs.S_start            = tmpsaccs.C_start(S_sel);
saccs.S_stop             = tmpsaccs.C_stop(:,S_sel);
tmp_S_forms              = tmpsaccs.C_forms(:,S_sel);
for i=1:length(saccs.S_start)
  indxs = [saccs.S_start(i):saccs.S_stop(i)];
  saccs.S_forms{i} = tmpsaccs.em(indxs);
end
saccs.S_corrective       = [10000 diff(saccs.S_times)] < CORRECTIVE_SACC_TIME;
saccs.S_blinktimes	 = OK_blinkindx;
saccs.S_blinkamps	 = blinkamps;
saccs.S_em               = tmpsaccs.em;
saccs.S_emvel            = tmpsaccs.emvel;

%
% use the monotonicity on the blink velocity to identify starts and stops
%
if numblinks
  k = ones(1,N_MONOTONIC)/N_MONOTONIC;
  for i=1:length(saccs.S_blinktimes)
     % this first part is a monotonicity criterion...convolution with a 'plateau' 	
     bt  	        = saccs.S_blinktimes(i);
     bbt		= bt-BLINK_HWIN;
     bbt 		= max(bbt,BLINK_HWIN)+1;	
     ebt		= bt+BLINK_HWIN;
     ebt		= min(ebt,length(saccs.S_em));	
     signform 		= sign(diff(saccs.S_em([bbt:ebt])));
     mono 		= doConv(signform,k);
     mono		= floor(mono*100)/100;
     startoffs 		= find(abs(mono)==1.0);
     stopoffs  		= find(abs(fliplr(mono))==1.0);
     if length(startoffs) & length(stopoffs)
        startoff		= startoffs(1);    % take the first one
        stopoff		= stopoffs(1);
        saccs.S_blinkstart(i)  = saccs.S_blinktimes(i)-BLINK_HWIN+startoff-N_MONOTONIC;
        saccs.S_blinkstop(i)   = saccs.S_blinktimes(i)+BLINK_HWIN-stopoff+N_MONOTONIC;
     else
        saccs.S_blinkstart(i)  = saccs.S_blinktimes(i)-70;
        saccs.S_blinkstop(i)   = saccs.S_blinktimes(i)+70;
     end
     saccs.S_blinkstart(saccs.S_blinkstart < 1) = 1;	
     saccs.S_blinkstop (saccs.S_blinkstop > length(saccs.S_em)) = length(saccs.S_em);	
  end
else 
     saccs.S_blinkstart = [];
     saccs.S_blinkstop  = [];
end


if SHOW_FIGURES

 figure(1)
 clf
 for i=1:length(saccs.S_forms)
  tmp = saccs.S_forms{i};	
  plot(tmp-tmp(1));
  hold on
 end
 set(gca,'YLim',[-1.2 1.2]);

 figure(2)
 bins = [-0.01:0.0002:0.01];
 plot(bins,hist(tmpsaccs.C_presaccfixslope,bins),'r');
 hold on
 plot(bins,hist(tmpsaccs.C_postsaccfixslope,bins),'b');
 title('PostSaccFixSlope');

 figure(3)
 bins = [-0.00:0.005:1.5];
 plot(bins,hist(tmpsaccs.C_presaccresid,bins),'r');
 hold on
 plot(bins,hist(tmpsaccs.C_postsaccresid,bins),'b');
 title('Residuals');

 figure(4)
 bins = [0.00:0.0006:0.03];
 bar(bins,hist(abs(tmpsaccs.C_transaccslope-tmpsaccs.C_presaccfixslope)...
	,bins)); 
 title('Transsaccadic slope');

 figure(5)
 clf
 plot(saccs.S_em);
 hold on
 plot(saccs.S_times, saccs.S_em(saccs.S_times),'k.');
 plot(saccs.S_times, -3, 'r*');
 text(saccs.S_times, saccs.S_em(saccs.S_times)-0.3,num2cell(saccs.S_stop-saccs.S_start),'FontSize',6);
 plot(saccs.S_times(~saccs.S_corrective), 0.2, 'r+');
 if sum(saccs.S_corrective)
    plot(saccs.S_times(saccs.S_corrective), 0.4, 'm+');
 end

 figure(6)
 clf
 for i=1:length(saccs.S_times)
   form = saccs.S_forms(:,i);
   start = saccs.S_times(i)-FORM_MS_HWIN;

   stop  = saccs.S_times(i)+FORM_MS_HWIN;
   plot([start:stop],saccs.S_em([start:stop]),'b'); 
   hold on
 end

end


    


