function sacc2d = saccadeOKNAnalysis2D(hsacc,vsacc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  FUNCTION:
%    saccadeAnalysis2D(hsacc,vsacc)
%
%  PURPOSE: 
%    This function takes two saccade lists (the horizontal and vertical
%    lists, respectively) generated in the saccadeAnalysis1D function and
%    finds matches between them.  It works by first identifying candidate 
%    matches, i.e. those having less than MAX_HV_MATCH_DELAY between saccades
%    detected in the two traces.  It then selects the saccade that has the 
%    largest amplitude from the possible candidates, and uses the mean 
%    peak velocity time between the two traces to identify the saccade 
%    midpoint.  Finally, the fixation periods are calculated by performing a 
%    a linear regression between successive saccades in H and V.
%
%   D.Leopold  15-APR-00
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


H_num_Vmatches = 0;
V_num_Hmatches = 0;

MAX_HV_MATCH_DELAY = 10;    %ms


if length(hsacc.S_times) < 2 | length(vsacc.S_times) < 2
   sacc2d = [];
   return;
end

% first we look for probably saccade matches between the two eyes
% within a particular acceptance window
for i=1:length(hsacc.S_times)
   %compare each saccade in the vertical to each in the horizontal  
   vh = vsacc.S_times-hsacc.S_times(i);
   VHdiff(:,i) = vh';
end

for i=1:length(vsacc.S_times)
   %compare each saccade int the vertical to each in the horizontal  
   hv = hsacc.S_times-vsacc.S_times(i);	
   HVdiff(:,i) = hv';
end


% Here we identify probable saccade matches with binary lists.  
% Because the identification in the individual lists is 
% some errors, there is a MAX_HV_MATCH_DELAY acceptance window.
%
Hboth_match = (abs(HVdiff') < MAX_HV_MATCH_DELAY);
Vboth_match = (abs(VHdiff') < MAX_HV_MATCH_DELAY);
    
	    
% Then we get the number of matches for each saccade.  These 
% lists shouldn't be all that different from one another.
%
if ~isempty(Hboth_match)
   H_num_Vmatches = sum(Hboth_match);
end

if ~isempty(Vboth_match)
   V_num_Hmatches = sum(Vboth_match);
end

H_val2Dsacc = H_num_Vmatches > 0;
V_val2Dsacc = V_num_Hmatches > 0;


H_val2Dsaccindx = find(H_val2Dsacc);
V_val2Dsaccindx = find(V_val2Dsacc);

multiple_matches = sum(H_num_Vmatches > 1) + sum(V_num_Hmatches > 1);
if multiple_matches
   % Get rid of multiple matches.
   h_multi = find(H_num_Vmatches > 1);
   v_multi = find(V_num_Hmatches > 1);

   % must get rid of all except one	
	
   for i=1:length(h_multi)
      h_indx           =  h_multi(i)				    % index of a problematic horiz saccade
      h_allvals        =  Hboth_match(:,h_indx);		    % vals for which there are multiple vert matches 
      v_matchindxs     =  find(h_allvals);   			    % indices for accessing the vert saccade list
      v_amps 	       =  vsacc.S_amps(v_matchindxs);		    % decide based on amplitude for now
      v_biggest	       =  v_amps == max(abs(v_amps));		    % the biggest vert sacc wins
      v_toosmall       =  v_matchindxs(find(~v_biggest));	    % the smallest ones lose
      V_val2Dsaccindx  =  setxor(V_val2Dsaccindx,v_toosmall);       % the smallest times are eliminated
   end	   

   for i=1:length(v_multi)
      v_indx           =  v_multi(i)				    % index of a problematic vertical saccade
      v_allvals        =  Vboth_match(:,v_indx);		    % vals for which there are multiple horiz matches 
      h_matchindxs     =  find(v_allvals);  			    % indices for accessing the horizontal saccade list
      h_amps 	       =  hsacc.S_amps(h_matchindxs);		    % decide based on amplitude for now
      h_biggest	       =  h_amps == max(abs(h_amps));		    % the biggest horizontal sacc wins
      h_toosmall       =  h_matchindxs(find(~h_biggest));	    % the smallest ones lose
      H_val2Dsaccindx  =  setxor(H_val2Dsaccindx,h_toosmall);       % the smallest times are eliminated
   end	   
end


Hmatch_times     =  hsacc.S_times(H_val2Dsaccindx);
Vmatch_times     =  vsacc.S_times(V_val2Dsaccindx);

%
% Last resort to get rid of multiple matches...
%
if length(Hmatch_times) ~= length(Vmatch_times) 
  if length(Hmatch_times) > length(Vmatch_times)
    htmp = [];
    vtmp = [];
    hi = [];
    vi = [];
    for i=1:length(Vmatch_times)
       % if nothing else works...
       t = Vmatch_times(i);
       [minval, minindx] = min(abs(t-Hmatch_times));
       htmp(i) = Hmatch_times(minindx); 
       hi(i)   = H_val2Dsaccindx(minindx);	
    end
    Hmatch_times = htmp;	
    H_val2Dsaccindx = hi;
  else	 
    for i=1:length(Hmatch_times)
       % if nothing else works...
       t = Hmatch_times(i);
       [minval, minindx] = min(abs(t-Vmatch_times));
       vtmp(i) = Vmatch_times(minindx); 
       vi(i)   = H_val2Dsaccindx(minindx);	
    end
    Vmatch_times = vtmp;
    V_val2Dsaccindx = vi;
 end
end
match_times      =  (Hmatch_times+Vmatch_times)/2.0;

% Found match... collect info about matching saccades
match_start  = min(hsacc.S_start(H_val2Dsaccindx'),vsacc.S_start(V_val2Dsaccindx'));
match_stop   = max(hsacc.S_stop(H_val2Dsaccindx'),vsacc.S_stop(V_val2Dsaccindx'));
	    

% Create the new saccade lists.  Start with the horizontal list (which contains the horizontal only
% information, and add the 'vertical only's to the end.  This is, in essence an 'or' operation.

tmptimes = hsacc.S_times;
tmptimes(H_val2Dsaccindx) = match_times;
new_S_times = [tmptimes vsacc.S_times(~V_val2Dsacc)];

tmpstart = hsacc.S_start;
tmpstart(H_val2Dsaccindx) = match_start;
new_S_start = [tmpstart vsacc.S_start(~V_val2Dsacc)];

tmpstop = hsacc.S_stop;
tmpstop(H_val2Dsaccindx) = match_stop;
new_S_stop = [tmpstop vsacc.S_stop(~V_val2Dsacc)];

tmpcorrective    	= hsacc.S_corrective;
new_S_corrective 	= [tmpcorrective vsacc.S_corrective(~V_val2Dsacc)];

[new_S_times, srt]  	= sort(new_S_times);
new_S_times 		= round(new_S_times);
new_S_start         	= round(new_S_start(srt));
new_S_stop          	= round(new_S_stop(srt));
Hstart_pos 		= hsacc.S_em(new_S_start);
Hstop_pos  		= hsacc.S_em(new_S_stop);
Vstart_pos 		= vsacc.S_em(new_S_start);
Vstop_pos  		= vsacc.S_em(new_S_stop);
Htime_pos  		= hsacc.S_em(new_S_times);	
Vtime_pos  		= vsacc.S_em(new_S_times);	
Hdiffs     		= Hstart_pos-Hstop_pos;
Vdiffs     		= Vstart_pos-Vstop_pos;
S_2damp	   		= sqrt((Hdiffs.^2)+(Vdiffs.^2));
Hpeak_vels 		= hsacc.S_emvel(new_S_times);	
Vpeak_vels 		= vsacc.S_emvel(new_S_times);	
S_2dvel    		= sqrt((Hpeak_vels.^2)+(Vpeak_vels.^2));	

sacc2d.S_times     	= new_S_times;
sacc2d.S_start     	= new_S_start;
sacc2d.S_stop      	= new_S_stop;
sacc2d.S_duration      	= new_S_stop-new_S_start;
sacc2d.S_startHpos 	= Hstart_pos;
sacc2d.S_stopHpos  	= Hstop_pos;
sacc2d.S_timeHpos  	= Htime_pos;
sacc2d.S_startVpos 	= Vstart_pos;
sacc2d.S_stopVpos  	= Vstop_pos;
sacc2d.S_timeVpos  	= Vtime_pos;
sacc2d.S_Hamps	    	= abs(Hdiffs);
sacc2d.S_Vamps	    	= abs(Vdiffs);
sacc2d.S_2damps	    	= S_2damp;
sacc2d.S_Hpeak_vels 	= Hpeak_vels;
sacc2d.S_Vpeak_vels 	= Vpeak_vels;
sacc2d.S_2dpeakvels 	= S_2dvel;
sacc2d.S_corrective    	= new_S_corrective;
sacc2d.S_blinktimes    	= vsacc.S_blinktimes;
sacc2d.S_blinkamps    	= vsacc.S_blinkamps;
sacc2d.S_blinkstart    	= vsacc.S_blinkstart;
sacc2d.S_blinkstop    	= vsacc.S_blinkstop;


%
% Calculate the slopes/y intercepts of the 'intersaccade' periods.
% 
totfix = 0;
for i = 1:(length(sacc2d.S_times)-1)
  % first check to see if there is a blink in the period
  beginper = sacc2d.S_stop(i);
  endper   = sacc2d.S_start(i+1);
  blinkstart = 0;
  for j=1:length(sacc2d.S_blinktimes)
     if (sacc2d.S_blinktimes(j) > beginper) & (sacc2d.S_blinktimes(j) < endper)
        blinkstart = sacc2d.S_blinkstart(j);
        blinkstop = sacc2d.S_blinkstop(j);
     end     
  end
  if blinkstart
    % if there is a blink within the fixation period...
    endper2 = endper;
    endper = blinkstart;
    beginper2 = blinkstop;
  end
  totfix = totfix+1;
  t = [beginper:endper];
  len = length(t); 
  [hp,hs] = polyfit(t,hsacc.S_em(t),1);   
  [vp,vs] = polyfit(t,vsacc.S_em(t),1);   
  allhslope(totfix)     = hp(1);
  allvslope(totfix)     = vp(1);
  allhintercept(totfix) = hp(2);	
  allvintercept(totfix) = vp(2);	
  allhresid(totfix)     = hs.normr;
  allvresid(totfix)     = vs.normr;
  allfixstart(totfix)   = beginper;
  allfixstop(totfix)    = endper;

  if blinkstart
     % ...calculate two fixation periods, one before the blink and one after.
     totfix = totfix+1;
     t = [beginper2:endper2];
     len = length(t); 
     [hp,hs] = polyfit(t,hsacc.S_em(t),1);   
     [vp,vs] = polyfit(t,vsacc.S_em(t),1);   
     allhslope(totfix)     = hp(1);
     allvslope(totfix)     = vp(1);
     allhintercept(totfix) = hp(2);	
     allvintercept(totfix) = vp(2);	
     allhresid(totfix)     = hs.normr;
     allvresid(totfix)     = vs.normr;
     allfixstart(totfix)   = beginper2;
     allfixstop(totfix)    = endper2;
  end

end

sacc2d.S_Hfixslope      = allhslope;
sacc2d.S_Vfixslope      = allvslope;
sacc2d.S_Hfixintercept  = allhintercept;
sacc2d.S_Vfixintercept  = allvintercept;
sacc2d.S_Hfixresid  	= allhresid;
sacc2d.S_Vfixresid  	= allvresid;
sacc2d.S_fixstart  	    = allfixstart;
sacc2d.S_fixstop  	    = allfixstop;
sacc2d.S_driftspeed	    = 1000.0*sqrt(allhslope.^2+allvslope.^2);         % drifting speed in deg/sec
sacc2d.S_driftdir	= 90-(180./3.14159265)*atan2(allvslope,allhslope);			 			

	

