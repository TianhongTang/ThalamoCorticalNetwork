function [returnResp,returnStim,returnImpulse] = ...
    erLinModel(file, type, cutindex, h_choice)

if nargin > 0
    [returnResp,returnStim,returnImpulse] = ...
    erLinModelFun(file, type, cutindex, h_choice);
else
% type = 
% controlstart, [stimon, stimoff], nrepeat, controlend, TR [s] %

%erLinearModel('test', [0 2*4 30*4 1 0 0.25], [0 0],3);

%erLinearModel('stim2_10_0.25', [0 2*4 10*4 1 0 0.25], [8 4],2);
%erLinearModel('stim2_30_0.5', [0 2*2 30*2 1 0 0.5], [4 0],2);
%erLinearModel('stim2_10_0.5', [0 2*2 10*2 1 0 0.5], [4 0],2);

%erLinearModel('stim4_11_0.25N8', [0 4*4 11*4 8 0 0.25], [0 4],2);
%erLinearModel('stim4_11_0.25N1', [0 4*4 11*4 1 0 0.25], [8 12],2);
%erLinearModel('stim6_9_0.25N8',  [0 6*4  9*4 8 0 0.25], [0 4],2);
%erLinearModel('stim12_18_2N8',  [0 12/2 18/2 8 0 2.0], [0 2],2);
%erLinearModel('stim12_18_4N8',  [0 12/4 18/4 8 0 4.0], [0 1],2);

%erLinearModel('stim4.8_9.6_0.6N8',  [0 8 16 8 0 0.6], [0 3],2);
%erLinearModel('stim4.8_9.6_0.6N4',  [0 8 16 4 0 0.6], [0 3],2);
%erLinearModel('stim4.8_9.6_2.4N8',  [0 2 4 8 0 2.4], [0 1],2);
%erLinearModel('stim4.0_8.0_0.5N8',  [0 8 16 8 0 0.5], [0 2],2);
%erLinearModel('stim4.0_8.0_0.5N1',  [0 8 16 1 0 0.5], [2 0],2);
%erLinearModel('stim4.0_8.0_1.0N8',  [0 4 8 8 0 1.0], [0 1],2);

%erLinearModel('stim8_8_6.0N4',  [0 8 8 4 0 6.0], [0 0],2);  % EPI13
%erLinearModel('stim8_8_6.0N2',  [0 16 16 2 0 6.0], [0 0],2);  % EPI13_HSE

%erLinearModel('stim3_10_1.0N28',  [3 3 10 28 0 1.0], [0 3],2);  % Biomotion
%erLinearModel('stim5_5_6.4N11',  [5 5 5 11 0 6.444], [0 0],2);  % Stelios

%erLinearModel('stim4_8_1.0N4',  [0 4 8 4 2 1.0], [0 0],2);  % George C00

end
return
%---------------------------------------------

function [returnResp,returnStim,returnImpulse] = ...
    erLinModelFun(file, type, cutindex, h_choice)
%
%% calculates convolution of stimulus with GAMMA function
%%
%%   file: output [! NOT implemented !]
%%   type = [40 40 48 4 1 0.25];    "number of images"
%%	     (controlstart,[stimon,stimoff],nrepeat,controlend, tr[sec])
%%   cutindex = [5 5]
%%   h_choice = 1: Friston 1994 Poisson model
%%              2: Cohen 1997 -- Gamma variate impulse model
%%
%%   June 2000, JP
%---------------------------------------------

outputDir = '//wks5/guest/mpidata/stim_txc/';
savename = file;

if ( h_choice < 0 )
   f_verbose = 0;
   h_choice = -1*h_choice;
else
   f_verbose = 1;
end
typeCont1 = type(1);
typeStimOn = type(2);
typeStimOff = type(3);
typeStimRep = type(4);
typeCont2 = type(5);
typeTR = type(6);

tr = typeTR;

   % ---  MODEL 1: Friston 1994 Poisson model
if ( h_choice == 1 )
   if mod(tr,1) ~= 0
      sprintf('--- Warning: tr not INT!!!')
      tr = round(tr)
   end

   n_h = 20./tr;
   lam = 7.37;
   time = (0:floor(n_h)-1)*tr;

   for itime= 1:n_h
      h1 =1.00;
      for i = 1:time(itime)
         h1 = h1*(lam/i);
      end 
      h(itime)=h1*exp(-lam);
   end

   % --- MODEL 2: Cohen 1997 -- 
   %     Gamma variate impulse model G(t, gama_a, gam_b, gam_c)
elseif ( h_choice == 2 )
   gam_a = 0.008;
   gam_b = 8.60;
   gam_c = 0.547;
   if ( f_verbose )
      disp( sprintf('Gamma model: %f %f',gam_b, gam_c));
      disp( sprintf('Gamma model: time-to-peak [s]= %f', gam_b*gam_c) );
      disp( sprintf('Gamma model: FWHM [s]~ %f', 2.35 * sqrt(gam_b) * gam_c) );
   end

   n_h = 20./tr;
   time = ( 0:floor(n_h)-1 ) * tr;
%%   h = gam_a .* (time.^gam_b).* exp(-time/gam_c);
   h = funGamma([gam_a gam_b gam_c],time);
   
   % --- MODEL 3: Friston SPM -- 
   %     Gamma variate impulse model with undershoot
elseif ( h_choice == 3 )
   gam_a = 0.008;
   gam_b = 8.60;
   gam_c = 0.547;
   if ( f_verbose )
      disp( sprintf('Gamma model: %f %f',gam_b, gam_c));
      disp( sprintf('Gamma model: time-to-peak [s]= %f', gam_b*gam_c) );
      disp( sprintf('Gamma model: FWHM [s]~ %f', 2.35 * sqrt(gam_b) * gam_c) );
   end

   n_h = 60./tr;
   time = ( 0:floor(n_h)-1 ) * tr;
   h = funGammaSpm(time);
else
   error(sprintf('convol_box(): undefined h_choice [1..2] = %d',h_choice))
end
%plot(time,h)

typeCont1 = type(1);
typeStimOn = type(2);
typeStimOff = type(3);
typeStimRep = type(4);
typeCont2 = type(5);

   % --- create stimuli%
n_stim = typeCont1 + (typeStimOn+typeStimOff)*typeStimRep + typeCont2;
timestim = (0:n_stim-1)*tr;
stim = zeros(1,n_stim);
spos = typeCont1;
for istim=0:typeStimRep-1
   stim(1,[spos+1:spos+typeStimOn]) = 1.;
   spos = spos+typeStimOn+typeStimOff;
end

   % ---convolution
resp = conv(stim,h);
resp = resp( (cutindex(1,1)+1) : (n_stim-cutindex(1,2)) );

returnResp = zeros(length(resp),2);
returnResp(:,1) = timestim((cutindex(1,1)+1) : (n_stim-cutindex(1,2)))';
returnResp(:,2) = resp';

returnStim = zeros(length(stim),2);
returnStim(:,1) = timestim';
returnStim(:,2) = stim';

returnImpulse = zeros(length(h),2);
returnImpulse(:,1) = time';
returnImpulse(:,2) = h';


if ( f_verbose )
   figure(1)
   plot(returnStim(:,1), returnStim(:,2), '*')
   title('time')
   hold on
   plot(returnImpulse(:,1), returnImpulse(:,2))
   hold off

   figure(2)
   plot(returnResp(:,1), returnResp(:,2)) 
   title('RESPONSE')
   disp( sprintf('--- number of RESPONSE time points = %d', length(returnResp(:,2))) )
end

key('--- save data: ')
data2write = resp';
save(strcat(outputDir, file, '.asc'), 'data2write', '-ascii')
fprintf('\n<%s> written\n', strcat(outputDir, file, '.asc'));

return
%--------------------------------------------
