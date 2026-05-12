function [resp,stim,timestim,h,time] = ...
      convolBox(file, tr, stimdef, cutindex, h_choice)
%
%% calculates convolution of stimulus with GAMMA function
%%
%%   file: output [! NOT implemented !]
%%   tr: repetition time [sec]
%%   stimdef = [40 1 40 48 2]; 
%%	     (controlstart,controlend,[stimon,stimoff],nrepeat)
%%   cutindex = [5 5]
%%   h_choice = 1: Friston 1994 Poisson model
%%              2: Cohen 1997 -- Gamma variate impulse model
%%
%%   June 2000, JP
%---------------------------------------------

if ( h_choice < 0 )
   f_verbose = 0;
   h_choice = -1*h_choice;
else
   f_verbose = 1;
end

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
   h = gam_a .* (time.^gam_b).* exp(-time/gam_c);
else
   error(sprintf('convol_box(): undefined h_choice [1..2] = %d',h_choice))
end
%plot(time,h)

   % --- create stimuli%
n_stim = sum(stimdef(1,1:2))+((sum(stimdef(1,3:4))*stimdef(1,5)));
timestim = (0:n_stim-1)*tr;
stim = zeros(1,n_stim);
spos = stimdef(1,1);
for istim=0:stimdef(1,5)-1
   stim(1,[spos+1:spos+stimdef(1,3)]) = 1.;
   spos = spos+stimdef(1,3)+stimdef(1,4);
end

   % ---convolution
resp = conv(stim,h);
resp = resp( (cutindex(1,1)+1) : (n_stim-cutindex(1,2)) );

if ( f_verbose )
   figure(1)
   plot(timestim, stim, '*')
   title('time')
   hold on
   plot(time, h)
   hold off

   figure(2)
   plot(timestim,resp) 
   title('RESPONSE')
   disp( sprintf('--- number of time points = %d', n_stim) )
end


%--------------------------------------------
% stimdef=[40 1 40 48 6]
% controlstart, controlend, [stimon, stimoff], nrepeat%

% convol_box('sli_x', 1,  [40 1 40 48 2], [0 0], 1);
% convol_box('sli_x', 0.11,  [0 0 1 119 1], [0 0], 2);

