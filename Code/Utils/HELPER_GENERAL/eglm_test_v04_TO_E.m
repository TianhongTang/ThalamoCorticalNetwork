%addpath('C:\Users\Takaya\Dropbox\TAKAYA_paper\SpikesAnalyses_Workflow\dependency\HELPER\HELPER_GENERAL');
%addpath('C:\Users\Takaya\Dropbox\TAKAYA_paper\SpikesAnalyses_Workflow\dependency\HELPER\HELPER_GENERAL\eglm_20220706');
clear
close all

ufractals = [6101 6102 6103 6104 6105 6201 6202 6203 6204 6205 6301 6302 6303 6304 6305 6401 6402 6403 6404 6405]';

ufractals_x = [1:20];

ufractals_er = [1 0.75 0.5 0.25 0 0 0.25 0.5 0.75 1 0.5 0.5 0.5 0.5 0.5 1 0.75 0.5 0.25 0];

ufractals_en = [1 0.75 0.5 0.25 0 1 0.75 0.5 0.25 0 1 0.75 0.5 0.25 0 0 0 0 0 0];

ufractals_sdr = [0 .433 .5 .433 0 0 .433 .5 .433 0 0 0 0 0 0 0 .433 .5 .433 0];
% sqrt(0.25*(1 - 0.25)^2 + 0.75*(-0.25)^2)
ufractals_sdn = [0 .433 .5 .433 0 0 .433 .5 .433 0 0 .433 .5 .433 0 0 0 0 0 0];

ufractalname = {'6101','6102','6103','6104','6105','6201','6202','6203','6204','6205' '6301','6302','6303','6304','6305','6401','6402','6403','6404','6405'};

%Dir = 'C:\Users\Takaya\Dropbox\RewPredNov\Workflow\dependency\';
Dir = '/Users/esb/Documents/work/proj/takayanovelty/';
load([Dir 'savedatabehavior_choice.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choices -- Column#1 - Offer1, #2 - Offer2, #3 - Chosen Offer, #4 - MonkeyID, #5 - Session types  

% MonkeyID -- 1 = Slayer, 2 = Sabbath, 3 = Zeppelin 

% Session types
%1 - Control
%2 - Control injection (DCZ)
%3 - Control injection (DMSO)
%4 - 1day post control injection
%6 - DCZ
%7 - 1day post DCZ
%8 - 3days post DCZ
%9 - Failed DCZ
%10 - all choice behavior without injections

% LateralityofCue -- Column#1 - Offer1, #2 - Offer2
% assign 1 if the offer is on the left, assign 0 if it's on the right, otherwise assign NaN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% If want to exclude block 3and4 %%%
% ch=choices;
% clear choices
% noblock3_1=find(ch(:,1)<6300 );
% noblock3_2=find(ch(:,2)<6300 );
% noblock34=intersect(noblock3_1,noblock3_2);
% ch=ch(noblock34,:);
% choices=ch;
% 
% LC = LateralityofCue;
% clear LateralityofCue
% LC = LC(noblock34, :);
% LateralityofCue = LC;

monkchoices = {};
monkchoices{end+1} =(choices(:,1:3));

xy = vertcat(monkchoices{:});

% convert original fractal ID#s in file to sequential fractal indexes into
% our lists of their properties (e.g. ufractals_er, _sdr, ...)
xy_new = nan(size(xy));
for ufi = 1:numel(ufractals)
    ok = xy == ufractals(ufi);
    xy_new(ok) = ufi;
end
assert(~any(isnan(xy_new(:))),'could not convert from original to new fractal ID#s!');
xy = xy_new;

addinfo = choices(:, 5:6);

% cols 1,2,3,4,5 =  offer1, offer2, chosen offer, monk ID, session types
xym = [xy addinfo];

model_type = 'per_fractal';
model_type = 'key factors';

plotmode = 'line';
plotmode = 'bar';

rescale_parameters = {};

switch model_type
    case 'per_fractal'
        xreg = {};
        plotparams = {};
        for ufi = 1:numel(ufractals)
            if ~strcmp(ufractalname(ufi),'6405')
                xreg{end+1} = struct('name',ufractalname{ufi},'terms',{{xym(:,1:2) == ufi}},'normalization','none');
                if ufractalname{ufi}(1) == 'I'
                    curcolor = 'r';
                else
                    curcolor = 'b';
                end
                plotparams{end+1} = {ufractalname{ufi}, 'color', curcolor, 'x',ufractals_x(ufi)};
            end
        end
        xreg{end+1} = struct('name','rgt','terms',{{[zeros(size(xym,1),1) ones(size(xym,1),1)]}});
        plotparams{end+1} = {'rgt'};
        
        %rescale_parameters = {'rescale betas',struct('name','big juice','xreg_name','NR100','xreg_change',1./0.4)};
    case 'key factors'
        xreg = {};
        plotparams = {};

%       xreg{end+1} = struct('name','ER','terms',{{ufractals_er(xym(:,1:2))}});
%       xreg{end+1} = struct('name','SDR','terms',{{ufractals_sdr(xym(:,1:2))}});
%       xreg{end+1} = struct('name','Nv','terms',{{ufractals_en(xym(:,1:2))}});
%       xreg{end+1} = struct('name','SDNv','terms',{{ufractals_sdn(xym(:,1:2))}});
%       xreg{end+1} = struct('name','ERxNv','terms',{{ ufractals_er(xym(:,1:2)), ufractals_en(xym(:,1:2)) }});
%       xreg{end+1} = struct('name','Lt','terms',{{LateralityofCue}}); % Laterality

% OLD
%       xreg{end+1} = struct('name','ERxLt','terms',{{ufractals_er(xym(:,1:2)), LateralityofCue}});
%       xreg{end+1} = struct('name','SDRxLt','terms',{{ufractals_sdr(xym(:,1:2)), LateralityofCue}});
%       xreg{end+1} = struct('name','NvxLt','terms',{{ufractals_en(xym(:,1:2)), LateralityofCue}});
%       xreg{end+1} = struct('name','ERxNvxLt','terms',{{ufractals_er(xym(:,1:2)), ufractals_en(xym(:,1:2)), LateralityofCue}});

      testmode = 5;
      switch testmode
          case 1
              % basic model
              xreg{end+1} = struct('name','ER','terms',{{ufractals_er(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDR','terms',{{ufractals_sdr(xym(:,1:2))}});
              xreg{end+1} = struct('name','Nv','terms',{{ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDNv','terms',{{ufractals_sdn(xym(:,1:2))}});

              xreg{end+1} = struct('name','ERxSDR','terms',{{ufractals_er(xym(:,1:2)),ufractals_sdr(xym(:,1:2))}});
              xreg{end+1} = struct('name','ERxNv','terms',{{ufractals_er(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','ERxSDNv','terms',{{ufractals_er(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDRxNv','terms',{{ufractals_sdr(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDRxSDNv','terms',{{ufractals_sdr(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});
              xreg{end+1} = struct('name','NvxSDNv','terms',{{ufractals_en(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});

              xreg{end+1} = struct('name','Lt','terms',{{LateralityofCue}}); % Laterality
          case 2
              % with laterality interactions (full terms)
              xreg{end+1} = struct('name','ER','terms',{{ufractals_er(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDR','terms',{{ufractals_sdr(xym(:,1:2))}});
              xreg{end+1} = struct('name','Nv','terms',{{ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDNv','terms',{{ufractals_sdn(xym(:,1:2))}});

              xreg{end+1} = struct('name','ERxSDR','terms',{{ufractals_er(xym(:,1:2)),ufractals_sdr(xym(:,1:2))}});
              xreg{end+1} = struct('name','ERxNv','terms',{{ufractals_er(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','ERxSDNv','terms',{{ufractals_er(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDRxNv','terms',{{ufractals_sdr(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDRxSDNv','terms',{{ufractals_sdr(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});
              xreg{end+1} = struct('name','NvxSDNv','terms',{{ufractals_en(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});

              xreg{end+1} = struct('name','Lt','terms',{{LateralityofCue}}); % Laterality
              
              xreg{end+1} = struct('name','LtxER','terms',{{LateralityofCue,ufractals_er(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxSDR','terms',{{LateralityofCue,ufractals_sdr(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxNv','terms',{{LateralityofCue,ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxSDNv','terms',{{LateralityofCue,ufractals_sdn(xym(:,1:2))}});

              xreg{end+1} = struct('name','LtxERxSDR','terms',{{LateralityofCue,ufractals_er(xym(:,1:2)),ufractals_sdr(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxERxNv','terms',{{LateralityofCue,ufractals_er(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxERxSDNv','terms',{{LateralityofCue,ufractals_er(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxSDRxNv','terms',{{LateralityofCue,ufractals_sdr(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxSDRxSDNv','terms',{{LateralityofCue,ufractals_sdr(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxNvxSDNv','terms',{{LateralityofCue,ufractals_en(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});
          case 3
              
              % with laterality interactions (No SDNv)
              xreg{end+1} = struct('name','ER','terms',{{ufractals_er(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDR','terms',{{ufractals_sdr(xym(:,1:2))}});
              xreg{end+1} = struct('name','Nv','terms',{{ufractals_en(xym(:,1:2))}});
              
              xreg{end+1} = struct('name','ERxSDR','terms',{{ufractals_er(xym(:,1:2)),ufractals_sdr(xym(:,1:2))}});
              xreg{end+1} = struct('name','ERxNv','terms',{{ufractals_er(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDRxNv','terms',{{ufractals_sdr(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              
              xreg{end+1} = struct('name','Lt','terms',{{LateralityofCue}}); % Laterality
              
              xreg{end+1} = struct('name','LtxER','terms',{{LateralityofCue,ufractals_er(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxSDR','terms',{{LateralityofCue,ufractals_sdr(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxNv','terms',{{LateralityofCue,ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxERxSDR','terms',{{LateralityofCue,ufractals_er(xym(:,1:2)),ufractals_sdr(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxERxNv','terms',{{LateralityofCue,ufractals_er(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxSDRxNv','terms',{{LateralityofCue,ufractals_sdr(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              
          case 4
              % with laterality interactions, and no ErxSDR, NxSDN, LtxErxSDR, LtxNxSDN
              xreg{end+1} = struct('name','ER','terms',{{ufractals_er(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDR','terms',{{ufractals_sdr(xym(:,1:2))}});
              xreg{end+1} = struct('name','Nv','terms',{{ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDNv','terms',{{ufractals_sdn(xym(:,1:2))}});

              xreg{end+1} = struct('name','ERxNv','terms',{{ufractals_er(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','ERxSDNv','terms',{{ufractals_er(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDRxNv','terms',{{ufractals_sdr(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDRxSDNv','terms',{{ufractals_sdr(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});

              xreg{end+1} = struct('name','Lt','terms',{{LateralityofCue}}); % Laterality
              
              xreg{end+1} = struct('name','LtxER','terms',{{LateralityofCue,ufractals_er(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxSDR','terms',{{LateralityofCue,ufractals_sdr(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxNv','terms',{{LateralityofCue,ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxSDNv','terms',{{LateralityofCue,ufractals_sdn(xym(:,1:2))}});

              xreg{end+1} = struct('name','LtxERxNv','terms',{{LateralityofCue,ufractals_er(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxERxSDNv','terms',{{LateralityofCue,ufractals_er(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxSDRxNv','terms',{{LateralityofCue,ufractals_sdr(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','LtxSDRxSDNv','terms',{{LateralityofCue,ufractals_sdr(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});

          case 5
              % without laterality, and no ErxSDR, no NxSDN
              xreg{end+1} = struct('name','ER','terms',{{ufractals_er(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDR','terms',{{ufractals_sdr(xym(:,1:2))}});
              xreg{end+1} = struct('name','Nv','terms',{{ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDNv','terms',{{ufractals_sdn(xym(:,1:2))}});

              xreg{end+1} = struct('name','ERxNv','terms',{{ufractals_er(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','ERxSDNv','terms',{{ufractals_er(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDRxNv','terms',{{ufractals_sdr(xym(:,1:2)),ufractals_en(xym(:,1:2))}});
              xreg{end+1} = struct('name','SDRxSDNv','terms',{{ufractals_sdr(xym(:,1:2)),ufractals_sdn(xym(:,1:2))}});
          otherwise
              error('unknown case');
      end
      

    otherwise
        error('unknown model');
end


yreg  = struct('name','chose Offer2','value',xym(:,3) == xym(:,2)); % chose offer2


% sessions type
%1 Control
%2 Control injection (DCZ)
%3 Control injection (DMSO)
%4 1day post control injection
%5 After surgery
%6 DCZ
%7 1day post DCZ
%8 3days post DCZ
%9 Failed DCZ
%10 previous all choice data (id =1, Slayer; id=2, Sabbath)



% Control
fit_Control = eglm_fit(xreg,yreg,'binary choices','oktrials',xym(:,4) == 1 & xym(:, 5) == 1,'name','Control',rescale_parameters{:});

% Control inejction (DCZ)
fit_ControlDCZ = eglm_fit(xreg,yreg,'binary choices','oktrials',xym(:,4) == 1 & xym(:, 5) == 2,'name','ControlDCZ',rescale_parameters{:});

% Control injection (DMSO)
fit_ControlDMSO = eglm_fit(xreg,yreg,'binary choices','oktrials',xym(:,4) == 1 & xym(:, 5) == 3,'name','ControlDMSO',rescale_parameters{:});

% Control injection (DCZ + DMSO)
clear temp temp_
temp = xym(:, 5) == 2 | xym(:, 5) == 3;
temp_ = xym(:,4) == 1 & temp;
fit_ControlInjection = eglm_fit(xreg,yreg,'binary choices','oktrials',temp_,'name','ControlDMSO',rescale_parameters{:});

% 1 day post control injection
fit_PostControl = eglm_fit(xreg,yreg,'binary choices','oktrials',xym(:,4) == 1 & xym(:, 5) == 4,'name','PostControl',rescale_parameters{:});

% DCZ
fit_DCZ = eglm_fit(xreg,yreg,'binary choices','oktrials',xym(:,4) == 1 & xym(:, 5) == 6,'name','DCZ',rescale_parameters{:});
% 1 day post DCZ
fit_1PostDCZ = eglm_fit(xreg,yreg,'binary choices','oktrials',xym(:,4) == 1 & xym(:, 5) == 7,'name','PostDCZ',rescale_parameters{:});
% 3days post DCZ
fit_3PostDCZ = eglm_fit(xreg,yreg,'binary choices','oktrials',xym(:,4) == 1 & xym(:, 5) == 8,'name','3PostDCZ',rescale_parameters{:});
% failed DCZ injection
fit_failedDCZ = eglm_fit(xreg,yreg,'binary choices','oktrials',xym(:,4) == 1 & xym(:, 5) == 9,'name','failedDCZ',rescale_parameters{:});

% all choices except sessions with injections
clear temp temp_
temp = xym(:, 5) == 1 | xym(:, 5) == 10;
temp_ = xym(:,4) == 1 & temp;
fit_ALlchoiceSl = eglm_fit(xreg,yreg,'binary choices','oktrials',temp_,'name','all choice Slayer',rescale_parameters{:});
fit_ALlchoiceSb = eglm_fit(xreg,yreg,'binary choices','oktrials',xym(:,4) == 2 & xym(:, 5) == 10,'name','all choice Sabbath',rescale_parameters{:});
fit_ALlchoiceZP = eglm_fit(xreg,yreg,'binary choices','oktrials',xym(:,4) == 3 & xym(:, 5) == 10,'name','all choice Zeppelin',rescale_parameters{:});
%% All choices without injections for each monkey
figuren;

nsubplot(3,1,1,1);
xname = fit_ALlchoiceSl.xname;
plotparams_fun = @(names,color) cellfun(@(z) {z,'color',color},names,'uniform',0);
plotparams = plotparams_fun(xname,'k');
eglm_plot_fit(fit_ALlchoiceSl,'params',plotparams,'plotmode',plotmode);
plot([-1 1]*10000, [0 0], '-k','xliminclude','off');
xtickangle(45);

nsubplot(3,1,2,1);
xname = fit_ALlchoiceSb.xname;
plotparams_fun = @(names,color) cellfun(@(z) {z,'color',color},names,'uniform',0);
plotparams = plotparams_fun(xname,'k');
eglm_plot_fit(fit_ALlchoiceSb,'params',plotparams,'plotmode',plotmode);
plot([-1 1]*10000, [0 0], '-k','xliminclude','off');
xtickangle(45);

nsubplot(3,1,3,1);
xname = fit_ALlchoiceZP.xname;
plotparams_fun = @(names,color) cellfun(@(z) {z,'color',color},names,'uniform',0);
plotparams = plotparams_fun(xname,'k');
eglm_plot_fit(fit_ALlchoiceZP,'params',plotparams,'plotmode',plotmode);
plot([-1 1]*10000, [0 0], '-k','xliminclude','off');
xtickangle(45);

% set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figure and orientation
% 
% destination = 'C:\Users\Takaya\Dropbox\RewPredNov\Workflow\Figure\';
% filename = 'GLMfit';
% print([destination filename], '-djpeg');



%% compare DCZ and control
figuren;

xname = fit_Control.xname;
plotparams_fun = @(names,color) cellfun(@(z) {z,'color',color},names,'uniform',0);

plotparams = plotparams_fun(xname,'k');
eglm_plot_fit(fit_Control,'params',plotparams,'plotmode',plotmode);
hold on;
plotparams = plotparams_fun(xname,[1 1 1]*.5);
eglm_plot_fit(fit_ControlInjection,'params',plotparams,'plotmode',plotmode);
plotparams = plotparams_fun(xname,'r');
eglm_plot_fit(fit_DCZ,'params',plotparams,'plotmode',plotmode);
% plotparams = plotparams_fun(xname,'b');
% eglm_plot_fit(fit_1PostDCZ,'params',plotparams,'plotmode',plotmode);
plot([-1 1]*10000, [0 0], '-k','xliminclude','off')
text(numel(xname)-1, 2, 'Control', 'color', 'k');
text(numel(xname)-1, 1.9, 'Control injection (DCZ+DMSO)', 'color', [0.5 0.5 0.5]);
text(numel(xname)-1, 1.8, 'DCZ', 'color', 'r');

xtickangle(45);

% set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figure and orientation
% 
% destination = 'C:\Users\Takaya\Dropbox\RewPredNov\Workflow\Figure\';
% filename = 'GLMfit';
% print([destination filename], '-djpeg');




%% -- Ethan's test with Ilya's code

%% -- setup
load([Dir 'datasave.mat']);

% Part 1: fit neural activity with the model

fracid = datasave(:,1);
y = datasave(:,2);

% convert original fractal ID#s in file to sequential fractal indexes into
% our lists of their properties (e.g. ufractals_er, _sdr, ...)
fracindex = nans(size(fracid));
for ufi = 1:numel(ufractals)
    ok = fracid == ufractals(ufi);
    fracindex(ok) = ufi;
end
assert(~any(isnan(fracindex(:))),'could not convert from original to new fractal ID#s!');

% use the same model as testmode == 5 above, 
% i.e.: without laterality, and no ErxSDR, no NxSDN
assert(testmode == 5);

xreg = {};
xreg{end+1} = struct('name','ER','terms',{{ufractals_er(fracindex)'}});
xreg{end+1} = struct('name','SDR','terms',{{ufractals_sdr(fracindex)'}});
xreg{end+1} = struct('name','Nv','terms',{{ufractals_en(fracindex)'}});
xreg{end+1} = struct('name','SDNv','terms',{{ufractals_sdn(fracindex)'}});

xreg{end+1} = struct('name','ERxNv','terms',{{ufractals_er(fracindex)',ufractals_en(fracindex)'}});
xreg{end+1} = struct('name','ERxSDNv','terms',{{ufractals_er(fracindex)',ufractals_sdn(fracindex)'}});
xreg{end+1} = struct('name','SDRxNv','terms',{{ufractals_sdr(fracindex)',ufractals_en(fracindex)'}});
xreg{end+1} = struct('name','SDRxSDNv','terms',{{ufractals_sdr(fracindex)',ufractals_sdn(fracindex)'}});

% ...but add a constant factor, as is necessary for fitting neural activity
xreg{end+1} = struct('name','1','terms',{{'flag',ones(size(fracid,1),1)}});

yreg = struct('name','rate','value',y);

% do the fit
fit_neuron = eglm_fit(xreg,yreg,'normal');


% Part 2: get total offer value (and its components) for the offers, according to the earlier big behavioral fit

% First, we make the single neuron's fit data structure
% be in same basic format as earlier big behavioral fit,
% so we can use 'apply existing fit' to apply that earlier fit to this
% data. Specifically:
% - duplicate the offer, to make it 'as if' there are 2 offers per trial
% - remove the constant term
% - make up dummy choice data (i.e. random choices)
dummy_xreg = xreg;
for xi = 1:numel(dummy_xreg)
    for ti = 1:numel(dummy_xreg{xi}.terms)
        if ~ischar(dummy_xreg{xi}.terms{ti})
            dummy_xreg{xi}.terms{ti} = [dummy_xreg{xi}.terms{ti} dummy_xreg{xi}.terms{ti}];
        end
    end
end
assert(strcmp(dummy_xreg{end}.name,'1'));
dummy_xreg = dummy_xreg(1:(end-1));
dummy_yreg = struct('name','DUMMY choice','value',rand < 0.5*ones(numel(fracindex),1));

% which per-offer effects to calculate. This is my stab at ones that may be
% interesting to correlate with neural actvity
peroffer_effects = struct( ...
    'value','*', ...
    'reward_only',{{'ER','SDR'}}, ...
    'novelty_only',{{'Nv','SDNv'}}, ...
    'reward_x_novelty',{{'ERxNv','ERxSDNv','SDRxNv','SDRxSDNv'}});

% NOTE: this code assumes the neuron was recorded during Slayer and hence 
% uses Slayer's big behavioral fit. Edit this as needed.
big_behav_fit = fit_ALlchoiceSl;

% evaluate the big behavioral model on the dummy data from this neuron's recording
fit_dummy = eglm_fit(dummy_xreg,dummy_yreg,'binary choices','apply existing fit',big_behav_fit,'calculate per-offer effects',peroffer_effects);

% get the offer effects to use for future analysis, 
% after removing the second, dummy duplicate copy of the offer
peroff = fit_dummy.peroffer;
fnames = fieldnames(peroff);
for fni = 1:numel(fnames)
    peroff.(fnames{fni}) = peroff.(fnames{fni})(:,1);
end

%% -- plot the results

figuren;
eglm_plot_fit(fit_neuron);

figuren;
ncol = numel(fnames);

for fni = 1:numel(fnames)
    fn = fnames{fni};
    
    nsubplot(3,ncol,1,fni);
    plot(fit_neuron.yreg.value,'ko');
    xlabel('Trial #');
    ylabel('Neuron rate');
    
    nsubplot(3,ncol,2,fni);
    plot(peroff.(fn),'ko');
    xlabel('Trial #');
    ylabel(fn,'interpreter','none');
    title(sprintf('GLM summed value of: "%s"',fn),'interpreter','none');

    nsubplot(3,ncol,3,fni);
    x = peroff.(fn);
    y = fit_neuron.yreg.value;
    
    plot(x,y,'ko');
    
    [rho,pval] = corr(x,y,'type','Spearman');
    [slope,yint] = type_2_regression(x,y);
    linexy([0 1],yint + slope.*[0 1],'k-');
    
    xlabel(fn,'interpreter','none');
    ylabel('Neuron rate');
    title('Neural-value correlation');
    
    etextn('lt',sprintf('rho %.2f\np %.4f\nn %d',rho,pval,numel(x)));
end