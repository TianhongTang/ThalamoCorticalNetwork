function [startinglowfreq,endinglowfreq,startinghighfreq,endinghighfreq,goodnessvalue,superficialchannel,deepchannel,highfreqmaxchannel,lowfreqmaxchannel,crossoverchannel] = FLIPAnalysis(nonnormpowmat,laminaraxis,freqaxis,setfreqbool,alphafreqrange, gammafreqrange)
arguments
    nonnormpowmat (:, :) double
    laminaraxis (1, :) double
    freqaxis (1, :) double
    setfreqbool (1, 1) double
    % set default values if not passed by user
    alphafreqrange (1, 2) double = [10 19];
    gammafreqrange (1, 2) double = [75 150];
end

% FLIPANALYSIS takes a non-normalized power matrix as an input and computes
% the spectrolaminar cross and returns pertinent values.
%
% Required Inputs:
%   nonnormpowmat: the non normalized power matrix that you want to analyze
%       using this script. The FLIP algorithm requires the first dimension
%       to represent frequency and the second dimension to represent
%       laminar position. An example power matrix for 32 probe contacts
%       with frequencies of interest 0-150 Hz (frequency bin size 1 Hz)
%       would be a 32x151 matrix.
%
%   laminaraxis: (1 x n) row vector where n represents the number of probe
%       contacts represented in the power matrix. Example vector for a 32
%       channel probe with 100 micron spacing: 0:0.1:3.1
%
%   freqaxis: (1 x n) row vector where n represents the number of frequency
%       bins represented in the power matrix. Example vector for a power
%       matrix using frequencies 0 Hz to 150 Hz with frequency bins of 1
%       Hz: 0:151
%
%   setfreqbool: boolean value that, if true, sets frequency window to the
%       optimal values for alpha/beta and gamma found in the Ubiquitous
%       Spectolaminar Motif paper (10-19 Hz for alpha/beta and 75-150 Hz
%       for gamma). If false, the algorithm will look at every possible
%       frequency range to optimize the alpha/beta and gamma goodness of
%       fit values (sacrificing efficiency).
%
%   Optional Parameters:
%       alphafreqrange: Gives the user the option to set the alpha/beta
%       frequency range to a desired value. Format is a 1x2 array where
%       element 1,1 is the lower bound and 1,2 is the upper bound. i.e.
%       [10 19] would utilize 10-19 Hz to determine the alpha/beta region.
%
%       gammafreqrange: Gives the user the option to set the gamma
%       frequency range to a desired value.Format is a 1x2 array where
%       element 1,1 is the lower bound and 1,2 is the upper bound. i.e.
%       [75 150] would utilize 75-150 Hz to determine the gamma region.
%
%   Example Usage:
%       FLIPANALYSIS(nonnormpowmat,laminaraxis,freqaxis,setfreqbool)
%       runs the FLIP algorithm on the optimal frequency ranges found in
%       the paper (on the condition that setfreqbool == TRUE).
%
%       FLIPANALYSIS(nonnormpowmat,laminaraxis,freqaxis,setfreqbool)
%       runs the FLIP algorithm on every possible range of frequencies to
%       determine which returns the best k value value (on the condition
%       that setfreqbool == FALSE).
%
%       FLIPANALYSIS(nonnormpowmat,laminaraxis,freqaxis,setfreqbool, alphafreqrange, gammafreqrange)
%       allows the user to determine the frequency ranges over which the
%       FLIP algorithm will analyze. FLIP expects these ranges to be input
%       in a 1x2 array format: [lower range, upper range] (on the condition
%       that setfreqbool == TRUE).

% validate user inputs and return necessary values for FLIP
% local function, throws an error if the inputs are incompatible with FLIP
[freqbinsize, minrange] = validateUserInputs(nonnormpowmat, freqaxis, laminaraxis, alphafreqrange, gammafreqrange);

% check to determine whether user wants to use a set frequency range
% if not, run the VFLIP analysis to optimize goodness value
if ~(setfreqbool)
    % Create a data structure to store the results of the VFLIP algorithm
    pooled_results = NaN(1,11,1);
    bki=1;
    nanboolean=~isnan(squeeze(nonnormpowmat));
    startrow= find(nanboolean ==1,1);
    if sum(nanboolean)==0
        startinglowfreq =  nan;
        endinglowfreq =  nan;
        startinghighfreq =  nan;
        endinghighfreq =  nan;
        goodnessvalue = nan;
        startrow=nan;
        superficialchannel=nan;
        deepchannel=nan;
        highfreqmaxchannel=nan;
        lowfreqmaxchannel=nan;
        crossoverchannel=nan;
    else

        powerspectra_image = nonnormpowmat;

        %setting different frequency bins to iterate over
        low_freq = ceil([1,10, 20, 30, 40, 50,60, 70] / freqbinsize);
        gamma_freq = ceil([30, 40, 50, 60, 70, 80,90,100,110,120,130,140,150 ] / freqbinsize);
        pooling_index = 0;

        pooled_results = NaN(1, 11, 1529);

        %iterating through all frequency combinations
        f = waitbar(0, "Finding Optimal Frequencies");
        for low_freq_startindex = 1:length(low_freq)-1
            for low_freq_endindex = low_freq_startindex+1:length(low_freq)

                for gamma_freq_startindex = 1:length(gamma_freq)-1
                    if gamma_freq(gamma_freq_startindex) < low_freq(low_freq_endindex)
                        continue
                    else

                        for gamma_freq_endindex = gamma_freq_startindex+1:length(gamma_freq)
                            alpha = low_freq(low_freq_startindex):low_freq(low_freq_endindex);
                            gamma = gamma_freq(gamma_freq_startindex):gamma_freq(gamma_freq_endindex);
                            % for each range of alpha/beta and gamma windows,
                            % calculate the optimal channel range that produces the highest goodness value
                            [goodness_k_value,upright_boolean,superficialchannel,deepchannel,gammamaxchannel,betamaxchannel,final_crossoverchannel] = evaluateFrequencyRange(powerspectra_image,alpha,gamma,minrange);
                            pooling_index = pooling_index+1;
                            %inputting the goodness value, frequency bin,
                            %information to the pooled_results matrix
                            pooled_results(bki,1,pooling_index) =low_freq(low_freq_startindex)*freqbinsize;
                            pooled_results(bki,2,pooling_index) =low_freq(low_freq_endindex)*freqbinsize;


                            pooled_results(bki,3,pooling_index) =gamma_freq(gamma_freq_startindex)*freqbinsize;
                            pooled_results(bki,4,pooling_index) =gamma_freq(gamma_freq_endindex)*freqbinsize;

                            pooled_results(bki,5,pooling_index) =goodness_k_value;
                            pooled_results(bki,6,pooling_index) =startrow;

                            %setting nans for non-identifiable probes
                            if isempty(superficialchannel)
                                superficialchannel=nan;
                            end
                            if isempty(deepchannel)
                                deepchannel=nan;
                            end

                            if isempty(final_crossoverchannel)
                                final_crossoverchannel = nan;
                            end


                            if isempty(gammamaxchannel)
                                gammamaxchannel = nan;
                            end

                            if isempty(betamaxchannel)
                                betamaxchannel = nan;
                            end
                            %inputting other information about each frequency
                            %pairing

                            pooled_results(bki,7,pooling_index) =superficialchannel;
                            pooled_results(bki,8,pooling_index) =deepchannel;
                            pooled_results(bki,9,pooling_index) =gammamaxchannel;
                            pooled_results(bki,10,pooling_index) =betamaxchannel;
                            pooled_results(bki,11,pooling_index) =final_crossoverchannel;

                        end
                    end
                end
            end
            waitbar(low_freq_startindex/(length(low_freq)-1), f, sprintf('Progress: %d %%', floor((low_freq_startindex/(length(low_freq)-1))*100)))
        end
        close(f);

        probe1results = squeeze(pooled_results(bki,:,:));


        if  isnan(max(abs(probe1results(5,:)))) || max(abs(probe1results(5,:)))==0
            startinglowfreq =  nan;
            endinglowfreq =  nan;
            startinghighfreq =  nan;
            endinghighfreq =  nan;
            goodnessvalue = nan;
            startrow=nan;
            superficialchannel=nan;
            deepchannel=nan;
            highfreqmaxchannel=nan;
            lowfreqmaxchannel=nan;
            crossoverchannel=nan;

        else
            %finidng the optimal frequency pairing that results in the highest
            %magnitude g value.


            indexnum = find(abs(probe1results(5,:))== max(abs(probe1results(5,:))));
            startinglowfreq =  probe1results(1,indexnum);
            endinglowfreq =  probe1results(2,indexnum);
            startinghighfreq =  probe1results(3,indexnum);
            endinghighfreq =  probe1results(4,indexnum);
            goodnessvalue = probe1results(5,indexnum);
            startrow=probe1results(6,indexnum);
            superficialchannel=probe1results(7,indexnum);
            deepchannel=probe1results(8,indexnum);
            highfreqmaxchannel=probe1results(9,indexnum);
            lowfreqmaxchannel=probe1results(10,indexnum);
            crossoverchannel=probe1results(11,indexnum);

            alphafreqrange = [startinglowfreq endinglowfreq];
            gammafreqrange = [startinghighfreq endinghighfreq];

            plot_result(nonnormpowmat, lowfreqmaxchannel, highfreqmaxchannel, crossoverchannel, goodnessvalue, setfreqbool, superficialchannel, deepchannel, alphafreqrange, gammafreqrange);
        end
    end

else
    alpha = alphafreqrange(1):alphafreqrange(2);
    gamma = gammafreqrange(1):gammafreqrange(2);
    [goodness_k_value,upright_boolean,superficialchannel,deepchannel,gammamaxchannel,betamaxchannel,final_crossoverchannel] = evaluateFrequencyRange(nonnormpowmat,alpha,gamma,minrange);
    startinglowfreq = alphafreqrange(1);
    endinglowfreq = alphafreqrange(2);
    startinghighfreq = gammafreqrange(1);
    endinghighfreq = gammafreqrange(2);
    goodnessvalue = goodness_k_value;
    highfreqmaxchannel = gammamaxchannel;
    lowfreqmaxchannel = betamaxchannel;
    crossoverchannel = final_crossoverchannel;

    % plot_result(nonnormpowmat, lowfreqmaxchannel, highfreqmaxchannel, crossoverchannel, goodnessvalue, setfreqbool, superficialchannel, deepchannel, alphafreqrange, gammafreqrange);
end
end

function [freqbinsize, minrange] = validateUserInputs(nonnormpowmat, freqaxis, laminaraxis, alphafreqrange, gammafreqrange)
% determine the size of the frequency bins in powerspectrum
freqSteps = diff(freqaxis);
if (all(freqSteps ~= freqSteps(1)))
    error("Error using FLIPAnalysis: frequency values must be consistently spaced");
elseif (length(freqaxis) ~= size(nonnormpowmat, 2))
    error("Error using FLIPAnalysis: lengths of power matrix and frequency axis are inconsistent")
elseif freqaxis(1) > freqaxis(length(freqaxis))
    error("Error using FLIPAnalysis: frequency axis must be in ascending order")
else
    freqbinsize = freqSteps(1);
end

% first, check to assure laminar axis is correct format for FLIP
% if laminar axis is in correct format, determine minimum channel range
laminarSteps = diff(laminaraxis);
if (all(laminarSteps ~= laminarSteps(1)))
    error("Error using FLIPAnalysis: laminar values must be consistently spaced");
elseif (length(laminaraxis) ~= size(nonnormpowmat, 1))
    error("Error using FLIPAnalysis: lengths of power matrix and laminar axis are inconsistent")
elseif laminaraxis(1) > laminaraxis(length(laminaraxis))
    error("Error using FLIPAnalysis: laminar axis must be in ascending order")
elseif laminaraxis(1) ~= 0
    error("First element of laminar axis must be 0")
else
    if max(laminaraxis) < 0.7
        error("Error using FLIPAnalysis: laminar axis must be larger than 0.6 mm")
    end
    ind = 1;
    while laminaraxis(ind) < 0.7
        ind = ind + 1;
    end
    minrange = ind;
end

% ensure both frequency ranges are 1x2 row vectors
if ~(isvector(alphafreqrange) && numel(alphafreqrange) == 2 && size(alphafreqrange, 1) == 1)
    error("Size of alpha/beta frequency range must be a 1x2 row vector")
elseif ~(isvector(gammafreqrange) && numel(gammafreqrange) == 2 && size(gammafreqrange, 1) == 1)
    error("Size of gamma frequency range must be a 1x2 row vector")
end
end

function [goodness_k_value,upright_boolean,superficialchannel,deepchannel,gammamaxchannel,betamaxchannel,final_crossoverchannel] = evaluateFrequencyRange(powerspectra_image,alphabeta,gamma,minrange)
set_pval=0.05;
individual_goodness_result = nan(3,1);
nanboolean=~isnan(powerspectra_image(:,1));
startrow= find(nanboolean ==1,1);
endrow = startrow+sum(nanboolean)-1;

for superficialchannel = startrow: endrow-minrange
    for deepchannel = superficialchannel+minrange:endrow
        powspec_window = powerspectra_image(superficialchannel:deepchannel,:);

        maxpow =squeeze(max(powspec_window,[],1));
        relpow = squeeze(powspec_window)./maxpow;
        lowband = mean(relpow(:,alphabeta),2);
        highband = mean(relpow(:,gamma),2);



        %extra normalization toggle
        %lowband = lowband-min(lowband);
        %lowband = lowband/max(lowband);
        %highband = highband-min(highband);
        %highband = highband/max(highband);

        [lowband_b,bint,r,rint,stats_low] = regress(lowband,horzcat([1:length(lowband)]',ones(length(lowband),1)));
        lowband_rsquared=stats_low(1);
        lowband_pval=stats_low(3);

        [highband_b,bint,r,rint,stats_high] = regress(highband,horzcat([1:length(highband)]',ones(length(highband),1)));
        highband_rsquared=stats_high(1);
        highband_pval=stats_high(3);
        reward = 0.04*(deepchannel-superficialchannel+1)+0.72;

        beta_peak_locations = find((lowband ==max(lowband)));
        gamma_peak_locations = find((highband ==max(highband)));

        if mean(beta_peak_locations)> length(lowband)/2
            beta_peak_index = max(find(lowband ==max(lowband)));
        else
            beta_peak_index = min(find(lowband ==max(lowband)));
        end

        if mean(gamma_peak_locations)> length(highband)/2
            gamma_peak_index= max(find(highband ==max(highband)));
        else
            gamma_peak_index= min(find(highband ==max(highband)));
        end

        if beta_peak_index==1
            beta_peak_max_check = [superficialchannel==startrow & lowband(beta_peak_index)>lowband(beta_peak_index+1)];
        elseif beta_peak_index==length(lowband)
            beta_peak_max_check =[deepchannel ==endrow & lowband(beta_peak_index)>lowband(beta_peak_index-1)];
        else
            beta_peak_max_check = [lowband(beta_peak_index)>lowband(beta_peak_index+1) & lowband(beta_peak_index)>lowband(beta_peak_index-1)];
        end

        if gamma_peak_index==1
            gamma_peak_max_check = [superficialchannel==startrow & highband(gamma_peak_index)>highband(gamma_peak_index+1)];
        elseif gamma_peak_index==length(highband)
            gamma_peak_max_check =[deepchannel ==endrow & highband(gamma_peak_index)>highband(gamma_peak_index-1)];
        else
            gamma_peak_max_check = [highband(gamma_peak_index)>highband(gamma_peak_index+1) & highband(gamma_peak_index)>highband(gamma_peak_index-1)];
        end

        if gamma_peak_max_check & beta_peak_max_check
        else
            goodness=nan;
            individual_goodness_result = horzcat(individual_goodness_result, [goodness; superficialchannel;deepchannel]);
            continue
        end


        if lowband_pval<set_pval & highband_pval<set_pval
            if lowband_b(1)>0 & highband_b(1)<0
                goodness = highband_rsquared*lowband_rsquared*reward;
            elseif lowband_b(1)<0 & highband_b(1)>0
                goodness = -highband_rsquared*lowband_rsquared*reward;
            else
                goodness=nan;
                continue
            end
            individual_goodness_result = horzcat(individual_goodness_result, [goodness; superficialchannel;deepchannel]);
        else
            continue
        end
    end
end


a = max((individual_goodness_result(1,:)));
b = min((individual_goodness_result(1,:)));
if a>=abs(b)

    goodness_k_value = a ;
    superficialchannel = individual_goodness_result(2,find(individual_goodness_result(1,:) ==a,1));
    deepchannel=individual_goodness_result(3,find(individual_goodness_result(1,:) ==a,1));
    upright_boolean=1;



else
    goodness_k_value = b ;
    superficialchannel = individual_goodness_result(2,find(individual_goodness_result(1,:) ==b,1));
    deepchannel=individual_goodness_result(3,find(individual_goodness_result(1,:) ==b,1));
    upright_boolean=0;
end
powspec_window = squeeze(powerspectra_image(superficialchannel:deepchannel,:));

maxpow =squeeze(max(powspec_window,[],1));
relpow = squeeze(powspec_window)./maxpow;
lowband = mean(relpow(:,alphabeta),2);
highband = mean(relpow(:,gamma),2);


beta_peak_locations = find(lowband ==max(lowband));
gamma_peak_locations = find(highband ==max(highband));

if mean(beta_peak_locations)> length(lowband)/2
    beta_peak_index = max(find(lowband ==max(lowband)));

else
    beta_peak_index = min(find(lowband ==max(lowband)));
end
if mean(gamma_peak_locations)> length(highband)/2
    gamma_peak_index= max(find(highband ==max(highband)));
else

    gamma_peak_index= min(find(highband ==max(highband)));
end

betamaxchannel = superficialchannel+beta_peak_index-1;
gammamaxchannel = superficialchannel+gamma_peak_index-1;


band_diff=highband-lowband;
band_diff_inverted = lowband-highband;
crossoverchannels=[];
for channel_index = 1:length(lowband)-2
    if goodness_k_value>0

        %gamma is greater than beta, and then it crosses over to
        %beta being greater than gamma
        if highband(channel_index)>lowband(channel_index) && lowband(channel_index+1)>highband(channel_index+1)
            if  abs(band_diff(channel_index))<= abs(band_diff(channel_index+1))
                crossoverchannels= [crossoverchannels channel_index];
            else
                crossoverchannels=[crossoverchannels channel_index+1];
            end


        elseif highband(channel_index)>lowband(channel_index) &&highband(channel_index+1)==lowband(channel_index+1) && lowband(channel_index+2)>highband(channel_index+2)
            crossoverchannels=[crossoverchannels channel_index+1];
        end


    elseif goodness_k_value<0
        if highband(channel_index)<lowband(channel_index) && lowband(channel_index+1)<highband(channel_index+1)
            if  abs(band_diff(channel_index))<= abs(band_diff(channel_index+1))
                crossoverchannels=[crossoverchannels channel_index];
            else
                crossoverchannels=[crossoverchannels channel_index+1];
            end


        elseif highband(channel_index)<lowband(channel_index) &&highband(channel_index+1)==lowband(channel_index+1) && lowband(channel_index+2)<highband(channel_index+2)
            crossoverchannels=[crossoverchannels channel_index+1];
        end
    end
end
%now, all valid cross over points have been collected

number_of_crosses = length(crossoverchannels);
crossover_ratings = horzcat(crossoverchannels',nan(number_of_crosses,1));
if number_of_crosses ==0
    final_crossoverchannel =nan;
elseif number_of_crosses==1
    final_crossoverchannel =crossoverchannels+superficialchannel-1;
else
    if goodness_k_value>0

        for index = 1: number_of_crosses
            crossover_choice = crossoverchannels(index);
            % calculates difference between highband and lowband on
            % each side of the crossover
            crossover_ratings(index,2) = sum(band_diff(1:crossover_choice))-sum(band_diff(crossover_choice:length(lowband)));
        end


    elseif goodness_k_value<0
        for index = 1: number_of_crosses
            crossover_choice = crossoverchannels(index);
            crossover_ratings(index,2) = sum(band_diff_inverted(1:crossover_choice))-sum(band_diff_inverted(crossover_choice:length(lowband)));
        end
    end
    %return index of ideal crossover point
    [~,indexval] = max(crossover_ratings(:,2));

    final_crossoverchannel = crossover_ratings(indexval,1)+superficialchannel-1;
end
end

function plot_result(nonnormpowmat, lowfreqmaxchannel, highfreqmaxchannel, crossoverchannel, goodnessvalue, setfreqbool, superficialchannel, deepchannel, lowfreqrange, highfreqrange)
arguments
    nonnormpowmat (:, :) double
    lowfreqmaxchannel (1, 1) double
    highfreqmaxchannel (1, 1) double
    crossoverchannel (1, 1) double
    goodnessvalue (1, 1) double
    setfreqbool (1, 1) double
    superficialchannel (1, 1) double
    deepchannel (1, 1) double
    % set default values if not passed by user
    lowfreqrange (1, 2) double = [10 19]; % default values
    highfreqrange (1, 2) double = [75 150]; % default values
end

maxpow =squeeze(max(nonnormpowmat,[],1));
relpow = squeeze(nonnormpowmat)./maxpow;

nchannels = size(nonnormpowmat, 1);

%%% Plot relpow map
figure;
subplot(1,3,[1 2]); hold on;
imagesc(relpow);set(gca, 'YDir', 'reverse');
xlim([1 size(relpow,2)]); ylim([1 size(relpow,1)]); title('LFP relative power; \beta-Freq Range:blue \gamma-Freq Range:red');
xlabel('Frequency (Hz)');ylabel('Channel Number');
cb=colorbar; caxis([0.3 1]); ylabel(cb, 'Relative Power');
% plot low freq range
line([lowfreqrange(1) lowfreqrange(1)],get(gca,'YLim'),'LineWidth',0.5,'Color','b','LineStyle','--');
line([lowfreqrange(2) lowfreqrange(2)],get(gca,'YLim'),'LineWidth',0.5,'Color','b','LineStyle','--');
% plot low freq range
line([highfreqrange(1) highfreqrange(1)],get(gca,'YLim'),'LineWidth',0.5,'Color','r','LineStyle','--');
line([highfreqrange(2) highfreqrange(2)],get(gca,'YLim'),'LineWidth',0.5,'Color','r','LineStyle','--');

line(get(gca, 'xlim'), [crossoverchannel crossoverchannel], 'LineWidth', 0.5, 'Color', 'k', 'LineStyle', '-.');
text(size(relpow,2), crossoverchannel, 'Crossover Channel', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
line(get(gca, 'xlim'), [lowfreqmaxchannel lowfreqmaxchannel], 'LineWidth', 0.5, 'Color', 'k', 'LineStyle', '-.');
text(size(relpow,2), lowfreqmaxchannel, 'Alpha/Beta Peak Channel', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
line(get(gca, 'xlim'), [highfreqmaxchannel highfreqmaxchannel], 'LineWidth', 0.5, 'Color', 'k', 'LineStyle', '-.');
text(size(relpow,2), highfreqmaxchannel, 'Gamma Peak Channel', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

%%% Plot banded relpow
subplot(1,3,3); hold on;
plot(mean(relpow(:,lowfreqrange(1):lowfreqrange(2)),2), 1:nchannels, 'b', 'LineWidth', 2); % alpha/beta
plot(mean(relpow(:,highfreqrange(1):highfreqrange(2)),2), 1:nchannels, 'r', 'LineWidth', 2); % gamma
set(gca, 'ydir', 'reverse')
xlim([0 1]); ylim([1 nchannels]); xlabel('Relative Power'); ylabel('Channel Number');

line(get(gca, 'xlim'), [crossoverchannel crossoverchannel], 'LineWidth', 0.5, 'Color', 'k', 'LineStyle', '-.');
text(1, crossoverchannel, 'Crossover Channel', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
line(get(gca, 'xlim'), [lowfreqmaxchannel lowfreqmaxchannel], 'LineWidth', 0.5, 'Color', 'k', 'LineStyle', '-.');
text(1, lowfreqmaxchannel, 'Alpha/Beta Peak Channel', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
line(get(gca, 'xlim'), [highfreqmaxchannel highfreqmaxchannel], 'LineWidth', 0.5, 'Color', 'k', 'LineStyle', '-.');
text(1, highfreqmaxchannel, 'Gamma Peak Channel', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
patch([0.95 1 1 0.95], [superficialchannel superficialchannel deepchannel deepchannel], 'yellow', 'FaceAlpha', 0.1)

title(['mean relpow; \beta:blue \gamma:red Regression Range:yellow']);

set(gcf,'Position',[200 200 1000 500]);
if (lowfreqrange == [10 19] & highfreqrange == [75 150])
    sgtitle(['Default Frequency Bin Values. G = ' num2str(goodnessvalue)]);
elseif setfreqbool == 1
    sgtitle(['Set Frequency Bin Values. G = ' num2str(goodnessvalue)]);
else
    sgtitle(['vFLIP Crossover. G = ' num2str(goodnessvalue)]);
end
end
% description of function outputs
%   startinglowfreq is the starting value of the optimal low frequency band
%   endinglowfreq is the ending value of the optimal low frequency band
%   startinghighfreq is the starting value of the optimal high frequency band
%   endinghighfreq is the ending value of the optimal high frequency band
%   goodnessvalue is the optimal G goodness value detected by v flip
%   superficial channel is the superficial channel of the optimal band
%   deep channel is the deep channel of the optimal band?
%   highfreqmaxchannel is the channel where the high freq band is at a max
%   lowfreqmaxchannel is the channel where the low freq band is at a max
%   crossoverchannel is the channel where the low and high bands cross over
