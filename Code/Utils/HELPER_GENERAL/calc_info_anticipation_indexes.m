function [stats] = calc_info_anticipation_indexes(rates_precue_postcue_preout,cond_info,cond_unc,nperms)
% [indexes,pvals,rocareas] = calc_info_anticipation_indexes(rates_precue_postcue_preout,cond_info,cond_unc,nperms)
%
% inputs:
%  - rates_precue_postcue_preout: Ntrials x 3 matrix of firing rates
%    Each row is a trial. The columns are the firing rates in the epochs:
%     column 1: pre-cue (a 500ms epoch ending at cue onset)
%     column 2: post-cue (a 500ms epoch starting 150ms after cue onset)
%     column 3: pre-outcome (a 500ms epoch ending at outcome onset)
%
%  - cond_info: Ntrials x 1 logical vector indicating each trial's infoness
%     true for Info trials, false for Noinfo trials
%
%  - cond_unc: Ntrials x 1 logical vector indicating each trial's uncertainty
%     true for Uncertain trials, false for Certain trials
%
%  - nperms: number of permutations to use for permutation tests
%
% output: data structure 'stats' with two fields, 
%   infocue: Info Cue Anticipation Index
%   uncout: Uncertain Outcome Anticipation Index
%  Each has several subfields with its various statistics:
%   index: the Index
%   p: its p-value (permutation test)
%   roc: a 1 x 2 vector with the Index's component ROCs
%    (the index is defined as roc(1) - roc(2)).
%   rocp: a 1 x 2 vector with the p-values for the two component ROCs
%
% Created 2019-04-29 by ESBM for use by Ilya & colleagues.

if nargin < 4
    nperms = 20000;
end


%% setup analysis parameters

default_struct = struct( ...
    'name','', ...
    'params', {{}}, ...
    'nperms', nperms, ...
    'index', nan, ...
    'p', nan, ...
    'roc', [nan nan], ...
    'rocp', [nan nan]);

% params{i} specifies which epoch & trials to use for roc i, which is
%  computed based on the contrast of trials specified in
%  params{i}{1} (signal) vs params{i}{2} (noise)
%
%  params{i}{j}{1} specifies the epoch ('precue', 'postcue', or 'preout')
%  params{i}{j}{2} specifies the infoness condition (0 or 1)
%  params{i}{j}{3} specifies the uncertainty condition (0 or 1)
%
% the index is then calculated as roc(1) - roc(2)

stats = struct();

% ROCprecue(Info Unc,Info Cer) - ROCprecue(Noinfo Unc, Noinfo Cer)
stats.infocue = default_struct;
stats.infocue.name = 'Info Cue Anticipation Index';
stats.infocue.params = { ...
    { {'precue',1,1}, {'precue',1,0} }, ...
    { {'precue',0,1}, {'precue',0,0} }, ...
    };

% ROCpreout(Noinfo Unc,Noinfo Cer) - ROCpostcue(Noinfo Unc,Noinfo Cer)
stats.uncout = default_struct;
stats.uncout.name  = 'Uncertain Outcome Anticipation Index';
stats.uncout.params = { ...
    { {'preout',0,1}, {'preout',0,0} }, ...
    { {'postcue',0,1}, {'postcue',0,0} }, ...
    };

%% sanity check the inputs
ntr = size(rates_precue_postcue_preout,1);
assert(isequaln(size(rates_precue_postcue_preout),[ntr 3]),'rates_precue_postcue_preout must be (Ntrials x 3) matrix!');
assert(isequaln(size(cond_info),[ntr 1]) & isequaln(size(cond_unc),[ntr 1]),'cond_info and cond_unc must be (Ntrials x 1) vectors!');

if ntr < 1
    disp('calc_info_anticipation_indexes: warning: inputs had no data, will return NaN results!');
    return;
end

if ~islogical(cond_info)
    assert(isnumeric(cond_info),'cond_info must be logical or at least numeric!');
    cond_info = cond_info == 1;
end
if ~islogical(cond_unc)
    assert(isnumeric(cond_unc),'cond_unc must be logical or at least numeric!');
    cond_unc = cond_unc == 1;
end

% need at least this many trials in each condition
min_ntr_in_cond = 1;

ntr_in_cond = [ ...
    sum(~cond_info & ~cond_unc) ...
    sum(~cond_info &  cond_unc) ...
    sum( cond_info & ~cond_unc) ...
    sum( cond_info &  cond_unc) ...
    ];

has_enough_trials = true;
if ~all(ntr_in_cond([1 2 4]) >= min_ntr_in_cond)
    has_enough_trials = false;
else
    if ntr_in_cond(3) == 0
        disp('calc_info_anticipation_indexes: warning: did not find any Info Certain trials but found all other trial types, will calculate Info Cue Anticipation Index using slightly different equation used for InfoMegaTask control task');

        % ROCprecue(Info Unc,Noinfo Cer) - ROCprecue(Noinfo Unc, Noinfo Cer)
        stats.infocue.params = { ...
            { {'precue',1,1}, {'precue',0,0} }, ...
            { {'precue',0,1}, {'precue',0,0} }, ...
            };
    else
        if ~(ntr_in_cond(3) >= min_ntr_in_cond)
            has_enough_trials = false;
        end
    end
end

if ~has_enough_trials
    disp('calc_info_anticipation_indexes: warning: must have at least %d trials per condition (Info or Noinfo) x (Uncertain or Certain), will return NaN results!');
    return;
end

%% do the analysis

fnames = fieldnames(stats);
for fni = 1:numel(fnames)
    fn = fnames{fni};
    
    % calculate component ROCs
    rates = cell(2,2);
    for roci = 1:2
        for di = 1:2
            % get rates from the specified epoch and trial conditions
            par = stats.(fn).params{roci}{di};
            switch par{1}
                case 'precue'
                    epochi = 1;
                case 'postcue'
                    epochi = 2;
                case 'preout'
                    epochi = 3;
                otherwise
                    error('unknown epoch "%s", expected ''precue'', ''postcue'', or ''preout''',par{1});
            end
            oktr = cond_info == par{2} & cond_unc == par{3};
            rates{roci,di} = rates_precue_postcue_preout(oktr,epochi);
        end
        
        % calculate ROC area and its p-value
        [stats.(fn).roc(roci),stats.(fn).rocp(roci)] = rocarea3(rates{roci,2},rates{roci,1});
    end
    
    % calculate Index and its p-value
    rate1sig = rates{1,1};
    rate1noi = rates{1,2};
    
    rate2sig = rates{2,1};
    rate2noi = rates{2,2};
    
    % construct proper input to permutation_test_general
    % to allow a permutation-based test of this index
    % (which shuffles the correct conditions to represent the null
    % hypothesis that the uncertainty signal does not change
    % between the two conditions+epochs being compared)
    addset = @(rate,lab,permgroup,newrate,newlab,newpermgroup) deal([rate ; newrate],[lab ; newlab*ones(numel(newrate),1)],[permgroup ; newpermgroup*ones(numel(newrate),1)]);

    rate = [];
    lab = [];
    permgroup = [];
    [rate,lab,permgroup] = addset(rate,lab,permgroup,rate1noi,1,1);
    [rate,lab,permgroup] = addset(rate,lab,permgroup,rate1sig,2,2);
    [rate,lab,permgroup] = addset(rate,lab,permgroup,rate2noi,3,1);
    [rate,lab,permgroup] = addset(rate,lab,permgroup,rate2sig,4,2);
    [stats.(fn).p,stats.(fn).index] = permutation_test_general(rate,lab,permgroup, ...
        @(rate,lab,permgroup) rocarea3(rate(lab==1),rate(lab==2))-rocarea3(rate(lab==3),rate(lab==4)), ...
        stats.(fn).nperms);
end
