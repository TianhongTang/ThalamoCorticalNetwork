function fit = eglm_fit(xreg, yreg, fitting_style, varargin)
% fit = eglm_fit(xreg, yreg, fitting_style[,...])
%
% Fit a GLM to data, using code originally meant for binary choices between
% offers that have multiple features that can contribute to the decision.
% The main inputs are data structures specifying the x (regressors) and y
% (variabled to be predicted), and the type of GLM to use.
%
% Note that this does NOT automatically add a regressor that is constant.
% If you want a constant regressor, you need to add it to the list of
% regressors in xreg. This is because in our usual type of models we
% already have a constant regressor built in (e.g. main effect of offer
% presentation order)
%
% inputs:
%  xreg - a Nregressors-length cell array of specs for each regressor. Each
%         entry xi should be a data structure with the following field:
%
%         xreg{xi}.terms: a cell array where each entry is either one of
%          the individual underlying terms whose product will be the 
%          regressor, or a 'flag' variable indicating that the following
%          terms should be treated as flag variables. The terms should all
%          be either (Ntrials x 1) vectors, or (Ntrials x 2) matrices where
%          the first and second columns represent the properties of the
%          first and second offers on that trial.
%
%         Optionally, it can have the following fields:
%
%         xreg{xi}.normalization: specifies what normalization style to
%          use. Options are:
%
%           'default' - apply z-scoring normaliation to all terms except the
%            flag variables.
%
%           'only interaction terms' - same as default, but do not 
%            normalize the first term in the list of terms. So all effects
%            will be essentially in units of 'effect of change in first 
%            listed variable, PER-STANDARDIZED-CHANGE in the other listed 
%            variables'.
%
%           'none' - do not normalize any of the terms.
%
%           {'oktrials', Z} - a vector Z (logical or numeric) of indexes of 
%            trials. Same effect as default, but when computing the mean 
%            and SD for the z-scoring normalization of each term, uses all 
%            trials specified by Z, NOT the 'oktrials' vector that is 
%            normally used to restrict the analysis.
%
%            This option is meant for cases where e.g. we are doing
%            separate regressions for each of 20 separate behavioral
%            sessions, but in each of those fits we want to normalize the
%            terms of each regressor based on their overall distribution in
%            ALL 20 of the sessions (to ensure that they are perfectly
%            comparable across sessions).
%
%         xreg{xi}.manually_specified_beta: forces the regression weight to
%          be a specific value. That beta weight's SE and p-value will be
%          set to NaN. Use this when you have a set of models in the same
%          general framework, but some of them operate by forcing specific
%          parameters to be fixed weights rather than freely fit weights.
% 
%  yreg - the variable to be predicted. Either a (Ntrials x 1) column
%         vector of its values on each trial, or a data structure with 
%         that stored in its yreg.value field.
%
%  fitting_style - specifies the type of GLM. Can be either a 2-element
%         cell array of strings in the format
%         {distribution_name,link_function_name}, or else a string
%         specifying one of the following common fitting styles for our
%         data:
%          'binary choices' = {'binomial','logit'}
%          'poisson' = {'poisson','identity'}
%          'normal' = {'normal','identity'}
%
% Optional inputs, specified as a string indicating the name of the input
%  parameter, followed by its value:
%
%  'name',name - give the fit a name. This will be stored in fit.name. May
%         be useful for reference e.g. when you have a large
%         collection of fits.
%
%  'oktrials',oktrials - a (Ntrials x 1) logical vector or a (k x 1) 
%         numeric vector (where k <= Ntrials) specifying a subset of 
%         trials. The analysis will be restricted to those trials.
%
%  'calculate predictions',mode - specify whether to calculate predicted
%         y values for each trial, log likelihoods of the data, and related
%         values. mode must be logical (true or false).
%
%  'crossvalidation',kfold - if kfold > 0, then in addition to the full
%         fit to the data, also do k-fold cross-validation and report
%         the results in the fit.crossval field of the output structure.
%         Note that this is ONLY done if the fit is in 'calculate
%         predictions' mode, since otherwise there are no predictions to
%         compare!
%
%  'shuffle',nshuffles - if nshuffles > 0, then in addition to the full
%         fit to the data, also do nshuffles fits on datasets in which the
%         association between x and y are randomly shuffled, and report
%         the results in the fit.shuffle field of the output structure.
%         Note that this is ONLY done if the fit is in 'calculate
%         predictions' mode, since otherwise there are no predictions to
%         compare!
%
%  'shuffle ids',shufids - same, but user directly specifies 'shufids' as
%         the (# trials) x (# shuffles) shuffling matrix to be used to do
%         the shuffling. Each column specifies a single shuffle, and must 
%         have the entries (1:(#trials)), (presumably in a random order).
%
%  'rescale betas' - rescales the beta weights and their SEs so that
%         instead of being in the raw units of the predicted variable
%         ('s transformation by the link function), they are in units
%         of 'effect of changing this regressor by +1, AS A FRACTION OF 
%         the effect that would occur by changing the regressor named 
%         xreg_name by an amount equal to xreg_change'. For instance, when
%         predicting behavior in a 'pay per view' experiment, you may want
%         to estimate the subject's willingness to pay to obtain each 
%         regressor, so you may want to scale the beta weights to be in 
%         units of 'effect of +1 mL offered juice' or 'effect of +$1
%         offered money'. Specify this as a struct with the fields:
%
%          'name' - name of the new units (e.g. 'Effect on log odds of
%           choice relative to +1 mL juice')
%          'xreg_name' - name of the regressor to use for rescaling (e.g.
%           'E[Reward]')
%          'xreg_change' - number indicating the change in the regressor to
%           use for setting the amount of rescaling (e.g. if the regressor
%           is in un-normalized units of uL juice, then you could use 
%           xreg_change = 1000 to scale beta weights to be relative to the
%           effect of 1000 x +1 uL = +1 mL of juice)
%         
%
%  'calculate per-offer effects',effectinfo - in a choice situation with two 'offers'
%          on each trial, calculate the net effects of subsets of variables 
%          for each offer, specified in effectinfo. The results are stored
%          in fit.peroffer.
%
%          effectinfo is a struct with fields which are all cell arrays (or '*'),
%          calculates separate subsets, each of which has its name given by
%          the corresponding struct field, and has its regressors specified
%          by the regressor names listed inside. E.g.:
%           struct( ...
%            'reward',{{'E[r],'SD[r]'}}, ...
%            'spatial_bias',{{'left','right'}}, ...
%            'value','*')
%          would compute three peroffer effects, one which is the summed
%          effect of those two reward terms, and one which is the summed
%          effect of those spatial bias terms, and one which uses ALL
%          regressors that are specified for both offers (specified by
%          '*'). Note that the latter will NOT include effects of
%          regressors that are not specified per-offer (e.g. certain ways
%          of including a constant factor in the GLM)
%
%          This is done by default, with effectinfo = struct('value','*').
%          To turn it off, specify effectinfo as empty (i.e. as struct()).
%
%
%  'apply existing fit',efit - instead of fitting the data, apply the result
%          of an existing fit data structure (efit) to this data. That is,
%          instead of actually fitting a new model, just force the beta
%          weights to be the same as the existing model. This is useful
%          if you want to see the results of 'calculate predictions',
%          'calculate per-offer effects', etc. on a new dataset (or on a
%          specific subset of the existing dataset).
% 
%
% output:
%  fit - a data structure with the fitting results. Fields include:
%
%  fit.warning - information about whether a warning was detected while
%   calling glmfit to fit the GLM. Will print out a message if a warning is
%   detected, in which case you should be wary of trusting the results!
%
%  fit.name - name of the fit, if specified by the user
%
%  fit.param - basic parameters of the fit, including:
%
%  fit.param.distribution - GLM distribution function
%  fit.param.link - GLM link function
%
%  fit.param.oktrials - which trials in the data were used to fit the model
%  fit.param.normalization{xi} - parameters use to normalize the xi-th
%   regressor, including whether it is normalized, and the mean & sd that
%   were used to do the z-scoring normalization
%
%  fit.param.betascale - parameters used to do scaling of the beta weights
%   (if scaling was done). 
%  fit.param.betascale.name - always holds the name of the units of the
%   beta weights (e.g. "Effect on logs odds of choosing offer 2")
%
%  fit.param.crossval.kfold - # folds to use for cross-validation
%  fit.param.shuffle.n - # of shuffles to use for shuffle control
%
%  fit.xreg, fit.yreg - data specifying x (regressors) 
%   and y (variable to be predicted)
%
%  fit.x, fit.y - regressors (x) and variable to be predicted (y)
%  fit.b, fit.stats - outputs of glmfit
%
%  fit.glmfun_inputs - data structure holding the raw inputs given to
%   glmfit / glmval, which are the same as fit.x and fit.b in typical
%   cases, but may be different if you set certain regressors to have
%   betas that are manually-specified instead of fitted by the GLM.
%
%  The GLMs predictions (actually, in the current implementation, 
%  postdictions) include:
%
%  fit.ypred - predicted value of y on each trial
%  fit.r - correlation between predicted and actual values of y
%  fit.r2 - r squared between predicted vs actual values
%  fit.loglik_tr - log likelihood of y on each trial
%  fit.loglik - total log likelihood of the entire dataset
%  fit.mean_loglik_per_trial - log likelihood divided by # of data points
%  fit.p_correct_prediction - for binomial data, the fraction of trials
%   where the model predicts y 'correctly' (i.e. with probability > 0.5).
%  fit.n_free_parameters - number of free parameters in the fit (excluding
%   manually specified parameters)
%  fit.aic - Akaike's Information Criterion (2*k - 2*loglik)
%  fit.bic - Bayesian Information Criterion (k*log(n) - 2*loglik)
%
%  fit.crossval - structure with all of the above fields, but calculated using
%   k-fold cross validation
%
%  fit.shuffle - structure with all of the above fields, but calculated using
%   shuffled datasets. Also saves the shuffling matrix that was used in
%   fit.shuffle.shufids
%
%  fit.calibration - data structure holding results of 'calibration', i.e.
%   splitting the data into bins with different predicted values of y, and
%   testing whether the mean value of y in each group is close to its
%   predicted value. Key fields are:
%  fit.calibration.nbin - number of bins
%  fit.calibration.ypredmean - each bin's mean predicted value 
%  fit.calibration.ymean - each bin's mean observed value
%  fit.calibration.ysd - SD of each bin's observed value
%  fit.calibration.yse - SE of each bin's observed value
%
%  fit.peroffer - struct holding per-offer effects of subsets of regressors


fit = struct();
fit.name = '';

% whether to calculate predictions for each trial, 
% and the log likelihood of the data given the fitted model
fit.param.calculate_predictions = true;

% parameters for calculating per-offer effects of subsets of regressors
% (default: compute total value of each offer, using ALL regressors)
fit.param.peroffer = struct('value','*');

% whether to use a subset of the original data
fit.param.oktrials = 'default';

% whether to apply an existing fit instead of fitting the current data
fit.param.apply_existing_fit = false;
existing_fit = [];

% whether to rescale the fitted beta weights so that they represent the
% effect of a specific amount of change in a specific xregressor's
% raw data values (e.g. scaling from "log odds of biasing choice" to
% "equivalent effect of adding +1 mL of juice reward")
fit.param.betascale.name = [];
fit.param.betascale.xreg_name = [];
fit.param.betascale.xreg_change = [];
fit.param.betascale.rescaling_factor = [];

% default number of folds for crossvalidation
fit.param.crossval.kfold = 0;
fit.param.shuffle.use_shuffle_ids = [];
fit.param.shuffle.n = 0;

% read in variable arguments.
vi = 1;
while vi <= numel(varargin)
    switch varargin{vi}
        case 'name'
            assert(numel(varargin) >= vi+1 && ischar(varargin{vi+1}),'variable argument "%s" must be followed by a string',varargin{vi});
            fit.name = varargin{vi+1};
            vi = vi + 1;
        case 'calculate predictions'
            assert(numel(varargin) >= vi+1 && islogical(varargin{vi+1}),'variable argument "%s" must be followed by a logical',varargin{vi});
            fit.param.calculate_predictions = varargin{vi+1};
            vi = vi + 1;
        case 'calculate per-offer effects'
            assert(numel(varargin) >= vi+1 && isstruct(varargin{vi+1}),'variable argument "%s" must be followed by a struct',varargin{vi});
            fit.param.peroffer = varargin{vi+1};
            vi = vi + 1;
        case 'oktrials'
            assert(numel(varargin) >= vi+1 && (isnumeric(varargin{vi+1}) || islogical(varargin{vi+1})),'variable argument "%s" must be followed by a numeric or logical',varargin{vi});
            fit.param.oktrials = varargin{vi+1};
            vi = vi + 1;
        case 'crossvalidation'
            assert(numel(varargin) >= vi+1 && (isnumeric(varargin{vi+1}) || islogical(varargin{vi+1})),'variable argument "%s" must be followed by a numeric or logical',varargin{vi});
            fit.param.crossval.kfold = varargin{vi+1};
            vi = vi + 1;
        case 'shuffle'
            % shuffle, using a newly-generated random shuffling matrix
            assert(numel(varargin) >= vi+1 && (isnumeric(varargin{vi+1}) || islogical(varargin{vi+1})),'variable argument "%s" must be followed by a numeric or logical',varargin{vi});
            fit.param.shuffle.use_shuffle_ids = [];
            fit.param.shuffle.n = varargin{vi+1};
            vi = vi + 1;
        case 'shuffle ids'
            % shuffle, with a user-specified shuffling matrix
            assert(numel(varargin) >= vi+1 && isnumeric(varargin{vi+1}),'variable argument "%s" must be followed by a numeric or logical',varargin{vi});
            fit.param.shuffle.use_shuffle_ids = varargin{vi+1};
            fit.param.shuffle.n = size(fit.param.shuffle.use_shuffle_ids,2);
            vi = vi + 1;
        case 'rescale betas'
            assert(numel(varargin) >= vi+1 && isstruct(varargin{vi+1}) && all(isfield(varargin{vi+1},{'name','xreg_name','xreg_change'})),'variable argument "%s" must be followed by a data structure with fields "name", "xreg_name", "xreg_change"',varargin{vi});
            fit.param.betascale = varargin{vi+1};
            vi = vi + 1;
        case 'apply existing fit'
            assert(numel(varargin) >= vi+1 && isstruct(varargin{vi+1}),'variable argument "%s" must be followed by a struct',varargin{vi});
            fit.param.apply_existing_fit = true;
            existing_fit = varargin{vi+1};
            vi = vi + 1;
        otherwise
            error('unknown variable argument "%s"',varargin{vi});
    end
    vi = vi + 1;
end


% fitting style specifies the probability distribution the data is assumed
% to follow, and the link function. It is either specified as:
% {distribution,link function}, or shorthand for a commonly used such pair.
if iscell(fitting_style)
    assert(numel(fitting_style) == 2 && ischar(fitting_style{1}) && ischar(fitting_style{2}),'if fitting_style is a cell array, must have two strings indicating the distribution and the link function');
else
    switch fitting_style
        case 'binary choices'
            fitting_style = {'binomial','logit'};
        case 'poisson'
            fitting_style = {'poisson','identity'};
        case 'normal'
            fitting_style = {'normal','identity'};
        otherwise
            error('unknown fitting style');
    end
end
fit.param.distribution = fitting_style{1};
fit.param.link = fitting_style{2};


% if applying an existing fit, it must have the same params as the current fit
if fit.param.apply_existing_fit
    assert(~isempty(existing_fit) && numel(existing_fit)==1,'did not receive a single unique "existing_fit" data structure');
    assert(numel(xreg)==numel(existing_fit.xreg),'existing and current fits must have same number of xregs');
    for xi = 1:numel(xreg)
        assert(isequaln(xreg{xi}.name,existing_fit.xreg{xi}.name),'existing and current fits must have same xreg names');
        current_msb = isfield(xreg{xi},'manually_specified_beta');
        existing_msb = isfield(existing_fit.xreg{xi},'manually_specified_beta');
        assert(current_msb == existing_msb,'existing and current fits must have same xregs with/without manually specified beta');
        if current_msb
            assert(isequaln(xreg{xi}.manually_specified_beta,existing_fit.xreg{xi}.manually_specified_beta),'existing and current fit xregs with manually specified beta must have the same manually specified setting');
        end
    end
    assert(isequaln(fit.param.distribution,existing_fit.param.distribution) && isequaln(fit.param.link,existing_fit.param.link),'existing and current fits must have same distribution and link function');
    
    % replace the current fit's normalization settings with the 
    % existing fit's normalization settings
    for xi = 1:numel(xreg)
        xreg{xi}.normalization = {'mean and sd', existing_fit.param.normalization{xi}.mean, existing_fit.param.normalization{xi}.sd};
    end
end


% convert raw regressor inputs into actual regressors
fit.xreg = xreg;
fit.yreg = yreg;
[fit.x,fit.y,fit.param.normalization,fit.glmfun_inputs] = eglm_make_regression_inputs(fit.xreg,fit.yreg,fit.param.oktrials);
fit.b_is_manually_specified = ~isnan(fit.glmfun_inputs.manually_specified_beta);
fit.nx = size(fit.x,2);
fit.n = numel(fit.y);

% get names of x and y variables
fit.xname = cellfun(@(z) z.name,fit.xreg,'uniform',0);
if isstruct(fit.yreg)
    fit.yname = fit.yreg.name;
else
    fit.yname = 'predicted variable';
end


% verify that the specified peroffer calculations are valid
fnames = fieldnames(fit.param.peroffer);
for fni = 1:numel(fnames)
    fval = fit.param.peroffer.(fnames{fni});
    if ischar(fval)
        assert(strcmp(fval,'*'),'if per-offer effects are specified as a char, must be the wildcard (''*'')');
    elseif iscell(fval)
        for regi = 1:numel(fval)
            assert(ismember(fval{regi},fit.xname),'if per-offer effects are specified as a cell, must contain only valid names of xregressors');
        end
    else
        error('per-offer effects must be specified as ''*'' or a cell array of regressor names');
    end
end


% track whether the fit gives a warning, and if so, what type
fit.warning.detected = false;
fit.warning.msg = [];
fit.warning.id = [];
[lastmsg_before_fit,lastid_before_fit] = lastwarn;
lastwarn('');

if fit.param.apply_existing_fit
    
    % set fitting results based on an existing fit
    fit.glmfun_inputs.b = existing_fit.glmfun_inputs.b;
    fit.stats = existing_fit.stats;
    
    % adjust fit so that our variables reflect the ones that were manually
    % specified to be fixed values
    fit = eglm_adjust_fit_for_manually_specified_beta(fit);
    
else

    % fit the model
    [fit.glmfun_inputs.b,~,fit.stats] = eglm_glmfit(fit.glmfun_inputs.x,fit.y,fit.param.distribution,'link',fit.param.link,'constant','off','offset',fit.glmfun_inputs.offset);

    % adjust fit so that our variables reflect the ones that were manually
    % specified to be fixed values
    fit = eglm_adjust_fit_for_manually_specified_beta(fit);


    % check if warnings were given
    [fit.warning.msg,fit.warning.id] = lastwarn;
    fit.warning.detected = ~isempty(fit.warning.msg) || ~isempty(fit.warning.id);
    if fit.warning.detected
        % print out a warning, if this type of warning is not being suppressed
        warnstate = warning('query',fit.warning.id);
        if strcmp(warnstate.state,'on')
            fprintf(' eglm_fit detected at least one warning when calling glmfit to fit the GLM! Be wary of trusting the fitting results!\n');
            fprintf('  last warning msg: %s\n  last warning  id: %s\n',fit.warning.msg,fit.warning.id);
        end
    end
    % if no warning was detected, restore prior state of Matlab's warning flags
    if ~fit.warning.detected
        lastwarn(lastmsg_before_fit,lastid_before_fit);
    end

end

% if required, calculate model's prediction about each trial
% and calculate the log likelihood
if fit.param.calculate_predictions
    % calculate model's prediction about each data point
    fit.ypred = glmval(fit.glmfun_inputs.b,fit.glmfun_inputs.x,fit.param.link,'constant','off','offset',fit.glmfun_inputs.offset);

    fit = calc_model_predictions(fit,fit);
%     % calculate log likelihood of each individual trial
%     % based on the model's fitted parameters and its assumption about
%     % the distribution the data follow
%     fit.loglik_tr = nans(size(fit.ypred));
%     switch fit.param.distribution
%         case 'binomial'
%             fit.loglik_tr = log((fit.y==1).*fit.ypred + (fit.y==0).*(1-fit.ypred));
%         case 'poisson'
%             fit.loglik_tr = log(poisspdf(fit.y,fit.ypred));
%         case 'normal'
%             fit.loglik_tr = log(normpdf(fit.y,fit.ypred,fit.stats.sfit));
%         otherwise
%             error('have not implemented calculation of log likelihood for GLM distribution "%s"',fit.param.distribution);
%     end
% 
%     % total log likelihood of entire dataset
%     fit.loglik = sum(fit.loglik_tr);
% 
%     % mean log likelihood per trial
%     fit.mean_loglik_per_trial = fit.loglik ./ fit.n;
%     
%     % correlation between predictions and data
%     fit.r = corr(fit.ypred,fit.y,'type','Pearson');
%     % r-squared
%     sse_total = sum((fit.y - mean(fit.y)).^2);
%     sse_resid = sum(fit.stats.resid.^2);
%     fit.r2 = sse_resid ./ sse_total;
%     
%     % if binomial data, calculate p(correct prediction)
%     if strcmp(fit.param.distribution,'binomial')
%         fit.p_correct_prediction = mean(sign(fit.ypred-0.5) == sign(fit.y-0.5));
%     else
%         fit.p_correct_prediction = nans(fit.n,1);
%     end
    
    if fit.param.crossval.kfold > 0
        assert(~fit.param.apply_existing_fit,'does not make sense to do cross-validation if applying an existing fit, since the fit parameters are not being fit here in the first place, they are just being applied from an existing fit');
        fit.crossval = struct();
        fit.crossval.varnames = {'ypred','loglik_tr','loglik','mean_loglik_per_trial','r','r2','p_correct_prediction'};
        for vni = 1:numel(fit.crossval.varnames)
            vn = fit.crossval.varnames{vni};
            fit.crossval.(vn) = nans(size(fit.(vn)));
        end
        fit.crossval.sfit = nans(size(fit.crossval.ypred));
        
        % divide trials into k folds whose sizes are as equal as possible.
        % we interleave the folds to average out any nonstationarities in
        % the data
        n_rep = ceil(fit.n ./ fit.param.crossval.kfold);
        cvid = repmat((1:fit.param.crossval.kfold)',n_rep,1);
        cvid = cvid(1:fit.n);
        
        % save certain results for each fold
        fit.crossval.fold_id = cvid; % which cross-validation fold each trial was classified into
        fit.crossval.trainfit.b = nans(fit.nx,fit.param.crossval.kfold);
        fit.crossval.trainfit.se = nans(fit.nx,fit.param.crossval.kfold);
        fit.crossval.trainfit.p = nans(fit.nx,fit.param.crossval.kfold);
        fit.crossval.trainfit.warning.msg = cell(fit.nx,1);
        fit.crossval.trainfit.warning.id = cell(fit.nx,1);
        fit.crossval.trainfit.warning.detected = nans(fit.nx,1);
        
        % for each fold...
        for cvi = 1:fit.param.crossval.kfold
            train_okx = cvid ~= cvi;
            test_okx = cvid == cvi;
            

            % track initial state of warnings
            [lastmsg_before_fit,lastid_before_fit] = lastwarn;
            lastwarn('');

            % train model on all but the current test set
            [train_b,~,train_stats] = eglm_glmfit(fit.glmfun_inputs.x(train_okx,:),fit.y(train_okx),fit.param.distribution,'link',fit.param.link,'constant','off','offset',fit.glmfun_inputs.offset(train_okx));

            % check if new warnings were given
            [fit.crossval.trainfit.warning.msg{cvi},fit.crossval.trainfit.warning.id{cvi}] = lastwarn;
            fit.crossval.trainfit.warning.detected(cvi) = ~isempty(fit.crossval.trainfit.warning.msg{cvi}) || ~isempty(fit.crossval.trainfit.warning.id{cvi});
%             if fit.crossval.trainfit.warning.detected(cvi)
%                 fprintf(' eglm_fit detected at least one warning when calling glmfit to fit the GLM! Be wary of trusting the fitting results!\n');
%                 fprintf('  last warning msg: %s\n  last warning  id: %s\n',fit.warning.msg,fit.warning.id);
%             end
            % if no warning was detected, restore prior state of Matlab's warning flags
            if ~fit.crossval.trainfit.warning.detected(cvi)
                lastwarn(lastmsg_before_fit,lastid_before_fit);
            end

            
            train_b_for_pred = train_b;
            % TEST
            %train_b_for_pred = train_b;
            %train_b_for_pred(train_stats.p > 0.05) = 0;
            
            fit.crossval.trainfit.b(:,cvi) = train_b;
            fit.crossval.trainfit.se(:,cvi) = train_stats.se;
            fit.crossval.trainfit.p(:,cvi) = train_stats.p;
            
            % get its prediction for the test set
            test_ypred = glmval(train_b_for_pred,fit.glmfun_inputs.x(test_okx,:),fit.param.link,'constant','off','offset',fit.glmfun_inputs.offset(test_okx));
            fit.crossval.ypred(test_okx) = test_ypred;
            fit.crossval.sfit(test_okx) = train_stats.sfit; % save fitted s parameter (needed to calc log likelihood for some distributions)
        end
        fit.crossval = calc_model_predictions(fit.crossval,fit,fit.crossval.ypred,fit.crossval.sfit);
    end
    
    if fit.param.shuffle.n > 0
        fit.shuffle = struct();
        fit.shuffle.varnames = {'b','se','p','ypred','loglik_tr','loglik','mean_loglik_per_trial','r','r2','p_correct_prediction'};
        for vni = 1:numel(fit.shuffle.varnames)
            vn = fit.shuffle.varnames{vni};
            if ismember(vn,{'se','p'})
                curnvar = size(fit.b,1);
            else
                curnvar = size(fit.(vn),1);
            end
            fit.shuffle.(vn) = nans(curnvar,fit.param.shuffle.n);
        end
        
        % if using pre-specified shuffle, check that it has the right
        % dimensions and properties
        if ~isempty(fit.param.shuffle.use_shuffle_ids)
            assert(isequaln(size(fit.param.shuffle.use_shuffle_ids),[fit.n fit.param.shuffle.n]),'user specified shuffling matrix has wrong dimensions');
            sorted_shuf = sort(fit.param.shuffle.use_shuffle_ids);
            sorted_shuf_is_ok = all(sorted_shuf == (1:fit.n)',1);
            assert(all(sorted_shuf_is_ok),'user-specified shufflings must have the elements of 1:fit.n in each column');
        end
        
        % save the specific shuffle ids that we use, for future reference
        fit.shuffle.shufids = nans(fit.n,fit.param.shuffle.n);
        
        for shufi = 1:fit.param.shuffle.n
            % run fit on shuffled data
            
            if ~isempty(fit.param.shuffle.use_shuffle_ids)
                % use pre-specified shuffle?
                shufids = fit.param.shuffle.use_shuffle_ids(:,shufi);
            else
                % generate new shuffle
                shufids = randperm(fit.n)';
            end
            
            % save shuffle IDs we are using
            fit.shuffle.shufids(:,shufi) = shufids;
            
            shuf_y = fit.y(shufids);
            
            
            % check if warnings will be generated
            [lastmsg_before_fit,lastid_before_fit] = lastwarn;
            lastwarn('');
            
            % do the fit
            [shuf_b,~,shuf_stats] = eglm_glmfit(fit.glmfun_inputs.x,shuf_y,fit.param.distribution,'link',fit.param.link,'constant','off','offset',fit.glmfun_inputs.offset);
            
            % check if new warnings were given
            [fit.shuffle.warning.msg{shufi},fit.shuffle.warning.id{shufi}] = lastwarn;
            fit.shuffle.warning.detected(shufi) = ~isempty(fit.shuffle.warning.msg{shufi}) || ~isempty(fit.shuffle.warning.id{shufi});
            % if no warning was detected, restore prior state of Matlab's warning flags
            if ~fit.shuffle.warning.detected(shufi)
                lastwarn(lastmsg_before_fit,lastid_before_fit);
            end
            
            shuf_ypred = glmval(shuf_b,fit.glmfun_inputs.x,fit.param.link,'constant','off','offset',fit.glmfun_inputs.offset);
            fit.shuffle.b(:,shufi) = shuf_b;
            fit.shuffle.se(:,shufi) = shuf_stats.se;
            fit.shuffle.p(:,shufi) = shuf_stats.p;
            fit.shuffle.ypred(:,shufi) = shuf_ypred;
            fit.shuffle = calc_model_predictions(fit.shuffle,fit,shuf_ypred,shuf_stats.sfit,shuf_y,[shufi fit.param.shuffle.n]);
        end
    end
end

% if required, calculate per-offer effects of subsets of regressors on the
% decision variable
fit.peroffer = struct();
fnames = fieldnames(fit.param.peroffer);
for fni = 1:numel(fnames)
    fn = fnames{fni};
    
    % which xregressors to include in computing this effect?
    okx = false(fit.nx,1);
    
    if ischar(fit.param.peroffer.(fn))
        assert(strcmp(fit.param.peroffer.(fn),'*'));
        % find all xregressors with settings specified for the max # offers
        % (i.e., if all xreg are specified for just 1 offer, we assume this
        % is a regression with just 1 offer, so we use all of them; but if
        % some are specified for 2, and others for 1, we only use the
        % regressors specified for both offers)
        x_noffers = nans(1,fit.nx);
        for xi = 1:fit.nx
            x_noffers(xi) = size(fit.glmfun_inputs.x_per_offer{xi},2);
        end
        max_noffers = max(x_noffers);
        for xi = 1:fit.nx
            okx(xi) = size(fit.glmfun_inputs.x_per_offer{xi},2) == max_noffers;
        end
        
    elseif iscell(fit.param.peroffer.(fn))
        for vi = 1:numel(fit.param.peroffer.(fn))
            curxname = fit.param.peroffer.(fn){vi};
            xi = find(strcmp(curxname,fit.xname));
            % for now, allow to use multiple regressors with same name,
            % whose effects will just be added together
            %assert(numel(xi)==1,'did not find one unique xregressor named "%s" in list of xnames!',curxname);
            okx(xi) = true;
        end
        
    else
        error('unknown per-offer effect specification');
    end
    
    % save which xregressors we are using for this effect computation
    fit.param.peroffer.xreg_used.(fn) = okx;
    fit.param.peroffer.nxreg_used.(fn) = sum(okx);
    fit.param.peroffer.noffers.(fn) = nan;
    
    if sum(okx) < 1
        % no xreg specified, so their effect is NaN
        fit.peroffer.(fn) = nans(fit.n,2);
    else
        % at least one xreg was specified, so compute effect
        % by summing up effects of each xreg
        fit.peroffer.(fn) = zeros(fit.n,2);
        x_noffers = nans(1,fit.nx);
        for xi = 1:fit.nx
            if okx(xi)
                x_noffers(xi) = size(fit.glmfun_inputs.x_per_offer{xi},2);
                for oi = 1:x_noffers(xi)
                    fit.peroffer.(fn)(:,oi) = fit.peroffer.(fn)(:,oi) + fit.b(xi)*fit.glmfun_inputs.x_per_offer{xi}(:,oi);
                end
            end
        end
        % if only N offers per trial were in the data, only compute that
        % many per-offer effects
        max_noffers = max(x_noffers(okx));
        if isempty(max_noffers)
            max_noffers = 0;
        end
        assert(inbounds(max_noffers,[1 2]),'expected 1 or 2 offers per trial, when computing per-offer effects');
        fit.param.peroffer.noffers.(fn) = max_noffers;
        fit.peroffer.(fn) = fit.peroffer.(fn)(:,1:max_noffers);
    end
    
end

% if required, rescale beta weights and their SEs (but no other stats)
% to be in units based on a fixed change in an original raw variable
% (e.g. "relative to the effect of adding +1 mL of juice")
if ~ischar(fit.param.betascale.name)
    % no scaling. Set the name of the beta weight units based on the name
    % of the y variable and its transformation by the link function
    curyname = fit.yname;
    if strcmp(fit.param.distribution,'binomial')
        if strcmp(fit.param.link,'logit')
            curyname = ['log odds of ' curyname];
        else
            curyname = ['p(' curyname ')'];
        end
    end
    
    fit.param.betascale.name = ['Effect on ' curyname];
else
    % scaling needs to be done
    xi = find(strcmp(fit.xname,fit.param.betascale.xreg_name));
    assert(numel(xi)==1,'did not find unique regressor with name "%s" to do rescaling of beta weights',fit.param.betascale.name);
    
    % get effect of adding the specified change to the regressor
    normparam = fit.param.normalization{xi};
    nterms = numel(normparam.sd);
    if any(normparam.is_normalized)
        assert(nterms==1,'cannot rescale beta weights using regressor "%s", can only rescale beta weights based on a non-normalized regressor or a regressor that was computed using a single term (not the product of multiple terms)',fit.param.betascale.name);

        % get SD used to normalize the regressor
        nsd = normparam.sd;
        % if was not normalized, we treat it as if normalization SD = 1
        % (i.e., no normalization by SD)
        if isnan(normparam.sd)
            nsd = 1;
        end

        % fitted beta (in units of change in predicted variable per addition 
        %  of +1 to the normalized variable)
        b_normalized = fit.b(xi);
        % fitted beta (in units of the change in predicted variable per
        %  addition of +1 to the ORIGINAL variable)
        b_orig = b_normalized ./ nsd;
    else
        % get effect of adding the specified change to the regressor,
        % (as a change in the product of the raw terms)
        b_orig = fit.b(xi);
    end
    
    % rescaling factor (i.e. factor to multiply all betas so that are in
    %  units of the equivalent change in predicted variable per addition of
    %  +xreg_change to the ORIGINAL variable)
    fit.param.betascale.rescaling_factor = 1 ./ (b_orig .* fit.param.betascale.xreg_change);
    
    % apply the rescaling
    fit.b = fit.b .* fit.param.betascale.rescaling_factor;
    fit.stats.beta = fit.stats.beta .* fit.param.betascale.rescaling_factor;
    fit.stats.se = fit.stats.se .* fit.param.betascale.rescaling_factor;
end

end

function [pfit] = calc_model_predictions(pfit,fit,ypred,sfit,y,fit_index_of_nfits)
    % pfit = data structure to which the model prediction results should be
    %        added as new fields
    % fit = the true "fit" data structure
    % ypred = predicted values of y
    % sfit = fitted SD used for the noise distribution of normal models
    % y = the y variable used when fitting the model
    %         (which may be different from the original y-variable, e.g.
    %         when fitting a shuffled dataset)
    % fit_index_of_nfits = [current_fit_index number_of_total_fits]
    %         used to say which fit (of possibly multiple fits) to store 
    %         the results in inside of pfit. E.g. store the single-trial 
    %         log likelihoods in pfit.loglik_tr, that matrix is of size
    %         ntrials x number_of_total_fits, and the results of this fit
    %         is stored in: in pfit.loglik_tr(:,current_fit_index)
    %
    % if calculating predictions of original fit, pfit=fit, so just call
    %  calc_model_predictions(fit)
    % if calculating cross-validated predictions, call with 
    %  calc_model_predictions(fit.crossval,fit,fit.crossval.ypred,fit.crossval.sfit)
    if nargin < 2 || isempty(fit)
        fit = pfit;
    end
    if nargin < 3 || isempty(ypred)
        ypred = fit.ypred;
    end
    if nargin < 4 || isempty(sfit)
        sfit = fit.stats.sfit;
    end
    if nargin < 5 || isempty(y)
        y = fit.y;
    end
    if nargin < 6 || isempty(fit_index_of_nfits)
        fit_index_of_nfits = [1 1];
    end
    pfiti = fit_index_of_nfits(1);
    nfits = fit_index_of_nfits(2);
    assert(pfiti <= nfits,'cannot save data for fit %d of %d fits, the first number must be smaller than the second!');
    
    % calculate log likelihood of each individual trial
    % based on the model's fitted parameters and its assumption about
    % the distribution the data follow
    pfit.loglik_tr(:,pfiti) = nans(size(ypred,1),1);
    switch fit.param.distribution
        case 'binomial'
            pfit.loglik_tr(:,pfiti) = log((y==1).*ypred + (y==0).*(1-ypred));
        case 'poisson'
            pfit.loglik_tr(:,pfiti) = log(poisspdf(y,ypred));
        case 'normal'
            pfit.loglik_tr(:,pfiti) = log(normpdf(y,ypred,sfit));
        otherwise
            error('have not implemented calculation of log likelihood for GLM distribution "%s"',fit.param.distribution);
    end

    % total log likelihood of entire dataset
    pfit.loglik(1,pfiti) = sum(pfit.loglik_tr(:,pfiti));

    % mean log likelihood per trial
    pfit.mean_loglik_per_trial(1,pfiti) = pfit.loglik(1,pfiti) ./ fit.n;
    
    % correlation between predictions and data
    pfit.r(1,pfiti) = corr(ypred,y,'type','Pearson');
    % r-squared
    cur_sse_total = sum((y - mean(y)).^2);
    resid = y - ypred;
    cur_sse_resid = sum(resid.^2);
    pfit.r2(1,pfiti) = 1 - (cur_sse_resid ./ cur_sse_total);
    
    % if binomial data, calculate p(correct prediction)
    if strcmp(fit.param.distribution,'binomial')
        pfit.p_correct_prediction(1,pfiti) = mean(sign(ypred-0.5) == sign(y-0.5));
    else
        pfit.p_correct_prediction(1,pfiti) = nan;
    end
    
    % calculate AIC and BIC
    pfit.n_free_parameters(1,pfiti) = sum(~fit.b_is_manually_specified);
    pfit.aic(1,pfiti) = 2*pfit.n_free_parameters(1,pfiti) - 2*pfit.loglik(1,pfiti);
    pfit.bic(1,pfiti) = pfit.n_free_parameters(1,pfiti)*log(fit.n) - 2*pfit.loglik(1,pfiti);
end