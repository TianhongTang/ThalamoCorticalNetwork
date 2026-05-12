function [x,y,normalization_params,glmfun_inputs] = eglm_make_regression_inputs(xreg,yreg,oktrials)
% [x,y,normalization_params,glmfun_inputs] = eglm_make_regression_inputs(xreg,yreg[,oktr])
% 
% internal function for eglm_fit to convert a specification of regression
% data in variables xreg and yreg to the actual input data to be
% given to the glm
%
% yreg directly specifies the to-be-predicted y variable
% xreg gives the raw materials used to compute the regressors. By default,
%  this function computes each x regressor by taking its component terms,
%  then (A) applies z-scoring normalization to each term, (B) multiplying
%  all of the terms with each other, (C) if terms were specified separately
%  for each of two 'offers' on each trial, converting them into a single
%  variable representing the difference of the regressor between the two 
%  offers.
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
%           {'mean and sd', M, S} - normalize as if the data had the 
%            specified mean and SD. M and S must be scalars.
%
%         xreg{xi}.manually_specified_beta: forces the regression weight to
%          be a specific value. That beta weight's SE and p-value will be
%          set to NaN. Use this when you have a set of models in the same
%          general framework, but some of them operate by forcing specific
%          parameters to be fixed weights rather than freely fit weights.
% 
%  yreg - either a column vector (Ntrials x 1) of the values to be
%         predicted, or a data structure with that stored in its
%         yreg.value field.
%
%  oktrials - optional argument, a (Ntrials x 1) logical vector or a 
%         (k x 1) numeric vector (where k <= Ntrials) specifying a subset 
%         of trials. The analysis will be restricted to those trials.
%
% outputs:
%  x, y - the variables x (regressors) and y (variable to be predicted)
%  normalization_params: cell array saying how regs. were normalized
%  glmfun_inputs: a data structure with certain actual inputs used to call
%   glm functions (like glmfit and glmval). These may be different from "x"
%   because they separate all the effects of regressors that have manually
%   specified weights into a separate 'offset' variable to be passed to the
%   glm functions. Key fields:
%
%    glmfun_inputs.x - Ntrials x Nactualregressors matrix of regressors to
%     be actually given to glmfit
%
%    glmfun_inputs.offset - Ntrials x 1 'offset' input including the
%     effects of the manually specified betas. If there are no manually
%     specified betas, this will be all 0s.
%
%    glmfun_inputs.manually_specified_beta - Noriginalregressors x 1
%     vector, with manually specified beta values for each regressor (or
%     NaN for typical regressors where beta is not manually specified and
%     hence will be fitted by the model)

% handle y
if isstruct(yreg)
    assert(isfield(yreg,'value'),'cannot find field yreg.value to get y variable for regression');
    y = yreg.value;
else
    y = yreg;
end
assert(iscolumn(y) && (isnumeric(y) || islogical(y)),'y data must be a column vector that is numeric or logical vector');

% if user specified it, only use a subset of trials (the 'ok' trials)
% for analysis. Otherwise, use all trials
if nargin < 3 || strcmp(oktrials,'default')
    oktrials = true(size(y));
end
y = y(oktrials,:);


% handle x
x = nans(size(y,1),numel(xreg));
normalization_params = cell(1,numel(xreg));
glmfun_inputs = struct('x',[],'offset',zeros(size(y,1),1),'manually_specified_beta',nans(numel(xreg),1));

x_per_offer = cell(1,numel(xreg)); % for saving & returning per-offer inputs to GLM fitting

for xi = 1:numel(xreg)
    cur = xreg{xi};
    assert(isfield(cur,'terms'),'xreg{%d}.terms must exist',xi);
    
    % get a raw list of the terms for this x regressor,
    % and which of them are 'flag' variables (i.e. the ones that come after
    % the term entry that is labeled 'flag'
    flagid = find(strcmp(cur.terms,'flag'));
    assert(numel(flagid) <= 1,'the string "flag" should only appear once in specification of xregressors, before the terms that are flag variables');
    if isempty(flagid)
        xorig = cur.terms;
        is_flag = false(1,numel(xorig));
    else
        xorig = cur.terms([(1:(flagid-1)) ((flagid+1):end)]);
        is_flag = (1:numel(xorig)) >= flagid;
    end
    nterms = numel(xorig);
    
    % select only the 'ok' trials for analysis
    xraw = cell(size(xorig));
    for ti = 1:nterms
        assert(numel(size(xorig{ti})) == 2 && (iscolumn(xorig{ti}) || ismember(size(xorig{ti},2),[1 2])),'expected all terms to have one or two columns!');
        xraw{ti} = xorig{ti}(oktrials,:);
    end
    
    %% if requested, normalize all usual variables (but not flag variables)
    % using a zscore-like normalization
    normparam = struct( ...
        'is_flag',is_flag, ...
        'is_normalized',false(1,nterms), ...
        'mean',nans(1,nterms), ...
        'sd',nans(1,nterms));
    
    if ~isfield(cur,'normalization')
        cur.normalization = 'default';
    end
    
    % normalization parameters can be specified by cell array or by a
    % string saying a commonly used method
    if ischar(cur.normalization)
        switch cur.normalization
            case {'default','only interaction terms'}
                % default: normalize each term based on its own mean
                % and SD
                normparam.is_normalized = ~is_flag;
                % if only normalize interaction terms, do not normalize the
                % first term
                if strcmp(cur.normalization,'only interaction terms')
                    normparam.is_normalized(1) = false;
                end
                for ti = 1:nterms
                    if normparam.is_normalized(ti)
                        normparam.mean(ti) = mean(xraw{ti}(:));
                        normparam.sd(ti) = std(xraw{ti}(:));
                    end
                end
            case 'none'
                % none: do not normalize anything
                normparam.is_normalized(:) = false;
            otherwise
                error('unknown normalization style, expected "default" or "none"');
        end
    else
        assert(iscell(cur.normalization),'field "normalization" must be a string or cell array');
        assert(ischar(cur.normalization{1}),'expected field "normalization", if a cell array, to have its first argument be a string indicating the type of normalization');
        
        switch cur.normalization{1}
            case 'oktrials'
                % 'oktrials': same as default, but normalize based on the data 
                % in a specific subset of trials of the original dataset 
                % (which can be different from the subset that was selected 
                %  for analysis by the 'oktrials' argument to the function)
                assert(numel(cur.normalization)==2);
                norm_oktrials = cur.normalization{2};

                normparam.is_normalized = ~is_flag;
                for ti = 1:nterms
                    if normparam.is_normalized(ti)
                        cur_xraw_to_use_to_normalize = xorig{ti}(norm_oktrials,:);
                        normparam.mean(ti) = mean(cur_xraw_to_use_to_normalize(:));
                        normparam.sd(ti) = std(cur_xraw_to_use_to_normalize(:));
                    end
                end
            case 'mean and sd'
                % 'mean and sd': normalize as if each term in the data 
                % had the specified mean and sd.
                assert(numel(cur.normalization)==3);
                M = cur.normalization{2}; % mean
                S = cur.normalization{3}; % sd
                assert(numel(M) == nterms && numel(S) == nterms,'cannot normalize regressor %d (which has %d terms) using specified mean and sd (which have %d and %d terms), number of terms do not match',xi,nterms,numel(M),numel(S));
                
                normparam.is_normalized = ~is_flag;
                for ti = 1:nterms
                    if normparam.is_normalized(ti)
                        assert(~isnan(M(ti)) && ~isnan(S(ti)),'cannot normalize with regressor %d term %d with NaN specified as the mean (%f) or sd (%f)!',xi,ti,M(ti),S(ti));
                        normparam.mean(ti) = M(ti);
                        normparam.sd(ti) = S(ti);
                    end
                end
                
            otherwise
                error('unknown normalization type');
        end

    end
    
    % apply the normalization
    for ti = 1:nterms
        if normparam.is_normalized(ti)
            xraw{ti} = (xraw{ti} - normparam.mean(ti)) ./ normparam.sd(ti);
        end
    end
    
    % construct the regressor as the product of its terms
    xcur = xraw{1};
    for ti = 2:nterms
        assert(isequal(size(xcur),size(xraw{ti})),'expected the matrices for all terms to be the same size!');
        xcur = xcur .* xraw{ti};
    end
    
    % if the regressor was specified separately for each offer, take the
    % difference between its values for the two offers
    x_per_offer{xi} = xcur;
    noffers = size(xcur,2);
    switch noffers
        case 1
            % xcur is alreday correct, do nothing
        case 2
            % xcur has separate data for two offers, convert it to a
            % differential regressor reflecting the difference in this
            % regressor between the two offers
            xcur = xcur(:,2) - xcur(:,1);
        otherwise
            error('expected xreg{%d} to have a single column or two columns reflecting its values for two separate offers',xi);
    end
    
    % save the results
    x(:,xi) = xcur;
    normalization_params{xi} = normparam;
    
    % handle regressors with manually-specified betas
    if isfield(cur,'manually_specified_beta')
        % record what it's specified beta is
        glmfun_inputs.manually_specified_beta(xi) = cur.manually_specified_beta;
        % apply its constant effect to the GLM's predictors, by adding it
        % to the special 'offset' predictor variable whose coefficient is
        % fixed at 1.0.
        glmfun_inputs.offset = glmfun_inputs.offset  + x(:,xi) .* cur.manually_specified_beta;
    end
end

% the 'x' used as input to fitting the GLM only consists of the regressors
% whose betas are NOT already forced to be a specific value by manual
% specifications
glmfun_inputs.x = x(:,isnan(glmfun_inputs.manually_specified_beta));
glmfun_inputs.x_per_offer = x_per_offer;
