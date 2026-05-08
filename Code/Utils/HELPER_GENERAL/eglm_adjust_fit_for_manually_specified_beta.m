function fit = eglm_adjust_fit_for_manually_specified_beta(fit)
% fit = eglm_adjust_fit_for_manually_specified_beta(fit)
%
% internal function for adjusting a GLM fit data structure to reflect the
% betas whose values were manually specified rather than fitted by the GLM
% fitting procedure.

% which betas are manually specified?
assert(isequaln(fit.b_is_manually_specified,~isnan(fit.glmfun_inputs.manually_specified_beta)),'inconsistency in which betas were fitted vs. manually specified in the model');
okman = fit.b_is_manually_specified;
okfit = ~okman;

% get betas from the raw betas from calling glmfit
% which should be saved in fit.glmfun_inputs.b
% manually specified betas are set to their specified values
assert(isfield(fit.glmfun_inputs,'b'),'was the fit carried out properly? fit.glmfun_inputs.b must exist to set betas based on raw fitted betas and manually specified betas!');
assert(isequaln(size(fit.glmfun_inputs.b),[sum(okfit) 1]),'was the fit carried out properly? fit.glmfun_inputs.b should be an N x 1 vector where N is the number of betas fitted by glmfit!');
fit.b = nans(fit.nx,1);
fit.b(okman) = fit.glmfun_inputs.manually_specified_beta(okman);
fit.b(okfit) = fit.glmfun_inputs.b;

% if there are any manually-specified variables, adjust the 'stats'
% structure
if any(okman)
    % adjust beta
    fit.stats.beta = fit.b;

    % adjust SE, t statistics, and p-values
    % manually specified betas have them set to NaN
    orig = fit.stats.se;
    fit.stats.se = nans(fit.nx,1);
    fit.stats.se(okfit) = orig;

    orig = fit.stats.t;
    fit.stats.t = nans(fit.nx,1);
    fit.stats.t(okfit) = orig;

    orig = fit.stats.p;
    fit.stats.p = nans(fit.nx,1);
    fit.stats.p(okfit) = orig;

    % adjust cov and coeffcorr matrices

    % est. covariance matrix of the betas
    % set it to 0 for all beta pairs involving a manually-specified beta
    orig_cb = fit.stats.covb;
    fit.stats.covb = nans(fit.nx,fit.nx);
    fit.stats.covb(okfit,okfit) = orig_cb;
    fit.stats.covb(:,okman) = 0;
    fit.stats.covb(okman,:) = 0;
    
    % est. correlation of fitted beta coefficients
    % set it to 0 for all beta pairs involving a manually specified beta, 
    % EXCEPT for correlations between a coefficient and itself, which are
    % NaN because it's a correlation between two variables which have no
    % variance.
    orig_cc = fit.stats.coeffcorr;
    fit.stats.coeffcorr = nans(fit.nx,fit.nx);
    fit.stats.coeffcorr(okfit,okfit) = orig_cc;
    fit.stats.coeffcorr(:,okman) = 0;
    fit.stats.coeffcorr(okman,:) = 0;
    okmanid = find(okman);
    for mid = 1:numel(okmanid)
        fit.stats.coeffcorr(okmanid(mid),okmanid(mid)) = NaN;
    end

end