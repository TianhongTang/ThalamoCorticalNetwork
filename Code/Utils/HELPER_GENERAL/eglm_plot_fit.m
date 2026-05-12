function eglm_plot_fit(fit,varargin)
% eglm_plot_fit(fit,...)
%
% basic plot of the beta weights from a fit, their standard errors, and
% their significance.
%
% required inputs are:
%  fit - the data structure resulting from running eglm_fit
%
% optional inputs are:
%
% 'axes',h - plots to the specified axes (default: current axes)
%
% 'params',plotparams - specifies a cell array instructing to plot a 
%  subset of beta weights, and (optionally) to use specific fomatting 
%  parameters when doing so (e.g. colors), and plotting them at specific 
%  locations on the x-axis. Define plotparams as a simple list of regressor
%  names, e.g.:
% 
%   {'regressor_name_1','regressor_name_2',...}
%
%  or as a more complex list of cell arrays, each starting with a regressor
%  name and then optionally specifying various plotting details, e.g.:
%  
%   { ...
%     { 'regressor_name', 'color', c, 'edgecolor', ec, 'facecolor', fc,
%       'facealpha', fa, 'marker', m, 'linestyle', ls,
%       'linewidth', lw, 'fontweight', fw, ...
%       'x', x, 'group', g}, ...
%     ...
%   }
%
%   ...where the only required parameter for a given beta weight is
%   regressor_name, and all the others are optional. Most of these specify
%   arguments to plotting functions. 'x' specifies the x-coordinate where
%   it should be plotted. 'group' specifies an arbitrary group ID; all
%   betas with the same group ID are plotted with the same parameters
%   (set by the first listed beta in the group), and in a line plot, are
%   plotted with a connected line. The default is for all beta weights to
%   have x coordinates in a sequence 1, 2, 3, ..., and to all be in the
%   same group.
%
% 'pvalmode',pvalmode - specifies the manner of plotting p-values. Give
%  a 3-element cell array with text elements specifying how to handle 
%  p-values that are significant, marginal, or other. Allowed 
%  specifications are 'star', 'text', and 'none'. If pvalmode is set to [], 
%  it uses the default settings:  {'star','text','none'}, indicating to 
%  plot asterisks for significant p-values, text for marginal p-values 
%  (0.05 <= p < 0.15), and nothing for other p-values.
%
% 'plotmode',plotmode - specifies the type of plot. Valid settings are 
%  'bar' (the default) and 'line' (useful if using 'params' to specify 
%  certain beta weights to be plotted at nearby x-coordinates, since in a 
%  bar plot the bars would be overlapping).
%
% 'barwidth',barwidth, specifies the width of the bars on a bar plot
% (default: 1)
%

%% handle variable input parameters
hax = [];
plotparams = {};
pvalmode = {'star','text','none'};
plotmode = 'bar';
barwidth = 1;

vi = 1;
while vi <= numel(varargin)
    assert(ischar(varargin{vi}),'expected next variable argument to be a string!');
    switch varargin{vi}
        case 'axes'
            hax = varargin{vi+1};
            vi = vi + 1;
        case 'params'
            plotparams = varargin{vi+1};
            vi = vi + 1;
        case 'pvalmode'
            pvalmode = varargin{vi+1};
            vi = vi + 1;
        case 'plotmode'
            plotmode = varargin{vi+1};
            vi = vi + 1;
        case 'barwidth'
            barwidth = varargin{vi+1};
            vi = vi + 1;
        otherwise
            error('unknown variable argument!');
    end
    vi = vi + 1;
end

% if parameters are not specified, apply defaults

% horrible hack to make a new figure+axes manually if one does not exist,
% to ensure that any such ones are created with my 'figuren' command and
% not the builtin 'figure', to ensure that 'hold' is 'on' so that I can
% plot multiple things without them deleting each other
% 
% note that we can't use 'gcf' for this since if no figure exists it will
% by default create one using 'figure'.
if isempty(get(groot,'CurrentFigure'))
    assert(isempty(hax),'you specified the plot to go onto a specified axes, but cannot find any current figure where that axes might reside...');
    figuren;
end

% horrible hack to plot to axes other than the current axes, by setting
% the specified axes to be current, and then un-setting it when this
% function is finished
% (this is slow because calls to axes can be very slow)
original_gca = gca;
if ~isempty(hax)
    axes(hax);
end

% generate easily usable data structure holding the information from the
% 'plotparams' argument telling us how to plot the data, using a fixed
% default color for variables whose colors are not user-specified
pps = eglm_generate_plotparams_data_structure(fit,plotparams,'fixed');

% keep track of the appropriate xticks and their labels as we generate the
% plots
raw_xtick = [];
raw_xticklabel = {};

% get groups of betas to plot
ppgroup = cellfun(@(z) z.group,pps);
ugroups = unique(ppgroup);

%% plot each group of beta weights
for ugi = 1:numel(ugroups)
    gi = ugroups(ugi);
    
    % get the individual beta weights in this group, and get their plotting
    % variables (x, y, SE, p-value, etc.)
    group_pxi = find(ppgroup==gi);
    
    x = nans(numel(group_pxi),1);
    y = nans(numel(group_pxi),1);
    yse = nans(numel(group_pxi),1);
    p = nans(numel(group_pxi),1);
    ptext = cell(numel(group_pxi),1);
    ptextparams = cell(numel(group_pxi),1);
    
    % for each individual beta weight in the group...
    for gpxi = 1:numel(group_pxi)
        pxi = group_pxi(gpxi);
        cur = pps{pxi};
        xi = cur.xi;
        
        % save reference to plotting params for the first beta weight in 
        % this group (we will use these params for all other members in 
        % the group, except for e.g. the name and x coordinate)
        if gpxi == 1
            grp = cur;
            switch grp.fontweight
                case 'normal'
                    fontweightstring = '\rm';
                case 'bold'
                    fontweightstring = '\bf';
                otherwise
                    error('unknown setting of grp.fontweight');
            end
        end
        
        % get its plot data
        x(gpxi) = cur.x;
        y(gpxi) = fit.stats.beta(xi);
        yse(gpxi) = fit.stats.se(xi);
        p(gpxi) = fit.stats.p(xi);
        
        % get xtick and xticklabel settings
        % do not overwrite xticks/labels that have already been set
        if ~any(raw_xtick == x(gpxi))
            raw_xtick(end+1) = x(gpxi);
            raw_xticklabel{end+1} = [texcolor(grp.color) fontweightstring cur.name];
        end
        
        % setup the text to indicate the p-value
        if p(gpxi) < 0.05
            ptype = 1;
        elseif 0.05 <= p(gpxi) < 0.15
            ptype = 2;
        else
            ptype = 3;
        end
        asterisk_char = '*';
        switch pvalmode{ptype}
            case 'star'
                if p(gpxi) < 0.001
                    curtext = repmat(asterisk_char,1,3);
                elseif p(gpxi) < 0.01
                    curtext = repmat(asterisk_char,1,2);
                elseif p(gpxi) < 0.05
                    curtext = repmat(asterisk_char,1,1);
                else
                    curtext = '';
                end
                if isempty(curtext)
                    ptext{gpxi} = '';
                else
                    ptext{gpxi} = ['\fontsize{24}' texcolor(grp.color) curtext];
                end
            case 'text'
                ptext{gpxi} = ['\fontsize{8}' texcolor(grp.color) roundstr(.001,p(gpxi))];
            otherwise
                ptext{gpxi} = '';
        end
        if y(gpxi) < 0
            % plot p-value text below negative betas
            ptextparams{gpxi} = {'ct',x(gpxi),y(gpxi) - 1.25*yse(gpxi)};
        else
            % plot p-value text above positive betas
            ptextparams{gpxi} = {'cb',x(gpxi),y(gpxi) + 1.25*yse(gpxi)};
        end
    end
    
    % main plot! Plot the bars (or lines) for the beta weights
    switch plotmode
        case 'bar'
            plot_errorbar(x,y,yse, ...
                {'bar',grp.facecolor,'fixedwidth',barwidth,'facealpha',grp.facealpha,'edgecolor',grp.edgecolor,'linestyle',grp.linestyle,'linewidth',grp.linewidth}, ...
                {'color',grp.color,'linestyle','-','linewidth',grp.linewidth});
        case 'line'
            plot_errorbar(x,y,yse, ...
                {'color',grp.color,'linestyle',grp.linestyle,'linewidth',grp.linewidth,'marker',grp.marker,'markerfacecolor',grp.markerfacecolor,'markeredgecolor',grp.markeredgecolor}, ...
                {'color',grp.color,'linestyle','-','linewidth',grp.linewidth});
        otherwise
            error('unknown plotmode');
    end

    % plot the p-value text for each beta weight
    for gpxi = 1:numel(ptext)
        if ~isempty(ptext{gpxi})
            etext(ptextparams{gpxi}{:},ptext{gpxi});
        end
    end
end

%% handle basic plot formatting:

% set x limits
setlim(gca,'xlim','tight',.05);

% set y limits and plot a horizontal line indicating beta = 0
setlim(gca,'ylim','tight',.1);
liney(0,'k');

% apply x ticks and their labels
% there may be redundant or out-of-order x values, so only take the unique,
% sorted x ticks and labels
[xtick,ixtick] = unique(raw_xtick,'sorted');
xticklabel = raw_xticklabel(ixtick);
set(gca,'xtick',xtick,'xticklabel',xticklabel);
xtickangle(gca,45);

% set y-label based on the name of the units of the beta weights
ylabel(fit.param.betascale.name);

% plot basic stats of the fit
basictext = {};
if ~isempty(fit.name)
    basictext{end+1} = fit.name;
end
if numel(pps) ~= fit.nx
    basictext{end+1} = sprintf('showing %d/%d regressors',numel(pps),fit.nx);
end
if ~isnan(fit.p_correct_prediction)
    basictext{end+1} = ['p(correct)=' roundstr(1,100*fit.p_correct_prediction) '%'];
end
basictext{end+1} = ['lik/tr = ' roundstr(.0001,exp(fit.mean_loglik_per_trial))];
basictext{end+1} = ['n=' num2str(fit.n) ' trials'];
etextn('lb',basictext,'color',[1 1 1]*.5);

% horrible hack to plot to axes other than the current axes, by setting
% the specified axes to be current, and then un-setting it when this
% function is finished
% (this is slow because calls to axes can be very slow)
if ~(original_gca == hax)
    axes(original_gca);
end
