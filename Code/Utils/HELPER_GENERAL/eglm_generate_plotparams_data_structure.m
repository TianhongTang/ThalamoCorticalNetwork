function pps = eglm_generate_plotparams_data_structure(fit,plotparams,default_color_mode)
% pps = eglm_generate_plotparams_data_structure(fit,plotparams,default_color_mode)
% 
% internal function used by eglm_plot functions to convert the
% user-specified 'plotparams' cell array into a data structure that can be
% used to easily specify plotting parameters

if isempty(plotparams)
    % plot all beta weights
    plotparams = fit.xname;
end
if iscell(plotparams) && any(cellfun(@ischar,plotparams))
    % specified a simple list of regressor names, so convert them into our 
    % more detailed format, so we can easily handle it with the code below
    assert(all(cellfun(@ischar,plotparams)),'plotparams must be a simple cell array of regressor names, or else in the more complex format (cell array of cell arrays)');
    plotparams = cellfun(@(z) {z},plotparams,'uniform',0);
    if strcmp(default_color_mode,'fixed')
        default_color = 'k';
    else
    end
else
    % specified params in our detailed format, so presumably we'll be
    % specifying the colors of some variables, so set the default color to
    % be less flashy (gray)
    default_color = [1 1 1]*.5;
end

switch default_color_mode
    case 'fixed'
        % already selected the default color
    case 'rotating'
        % use a rotating selection of colors
        default_color_set = {'r',interpcolor('g','k',.3),'b','c',interpcolor('m','k',.25),[1 .5 0],[1 1 1]*.5};
        default_color_set = [ ...
            default_color_set, ... % original colors
            cellfun(@(z) interpcolor(z,'w'),default_color_set,'uniform',0), ... % lighter versions of colors
            cellfun(@(z) interpcolor(z,'k'),default_color_set,'uniform',0), ... % darker versions of colors
            ];
        % convert them to RGB format
        default_color_set = cellfun(@colorspec_to_rgb,default_color_set,'uniform',0);
        
        
        % find all user-specified colors, and exclude all default colors
        % that are very close to them
        user_specified_colors = {};
        for pxi = 1:numel(plotparams)
            vi = 2;
            while vi <= numel(plotparams{pxi})
                if ischar(plotparams{pxi}{vi})
                    if ismember(plotparams{pxi}{vi},{'color','edgecolor','facecolor'}) && vi+1 <= numel(plotparams{pxi})
                        % get user-specified color and convert to RGB format
                        user_specified_colors{end+1} = colorspec_to_rgb(plotparams{pxi}{vi+1});
                    end
                    vi = vi + 2;
                else
                    error('expected plotparams{%d}{%d} to be a string',pxi,vi);
                    vi = vi + 1;
                end
            end
        end
        min_color_distance = inf*ones(size(default_color_set));
        for uci = 1:numel(user_specified_colors)
            for dci = 1:numel(default_color_set)
                cur_color_distance = sqrt(sum((user_specified_colors{uci} - default_color_set{dci}).^2));
                min_color_distance(dci) = min(min_color_distance(dci),cur_color_distance);
            end
        end
        % remove all default colors that are very close to a user-specified color
        ok_default_colors = min_color_distance >= 0.15;
        default_color_set = default_color_set(ok_default_colors);
        
        % if all have been excluded, force a single arbitrary default color to be used
        if isempty(default_color_set)
            default_color_set = {'k'};
        end
        
        % start with #1, continue until the end, then rotate colors
        default_color = { ...
            default_color_set, 
            1, ...
            };
    otherwise
        error();
end

% validate plot parameters and store in a data structure
pps = cell(numel(plotparams),1);
for pxi = 1:numel(plotparams)
    curpp = plotparams{pxi};
    assert(~isempty(curpp) && ischar(curpp{1}));
    xname = curpp{1};
    xi = find(strcmp(fit.xname,xname));
    if numel(xi) ~= 1
        fprintf('warning: did not find a unique regressor with name "%s" from %d-th entry of plotparams, skipping\n',xname,pxi);
        continue;
    end
    
    pps{pxi} = struct('name',xname,'xi',xi, ...
        'color',[],'facecolor',[],'edgecolor',[],'facealpha',1,'markerfacecolor',[],'markeredgecolor',[],'marker','o','linestyle','-','linewidth',1,'fontweight',[], ...
        'x',[],'group','default');
    
    vi = 2;
    while vi <= numel(curpp)
        assert(ischar(curpp{vi}),'expected entry %d of plotparams{%d} to be a string!',vi,pxi);
        assert(vi+1 <= numel(curpp),'expected entry %d of plotparams{%d} to be followed by a value for that parameter!',vi,pxi);
        % set the specified parameter to its specified value
        pps{pxi}.(curpp{vi}) = curpp{vi+1};
        vi = vi + 2;
    end
    % apply defaults to un-specified parameters
    
    % if main 'color' is not already specified, try picking a color by
    % looking in these variables one by one:
    %  edgecolor, facecolor, markerfacecolor, markeredgecolor
    if isempty(pps{pxi}.color)
        candidate_colors = { pps{pxi}.edgecolor, pps{pxi}.facecolor, pps{pxi}.markerfacecolor, pps{pxi}.markeredgecolor };
        for cci = 1:numel(candidate_colors)
            if ~isempty(candidate_colors{cci}) && ~ismember(candidate_colors,{'auto','none'})
                pps{pxi}.color = candidate_colors{cci};
                break;
            end
        end
    end
    % if still no 'color' has been found, use the next default color
    if isempty(pps{pxi}.color)
        % apply the current default color, and if we are using a rotating
        % default color roster, move to the next default color
        if iscell(default_color)
            curcolor = default_color{1}{default_color{2}};
            default_color{2} = default_color{2} + 1;
            if default_color{2} > numel(default_color{1})
                default_color{2} = 1;
            end
        else
            curcolor = default_color;
        end
        
        pps{pxi}.color = curcolor;
    end
	
    % default edge color = same as main color
    if isempty(pps{pxi}.edgecolor)
        pps{pxi}.edgecolor = pps{pxi}.color;
    end
    % default face color = faded version of its color
    if isempty(pps{pxi}.facecolor)
        pps{pxi}.facecolor = interpcolor(pps{pxi}.color,'w',.85);
    end
    % default marker edge color = 'auto'
    if isempty(pps{pxi}.markeredgecolor)
        pps{pxi}.markeredgecolor = 'auto';
    end
    % default marker face color = 'none'
    if isempty(pps{pxi}.markerfacecolor)
        pps{pxi}.markerfacecolor = 'none';
    end
        
    % default font weight = bold if it is a thick line
    if isempty(pps{pxi}.fontweight)
        if pps{pxi}.linewidth > 1
            pps{pxi}.fontweight = 'bold';
        else
            pps{pxi}.fontweight = 'normal';
        end
    end
end

% remove plotparams for beta weights that we didn't find in the data
pps = pps(cellfun(@(z) ~isempty(z),pps));

% set default x coordinates for each plotted variable that wasn't already
% specified
xnext = 1;
for pxi = 1:numel(pps)
    if isempty(pps{pxi}.x)
        % x coord not specified, so set to our default
        pps{pxi}.x = xnext;
        % set next variable's x coordinate to come after this variable's
        % x coordinate
        xnext = xnext + 1;
    else
        % set next variable's x coordinate to come after all x coords
        % that have been specified thus far
        xnext = max(xnext,pps{pxi}.x + 1);
    end
end

% set default 'groups' for plotting
ppgroup = cellfun(@(z) z.group,pps,'uniform',0);
pp_with_default_group = strcmp(ppgroup,'default');

non_default_groupid = cellfun(@(z) z.group,pps(~pp_with_default_group));
% start with a unique group ID that is higher than all the user-specified
% groups
next_default_group = max(non_default_groupid) + 1;
if isempty(next_default_group)
    next_default_group = 1;
end

% by default, put each successive plotted beta weight in the same group
% if it has the same plot parameters (except for its name, regressor index, 
% x-coordinate, etc.) as the previous beta weight; otherwise create a new
% unique group.
prev_pxi = [];
for pxi = 1:numel(pps)
    if pp_with_default_group(pxi)
        cur = pps{pxi};
        % if we previously assigned a beta weight to a default group...
        if ~isempty(prev_pxi)
            % was it basically the same plot params as the current beta
            % weight?
            prev = pps{prev_pxi};
            fnames = setdiff(fieldnames(cur),{'name','xi','x','group'});
            is_identical = true;
            for fni = 1:numel(fnames)
                fn = fnames{fni};
                is_identical = is_identical & isequaln(cur.(fn),prev.(fn));
                if ~is_identical
                    break;
                end
            end
            
            if ~is_identical
                next_default_group = next_default_group + 1;
            end
        end
        pps{pxi}.group = next_default_group;
        prev_pxi = pxi;
    end
end
