function varargout = dgz_explore(file)
% DGZ_EXPLORE  A tool for interactively exploring the data in DGZ files
% Usage:
%     dgz_explore
%     dgz_explore(filename)
%     dgz_explore(directory)
%     dgz_explore(dgz_struct)
%     dgz = dgz_explore(...)
% dgz_explore, with no inputs, opens up a dialog to select a DGZ file.
% dgz_explore(filename) opens the given file.
% dgz_explore(directory) opens a dialog to select a DGZ file starting in
%     the specified directory.
% dgz_explore(dgz_struct) explores the DGZ structure from an earlier call
%   to DG_READ.
% If an output argument is requested, dgz_explore returns the DGZ
% structure.
%
% INTERACTIVELY EXPLORING A DGZ FILE:
%     When debugging ESS code, it can be useful to be able to quickly and
% easily see what events are being written to file. DGZ_EXPLORE allows you
% to step through the DGZ file observation by observation, seeing which
% events were written, when they were written, and with what parameters.
% DGZ_EXPLORE begins by displaying the first observation, then presenting
% you with a prompt "DGZ>>". At this prompt, you can step through the
% observations by typing "n" or "next", or "p" or "prev" (for previous).
% You can also go to a specific observation by typing its number. When
% finished, type "q" or "quit".

persistent filepath

if nargin < 1 || isempty(file)
    [filename, filepath] = uigetfile('*.dgz','Choose DGZ File',filepath);
    file = fullfile(filepath,filename);
elseif ischar(file) && isdir(file)
    [filename, filepath] = uigetfile('*.dgz','Choose DGZ File', file);
    file = fullfile(filepath,filename);
end

if ischar(file)
    try
        dgz = dg_read(file);
    catch % Maybe this is the stupid dg_read bug, try again:
        dgz = dg_read_local(file);
    end
elseif isstruct(file)
    dgz = file;
else
    error('Can only explore files or structs.');
end

evtNames = cellstr(dgz.e_names);

nObs = numel(dgz.obs_times);

helptext='Valid commands: (n)ext, (p)rev, (q)uit, or observation number.';

printObs(1);
curObs = 1;
done = false;
while ~done
    user_entry = input('DGZ>> ', 's');
    switch lower(user_entry)
        case {'n', 'next'}
            curObs = curObs + 1;
            if curObs > nObs
                curObs = 1;
            end
            printObs(curObs);
        case {'p', 'prev', 'previous'}
            curObs = curObs - 1;
            if curObs < 1
                curObs = nObs;
            end
            printObs(curObs);
        case {'q', 'quit', 'exit'}
            done = true;
        case {'?', 'help'}
            fprintf('%s\n',helptext);
        otherwise
            % Was it a number?
            [n, tf] = str2num(user_entry); %#ok<ST2NM>
            if tf % Yes, it was.
                curObs = n;
                printObs(curObs);
            else % No, unrecognized command.
                fprintf('%s\n',helptext);
            end
    end
end

if nargout > 0
    varargout{1} = dgz;
end

    function printEvent(obs,evt)
        evtTime = dgz.e_times{obs}(evt);
        evtType = dgz.e_types{obs}(evt);
        evtName = evtNames{evtType+1};
        evtSubtype = dgz.e_subtypes{obs}(evt);
        evtParams = dgz.e_params{obs}{evt};
        if ischar(evtParams)
            parStr = evtParams;
        else
            parStr = sprintf('%g ',evtParams);
            parStr = parStr(1:end-1);
        end
        typeStr = sprintf('%s (%d)', evtName, evtSubtype);
        if mod(evt,2)
            typeStr(end+1)    = ' ';
            typeStr(end+1:26) = '.';
        else
            typeStr(end+1:26) = ' ';
        end
        fprintf('%3d)%6d: %s: %s\n', evt, evtTime, typeStr, parStr);
    end
    
    function printObs(obs)
        nEvts = numel(dgz.e_times{obs});
        fprintf('Exploring %s, observation %d/%d at time %g with %d events:\n', ...
            dgz.filename, obs, nObs, dgz.obs_times(obs), nEvts);
        fprintf('    Time   Event Name (Subtype)        Parameters\n');
        for k=1:nEvts
            printEvent(obs,k);
        end
    end
end