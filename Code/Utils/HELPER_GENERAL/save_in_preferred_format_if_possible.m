function [format_used] = save_in_preferred_format_if_possible(varargin)
% [FORMAT_USED] = save_in_preferred_format_if_possible(FORMAT_LIST, ...)
%
% where FORMAT_LIST = {format1,format2,...}
%
% Calls save in the caller's workspace with the given (variable number of)
% character arguments (...),
%
% First attempts to save in format1; if that gives a warning message, then
% attempts to save in format2; and so on until has tried all formats in
% FORMAT_LIST.
%
% Currently allowed formats are:
%  'v4', 'v6', 'v7', 'v7.3'
% (note that they are specified without a dash, e.g. 'v7' not '-v7')
%
% If FORMAT_LIST is not specified (i.e. if the first argument is not a cell
% array), it defaults to {'v7','v7.3'}, i.e. first trying to save in the
% old fast & space-efficient format (v7), and then if it fails, trying to
% save in the new space-inefficient format (v7.3).
%
% Returns FORMAT_USED, the format used to save the matfile (or an empty 
% string, if all formats in the list failed to save without giving a 
% warning).

if iscell(varargin{1})
    format_list = varargin{1};
    varargin = varargin(2:end);
else
    format_list = {'v7','v7.3'};
end
assert(all(cellfun(@ischar,format_list)),'format_list must be a cell array of strings');
assert(all(cellfun(@ischar,varargin)),'arguments to save must all be of class char');

format_list_string = strjoin(format_list,' , ');
%disp(['save_in_preferred_format_if_possible: using format list = ' format_list_string]);

% save old last warning state (so can restore it afterward)
[old_msg,old_msgid] = lastwarn();
lastwarn('','');

% join all string arguments together so that the strings, e.g.:
%   savefile.mat
%   var1
%   var2
% are converted into a single string that can be used for
% the call to "evalin", e.g.:
%   'savefile.mat' , 'var1' , 'var2'
save_arg_string = ['''' strjoin(varargin,''' , ''') ''''];

format_used = '';

needed_to_try_other_formats = false;

for fi = 1:numel(format_list)
    save_format_string = ['''-' format_list{fi} ''''];
    save_command_string = ['save(' save_arg_string ' , ' save_format_string ')'];
    %disp(['save_in_preferred_format: attempting to save in format ' format_list{fi} ' using command:']);
    %disp([' ' save_command_string]);
    timer = tic;
    evalin('caller',save_command_string);
    time_taken = toc(timer);
    time_string = ['took ' num2str(time_taken) ' sec'];

    [new_msg,new_msgid] = lastwarn();
    if isempty(new_msgid)
        if needed_to_try_other_formats
            disp(['save_in_preferred_format_if_possible: ...' time_string ', no warning detected, save is presumed successful!']);
        end
        format_used = format_list{fi};
        break;
    else
        disp(['save_in_preferred_format_if_possible: ...' time_string ', warning detected, will try further formats']);
        needed_to_try_other_formats = true;
        lastwarn('','');
    end
end

if isempty(new_msgid)
    % if we didn't receive any warnings, restore old last warning state
    lastwarn(old_msg,old_msgid);
else
    disp('warning: save_in_preferred_format_if_possible: no file formats allowed a successful save without warnings!');
end

