function varargout = WorkInStruct(S,varargin)
% WORKINSTRUCT Treat a struct as if it were a workspace
% Usage:
%    WorkInStruct S
%        statements
%    WorkInStruct -end
%
%    WorkInStruct S -option arguments
%        statements
%    WorkInStruct -end -option arguments
%
%    WorkInStruct(S,...)
%        statements
%    T = WorkInStruct('-end',...)
%
% WORKINSTRUCT effectively steps into the "workspace" of a (scalar) struct.
% In this (pseudo) workspace, the "variables" are the fields of the struct.
% This may be particularly useful for dealing with heirarchical structs, in
% which accessing a single data member may require rather cumbersome
% expressions like "structVar.structField1.structField2.field". If many
% values need to be manipulated at the same level of the hierarchy, it may
% be simpler to conceptualize working within the "workspace" of that level.
% In such instances, WORKINSTRUCT may simplify the code considerably by
% removing the need to constantly rereference the whole hierarchy.
%   For example, suppose S1.S2 is a struct with fields 'x', 'y', and 'z'.
% Issuing the command "WorkInStruct S1.S1" will clear the workspace and
% replace it with variables 'x', 'y', and 'z', with values taken from the
% fields of S1.S2. Effectively, you are now working "within" the structure
% S1.S2. A second call of "WorkInStruct -end" will "rewrap" the S1.S2
% workspace back into the struct S1.S2 and restore the original workspace.
% So the following code:
%    WorkInStruct S1.S2
%        q = x + y;
%    WorkInStruct -end
% effectively is the same as writing:
%    S1.S2.q = S1.S2.x + S1.S2.y;
% This is, of course, most useful if there are many statements like this,
% rather than just one. 
%
% Note that the "rewrapping" that occurs at "WorkInStruct -end" will occur
% only if the struct expression supplied at the initiating statement of the
% block is also a valid target of assignment. For example, "S1.S2" will
% successfully rewrap, but the output of a function will not. See 
% "POTENTIAL SOURCES OF BUGS" below.
%
% THE FUNCTION-WITH-RETURN SYNTAX:
%    In the version "T = WorkInStruct('-end')", WORKINSTRUCT rewraps the
% struct not into the original struct S, but into the new struct T. If the
% "-copy" flag is specifed (see OPTIONS, below), this results in both the
% original struct S and a new, modified struct T being present in the
% workspace. Otherwise, the original struct S is lost. The exception is
% that if S is not a variable but instead an expression like "S1.S2", the
% "-copy" behavior is forced when "WorkInStruct('-end')" has an output.
%
% OPTIONS:
%    The initiating and ending calls to WORKINSTRUCT take different option
% flags. Valid options for the initiating call are:
%    -with    : Enters the struct's workspace "with" all variables
%               specified after the flag. These variables remain
%               accessible, but are still part of the calling workspace.
%               Note that an expression like "WorkInStruct S -with S" is
%               permissible; S is unchanged while working in the struct,
%               but might be overwritten by the ending call. Also note that
%               in the event of a naming conflict between a -with variable
%               and a field in S, the fields of S take precedence.
%    -withall : Enters the struct's workspace "with" all variables in the
%               calling workspace, except for the variable specifying the
%               struct. To include the struct as well, use the expression
%               "WorkInStruct S -with S -withall".
%    -copy    : Works in a copy of the struct rather than in the original
%               struct. At "WorkInStruct -end", the struct is restored to
%               its original state.
% Valid options for the ending call are:
%    -clear   : Clears all variables specified after the flag. This is
%               a shorthand for calling CLEAR on those variables just
%               before the ending call.
%    -export  : "Exports" the variables specified after the flag to the
%               calling workspace. The exported variables are not stored as
%               fields in the struct.
% 
% POTENTIAL SOURCES OF BUGS:
% Each "WorkInStruct" call MUST be followed by a call to
% "WorkInStruct -end". The "WorkInStruct ... WorkInStruct -end"
% construction should be thought of as delineating a block of code in the
% same way that "if ... end" or "for ... end" does. Unlike these cases, the
% MATLAB parser will not generate an error if "WorkInStruct -end" is not
% called, and failure to call "WorkInStruct -end" appropriately may lead to
% unexpected behavior. Some examples to illustrate the point:
%   Wrong:                        Right:
%     function y = F(x)             function y = F(x)
%     WorkInStruct x                WorkInStruct x
%     statements                        statements
%                                   WorkInStruct -end
%   Wrong:                        Right:
%     WorkInStruct s                WorkInStruct s
%     if expr                           if expr
%         statements                        statements
%         WorkInStruct -end             end
%     end                           WorkInStruct -end
%   Wrong:                        Right, but may be inefficient:
%     for expr                      for expr
%         WorkInStruct s                WorkInStruct s
%         statements                        statements
%     end                               WorkInStruct -end
%     WorkInStruct -end             end
% In short, treat "WorkInStruct -end" the same way you would treat the
% "end" of an if, for, while, switch, or any other statement that begins a
% discrete block of code. Following good indentation practice will help
% keep this straight.
%    However, it is important to note that if a WORKINSTRUCT block is
% defined within a code block like while...end that can be escaped with
% BREAK, or within a function that uses RETURN, "WorkInStruct -end" *must*
% be called prior to BREAK or RETURN. An example:
%    function y = f(x)
%    WorkInStruct(x)
%        if condition
%            WorkInStruct -end   % <-- Add this line
%            return
%        end
%    WorkInStruct -end
% This should be done if and only if the BREAK or RETURN statement would
% cause the program flow to move beyond the "WorkInStruct -end" statement.
% Similarly, if a WORKINSTRUCT block is defined within a try...catch block,
% then the "catch" portion of the block should call "WorkInStruct -end" as
% necessary.
%
% As a side effect, "WorkInStruct -end" produces a struct in which the
% field names are alphabetized. This may be undesirable in some cases. If
% the "-copy" flag is specified, then only the copy is alphabetized; the
% original is left pristine.
%
% If WORKINSTRUCT is called on a struct indexed from an array of structs
% (e.g., S(2)), WORKINSTRUCT may not be able to rewrap the struct into the
% original array element. In particular, if the statements executed in the
% WORKINSTRUCT block add or remove fields to the struct, the "rewrapping"
% will generate a "Subscripted assignment between dissimilar structures"
% error. WORKINSTRUCT will catch the error and handle it, but the original
% struct will not be reassigned. In effect, this is a forced "-copy".
% Additionally, even if no fields are added or removed, the alphabetization
% side effect mentioned above means that the result of "WorkInStruct -end"
% is effectively dissimilar to the original struct unless the original's
% fields were already alphabetized; again, this is like a forced "-copy".
% 
% PERFORMANCE:
% WORKINSTRUCT involves a certain amount of overhead each time it is
% called, and so may not be appropriate for applications in which speed is
% at a premium.

% Basic flow pattern:
% On WorkInStruct S:
%    1. If S is a string and evaluates to a struct, copy it to the
%       structName field of the database. If it is a struct and has an
%       inputname, copy that to the structName field. Otherwise leave the
%       structName field empty.
%    2. Check for -with variables. If any -with variables are also field
%       names in S, remove them from the keep list and give a warning.
%    3. If -copy is specified, copy S to the struct field of the database.
%    4. If S is a variable name or has an inputname, clear it from the
%       (calling) workspace.
%    4. Copy all variables not in the keep list to the ws field of the
%       database and then clear them from the workspace.
%    5. Copy all the fields of S to variables in the workspace.
% On WorkInStruct -end:
%    1. Check for -export and -clear variables
%    2. Clear the variables in the clear list
%    3. If the structName field exists, copy the values of all variables
%       not in the export list to fields in a struct. Then clear them from
%       the workspace.
%    4. If possible, assign the value of that struct to the structName
%       expression (should be either a variable name or a valid target of
%       subscripted assignment).
%    5. Assign variables in the workspace according to the persistent
%       database.

% Define the persistent WorkInStruct database:
persistent WISdb
if isempty(WISdb)
    WISdb = struct('structName',{}, 'ws',{}, 'struct',{}, 'keep',{});
end

endFlag = ischar(S) && strcmpi(S,'-end');
if ~endFlag % Start of block
    % Usages:
    % WorkInStruct S (options)
    % WorkInStruct(S,(options))
    error(nargoutchk(0,0,nargout)); % No output arguments for this syntax
    if ischar(S)
        % Command form:
        % WorkInStruct S (options)
        wsStructName = S;
        % We will use evalin without checking isvarname because we would
        % like to allow the user to be able to use any MATLAB expression
        % that evaluates to a struct. Hopefully that expression will also
        % be a valid target for assignment, such as "S1.S2".
        wsStruct = evalin('caller',S);
    else
        % Function form:
        % WorkInStruct(S,(options))
        wsStruct = S;
        wsStructName = inputname(1);
    end
    if ~isstruct(wsStruct)
        error('WorkInStruct:nonStruct','Input must be a struct!');
    end
    if ~isscalar(wsStruct)
        error('WorkInStruct:nonScalar','Input must be scalar!');
    end
    
    % Get the names of the variables in the calling workspace:
    callerVars = evalin('caller','who');
    if isvarname(wsStructName) && ~ismember(wsStructName,callerVars)
        % wsStructName is a valid variable name, but is not a variable in
        % the calling workspace, which probably means it is a function. So
        % don't actually store its name.
        wsStructName = '';
    end
    % Parse the options:
    [keepVars,copyFlag] = parse_start(varargin,callerVars,wsStructName);
    
    % Determine what variables are going to be defined by the struct:
    wsStructVars = fieldnames(wsStruct);
    % Are any of the -with variables also defined by the struct?
    [tf,loc] = ismember(wsStructVars,keepVars);
    if any(tf) % If so, remove them. The struct's fields take precedence.
        warning('WorkInStruct:keepOverriden', ...
            ['One or more of the variables specified with -with are ' ...
            'defined by the workspace structure, and so cannot be kept.']);
        keepVars = keepVars(setdiff(1:end,loc));
    end
    
    % Determine which variables need to be stored in the WISdb (all that
    % are not specified with -with except the workspace struct)
    backupVars = setdiff(callerVars, [keepVars(:);{wsStructName}]);
    % Store those variables into a struct:
    backupStruct = struct;
    for k=1:numel(backupVars)
        backupStruct.(backupVars{k}) = evalin('caller',backupVars{k});
    end
    
    % Create a new record in the WISdb:
    newEntry.structName = wsStructName;
    newEntry.ws = backupStruct;
    if copyFlag
        newEntry.struct = wsStruct;
    else
        newEntry.struct = struct([]);
    end
    newEntry.keep = keepVars;
    if isempty(WISdb)
        mlock;
    end
    WISdb(end+1) = newEntry;
    
    % Clear the workspace struct and all variables that have been backed
    % up:
    if isvarname(wsStructName)
        evalin('caller',['clear ' wsStructName]);
    end
    if ~isempty(backupVars)
        clearCmd = ['clear', sprintf(' %s',backupVars{:})];
        evalin('caller',clearCmd);
    end
    
    % Load the new workspace variables:
    for k=1:numel(wsStructVars)
        assignin('caller',wsStructVars{k},wsStruct.(wsStructVars{k}));
    end
    
    
else % End of block
    % Usages:
    % WorkInStruct -end (options)
    % WorkInStruct('-end',(options))
    % S = WorkInStruct('-end',(options))
    error(nargoutchk(0,1,nargout)); % Maximum of 1 output
    
    if isempty(WISdb)
        error('WorkInStruct:endWithoutStart', ...
            'Cannot use -end option without first starting WorkInStruct.');
    end
    
    % Get the names of the variables in the "struct workspace":
    callerVars = evalin('caller','who');
    % Determine which, if any, variables are to be exported or discarded:
    [exportVars,discardVars] = parse_end(varargin,callerVars);
    if ~isempty(discardVars)
        % Discard the discards:
        clearCmd = ['clear' sprintf(' %s',discardVars{:})];
        evalin('caller',clearCmd);
        callerVars = setdiff(callerVars,discardVars);
    end
    % Load the entry on the top of the WISdb stack:
    curEntry = WISdb(end);
    copyFlag = ~isempty(curEntry.struct);
    exportVars = [exportVars(:);curEntry.keep(:)];
    % Determine which variables are not being exported:
    noExportVars = setdiff(callerVars,exportVars);
    if (~isempty(curEntry.structName)&&~copyFlag) || nargout > 0
        % We need to capture the "struct workspace" to rewrap back into a
        % struct (export vars are not rewrapped):
        wsStruct = struct;
        for k=1:numel(noExportVars);
            wsStruct.(noExportVars{k}) = evalin('caller',noExportVars{k});
        end
    end
    if ~isempty(noExportVars)
        % Remove everything that's not being exported:
        clearCmd = ['clear' sprintf(' %s',noExportVars{:})];
        evalin('caller',clearCmd);
    end
    % Copy the original workspace back:
    origVars = fieldnames(curEntry.ws);
    restoreVars = setdiff(origVars, exportVars);
    for k=1:numel(restoreVars)
        assignin('caller',restoreVars{k},curEntry.ws.(restoreVars{k}));
    end
    if nargout > 0
        % If there is an output argument, then the workspace struct is
        % returned:
        varargout{1} = wsStruct;
    end
    if copyFlag
        wsStruct = curEntry.struct;
    end
    if (~isempty(curEntry.structName) && nargout == 0) || copyFlag
        % If there is no output argument but the original structName can be
        % resolved, attempt to store the rewrapped workspace struct in its
        % original value:
        if isvarname(curEntry.structName)
            % This is easy if it was a variable name:
            assignin('caller',curEntry.structName,wsStruct);
        else
            % It was probably the result of a subscripting, e.g.,
            % S.field. So we can still try to reassign it.
            % Create a temporary variable in the parent workspace to hold
            % it:
            tmpname = evalin('caller','genvarname(''tmp'',who)');
            assignin('caller',tmpname,wsStruct);
            % Attempt to use the assignment operator on the expression
            % (perfectly valid if it was a subscripting expression):
            assgnCmd = sprintf('%s = %s;', curEntry.structName, tmpname);
            try
                evalin('caller',assgnCmd);
            catch
                % It didn't work. This can happen if the struct was the
                % result of a subscripting expression like S(3), and the
                % names of the fields changed.
                warning('WorkInStruct:noRewrap', ...
                    ['Unable to rewrap the struct workspace into the ' ...
                    'original struct.']);
            end
            % Clean up the temporary variable:
            evalin('caller',['clear ' tmpname]);
        end
    else
        % The user didn't try to capture the rewrapped workspace struct,
        % and the database didn't have a name for the struct. Warn the user
        % and move on.
        warning('WorkInStruct:noRewrap', ...
            ['Unable to rewrap the struct workspace into the ' ...
            'original struct.']);
    end
    
    % Pop the last element from the WISdb stack.
    WISdb = WISdb(1:end-1);
    if isempty(WISdb)
        munlock;
    end
end

function [keepVars, copyFlag] = parse_start(args,callerVars,structName)
if ~iscellstr(args)
    error('WorkInStruct:badArguments', ...
        'Additional options must be flags or variable names.');
end
keepVars = {};
copyFlag = false;
mode = '';
for k=1:numel(args)
    if strcmpi(args{k},'-withall')
        if ismember(structName,keepVars)
            keepVars = callerVars;
        else
            ind = ~strcmp(structName,callerVars);
            keepVars = callerVars(ind);
        end
        mode = '';
    elseif strcmpi(args{k},'-with')
        mode = 'with';
    elseif strcmpi(args{k},'-copy')
        copyFlag = true;
        mode = '';
    else
        switch mode
            case 'with'
                if ismember(args{k},callerVars)
                    if ~ismember(args{k},keepVars)
                        keepVars{end+1} = args{k};
                    end
                else
                    warning('WorkInStruct:keepNotDefined', ...
                        'Undefined variable %s, ignoring.', args{k});
                end
            otherwise
                error('WorkInStruct:unrecognizedSyntax', ...
                    'Unrecognized calling syntax.');
        end
    end
end

function [exportVars,discardVars] = parse_end(args,callerVars)
if ~iscellstr(args)
    error('WorkInStruct:badArguments', ...
        'Additional options must be flags or variable names.');
end
exportVars = {};
discardVars = {};
mode = '';
for k=1:numel(args)
    if strcmpi(args{k},'-export')
        mode = 'export';
    elseif strcmpi(args{k},'-clear')
        mode = 'discard';
    else
        switch mode
            case 'export'
                if ismember(args{k},callerVars)
                    if ~ismember(args{k},exportVars)
                        exportVars{end+1} = args{k};
                    end
                else
                    warning('WorkInStruct:exportNotDefined', ...
                        'Undefined variable %s, ignoring.', args{k});
                end
            case 'discard'
                if ismember(args{k},callerVars)
                    if ~ismember(args{k},discardVars)
                        discardVars{end+1} = args{k};
                    end
                else
                    warning('WorkInStruct:discardNotDefined', ...
                        'Undefined variable %s, ignoring.', args{k});
                end
            otherwise
                error('WorkInStruct:unrecognizedSyntax', ...
                    'Unrecognized calling syntax.');
        end
    end
end