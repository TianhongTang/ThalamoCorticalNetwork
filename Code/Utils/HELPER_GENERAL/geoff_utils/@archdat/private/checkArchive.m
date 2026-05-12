function checkArchive(filename, varname)

classname = classConst('classname');

if exist(filename,'file')~=2
    error([classname ':fileNotFound'], 'Archive file "%s" not found.', ...
        filename);
end

if nargin > 1
    if ischar(varname)
        if isempty(who(varname,'-file',filename))
            error([classname ':noSuchVariable'], ...
                'Variable "%s" in not found in archive file "%s".', ...
                varname, filename);
        end
    elseif iscell(varname)
        fileVars = who(varname{:},'-file',filename);
        inFile = ismember(varname,fileVars);
        if ~all(inFile)
            missingVars = varname(~inFile);
            error([classname ':noSuchVariable'], ...
                'Variable "%s" in not found in archive file "%s".', ...
                missingVars{1}, filename);
        end
    end
end