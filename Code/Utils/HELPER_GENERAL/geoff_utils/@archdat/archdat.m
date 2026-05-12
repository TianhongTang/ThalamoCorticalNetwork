function ob = archdat(varargin)
% ARCHDAT  Archived data reference class
% archdat objects are references to large datasets stored on disk in MAT
% files. In general, accessing archdat data is slow, because it requires
% disk access. However, when handling many large datasets at once, it may
% be necessary to control exactly what data MATLAB is storing in memory.
% archdat objects make it easier to juggle these large datasets.
%
% Usage:
%   ad = archdat
%   ad = archdat(n)
%   ad = archdat(sz)
%   ad = archdat(filename)
%   ad = archdat(filename, varname)
%   ad = archdat(filename, varname, data)
%   ad = archdat(copyobj)
%
% archdat, with no inputs, returns a "null" archdat referring to no data.
%   The null archdat cannot have data assigned with the "set" method; it
%   can only be replaced with a non-null archdat.
% archdat(n), where n is a positive integer scalar, creates a column vector
%   of null archdat objects with n elements
% archdat(sz), where sz is a size vector, creates an array of null archdat
%   of size sz.
% archdat(filename) creates a column vector of archdat objects, one
%   referring to each variable in the file.
% archdat(filename, varname) has three usages. If filename and varname are
%   both strings, it creates an archdat object referring to the variable
%   varname in the file filename. If varname is a cell array of strings, it
%   creates an array of archdat objects referring to the variables
%   specified. If filename and varname are both cell arrays of strings with
%   the same number of elements, it creates an array of archdat objects
%   referring to the variables in the filename-varname pairs. Any other
%   usage is an error.
% archdat(filename, varname, data) stores the value of data in the variable
%   varname in the MAT-file filename and creates an archdat referring to
%   it. If filename does not exist, it is created; if it already exists and
%   does not contain the variable varname, varname is appended to the file.
%   If varname already exists in the file, it is overwritten.
% archdat(copyobj), where copyobj is an archdat, is the copy constructor
%   usage, and returns an archdat identical to copyobj.
% 
% The big idea:
%    An archdat points to data stored in a MAT-file. The archdat object may
% be used to retrieve this data in a quick and easy way, as if it were
% actually available in the workspace. This is achieved with the GET
% method, "data = get(ad)", where ad is a (scalar) archdat object.
% 
% archdat methods:
%   file = filename(ad)   get the filename of the archive
%   var = varname(ad)     get the variable name in the archive
%   data = get(ad)        get archdat data (from disk or memory)
%   ad = set(ad,data)     change archdat data value (on disk and in memory)
%   sz = datasize(ad)     get the dimensions of the data
%   class = dataclass(ad) get the class of the data
%   ad = ready(ad)        bring data out of archive and into memory
%   ad = release(ad)      clear data from memory, leaving archive untouched
%   tf = isready(ad)      get the readiness state of the archdat
%   tf = isnull(ad)       determine whether this is a null reference
%   infoStr = info(ad)    get information about the archdat
%
% See also: ARCHDAT/FILENAME, ARCHDAT/VARNAME, ARCHDAT/GET, ARCHDAT/SET,
%           ARCHDAT/DATASIZE, ARCHDAT/DATACLASS, ARCHDAT/READY,
%           ARCHDAT/RELEASE, ARCHDAT/ISREADY, ARCHDAT/ISNULL, ARCHDAT/INFO
%
% Written by GKA, Nov. 2007

% Make constants persistent for negigible speed boost
persistent classname unknownConstructorErr empty_archdat

% Define constants
if isempty(classname)
    classname = classConst('classname');
    
    unknownConstructorErr.message = ...
        sprintf('Unrecognized constructor usage, try "help %s"',classname);
    unknownConstructorErr.identifier = [classname ':unknownConstructor'];
    
    empty_archdat.file = '';
    empty_archdat.var  = '';
    empty_archdat.data = [];
    empty_archdat.inmem = false;
end

switch nargin
    case 0
        % Null constructor usage
        ob = class(empty_archdat,classname);
    case 1
        if isa(varargin{1},classname)
            % Copy constructor usage
            ob = varargin{1};
        elseif isnumeric(varargin{1})
            % Multiple null constructor usage
            sz = varargin{1};
            szsz = size(sz);
            if szsz(1)~=1 || ndims(sz) > 2
                error([classname ':nonRowSize'], ...
                    'Size vector must be a row vector.');
            end
            if isscalar(sz)
                sz = [sz 1];
            end
            ob = repmat(empty_archdat,sz);
            ob = class(ob,classname);
        elseif ischar(varargin{1})
            % Whole file constructor usage
            filename = varargin{1};
            checkArchive(filename);
            varname  = who('-file',filename);
            ob = repmat(empty_archdat,size(varname));
            [ob.file] = deal(filename);
            [ob.var] = deal(varname{:});
            ob = class(ob,classname);
        else
            error(unknownConstructorErr);
        end
        
    case 2
        % Reference to existing data in file
        filename = varargin{1};
        varname  = varargin{2};
        if ischar(filename) && ischar(varname)
            % Scalar constructor usage
            checkArchive(filename,varname);
            ob = empty_archdat;
            ob.file = filename;
            ob.var = varname;
        elseif ischar(filename) && iscellstr(varname)
            % Multiple variables in a file usage
            checkArchive(filename,varname);
            ob = repmat(empty_archdat,size(varname));
            [ob.file] = deal(filename);
            [ob.var] = deal(varname{:});
        elseif iscellstr(filename) && iscellstr(varname) ...
                && numel(filename)==numel(varname)
            % Multiple variables and files usage
            for k=1:numel(filename)
                checkArchive(filename{k},varname{k});
            end
            ob = repmat(empty_archdat,size(filename));
            [ob.file] = deal(filename{:});
            [ob.var]  = deal(varname{:});
        else
            error(unknownConstructorErr);
        end
        ob = class(ob,classname);
    case 3
        % Archive creation usage
        filename = varargin{1};
        varname  = varargin{2};
        data     = varargin{3};
        ob = empty_archdat;
        if ~(isa(filename,'char') && isa(varname,'char'))
            error(unknownConstructorErr);
        end
        opts = {};
        if exist(filename,'file')==2
            opts = [opts {'-append'}];
        end
        
        s.(varname) = data;
        ob.var = varname;
        save(filename,'-struct','s',opts{:});
        ob.file = filename;
        ob.data = data;
        ob.inmem = true;
        ob = class(ob,classname);
    otherwise
        error(unknownConstructorErr);
end
