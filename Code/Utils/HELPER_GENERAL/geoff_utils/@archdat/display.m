function display(ob)
% ARCHOBJ/DISPLAY  Display method for archdat objects

classname = classConst('classname');

% Get the name of the variable ob in the calling workspace:
inName = inputname(1);
if ~isempty(inName)
    if isequal(get(0,'FormatSpacing'),'compact')
        disp([inName ' =']);
    else
        disp(' ');
        disp([inName ' =']);
        disp(' ');
    end
end

if numel(ob)~=1
    sizeStr = size2str(size(ob));
    msg = sprintf('%s %s array', sizeStr, classname);
    disp(msg);
    disp(' ');
elseif isnull(ob)
    disp(sprintf('%s null reference', classname));
    disp(' ');
else
    if ob.inmem
        memState = '"ready"';
    else
        memState = '"not ready"';
    end
    disp(sprintf('%s archive reference (data currently %s):', ...
        classname, memState));
    disp(sprintf('        File: %s', ob.file));
    
    varStr = ob.var;
    disp(sprintf('    Variable: %s', varStr));
    try
        dsize  = datasize(ob);
        dclass = dataclass(ob);
        sizeStr = size2str(dsize);
        disp(sprintf('        Data: %s %s', sizeStr, dclass));
    catch
        disp('        Data: (Unable to read archive!)');
    end
    disp(' ');
end

function sizeStr = size2str(sz)
% Generate a string representing a size vector. Creates strings of the
% format "MxN" for scalar, vector, or matrix sizes, "MxNxP" for
% 3-dimensional arrays, and "N-D" for N>3 dimensional arrays.
dim = numel(sz);
if dim<4
    sizeStr = num2str(sz,'%dx');
    sizeStr = sizeStr(1:end-1);
else
    sizeStr = num2str(dim,'%d-D');
end