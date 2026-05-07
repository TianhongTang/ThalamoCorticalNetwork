function bytes = bytesize(thing, units)
% BYTESIZE   Get the amount of memory used by a value
%      bytes = bytesize(val)
%      bytes = bytesize(val, unit)
%    Typically, it's easiest to find the amount of memory used by a MATLAB
%    variable by using WHOS. However, if you're curious about the memory
%    usage of something other than a variable (eg. one field of a
%    structure), WHOS does not readily give that information from the
%    command prompt.
%    The optional "unit" argument returns the memory usage in the specified
%    units, rather than in bytes. Valid options are:
%      'b': bytes
%      'k': kibibytes   (1024 bytes)
%      'm': mebibytes   (1024^2 bytes)
%      'g': gibibytes   (1024^3 bytes)
%      't': tebibytes   (1024^4 bytes)
%
%    Examples:
%      bytesize(a,'m')
%      bytesize(structure.field)/bytesize(structure)
%    
%    GKA July 2007

if nargin<2
    units = 'b';
end

switch units
    case 'b'
        convert = 1;
    case 'k'
        convert = 1024;
    case 'm'
        convert = 1024^2;
    case 'g'
        convert = 1024^3;
    case 't'
        convert = 1024^4;
    otherwise
        warning('GKA:huhWhat', ...
            'Unrecognized unit character, returning bytes.');
        convert = 1;
end

info = whos('thing');
bytes = info.bytes/convert;