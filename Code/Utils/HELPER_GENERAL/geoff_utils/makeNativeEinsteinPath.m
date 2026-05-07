function nativepath = makeNativeEinsteinPath(path)
% MAKENATIVEEINSTEINPATH  Locally native usage of paths to Einstein
% Usage:
%   nativepath = makeNativeEinsteinPath(path)
% Converts a path to a file on Einstein, expressed for either Linux or
% Windows, to the appropriate path for the current operating system.
% Examples:
%   Input:       'Z:\some\file\path'
%   Output:
%     Windows:   'Z:\some\file\path'
%     Linux/Mac: '/einstein0/some/file/path'
%
%   Input:       'L:\some\file\path'
%   Output:
%     Windows:   'L:\some\file\path'
%     Linux/Mac: '/einstein0/USRlab/some/file/path'
%
%   Input:       '/einstein0/USRlab/some/file/path'
%   Output:
%     Windows:   'Z:\USRlab\some\file\path'
%     Linux/Mac: '/einstein0/USRlab/some/file/path'
%
%   Input:       'some/relative/path' or 'some\relative\path'
%   Output:
%     Windows:   'some\relative\path'
%     Linux/Mac: 'some/relative/path'
%
% Written by GKA, Nov 2007

winZ = 'Z:\';
winL = 'L:\';
unixZ = '/einstein0/';
unixL = [unixZ 'USRlab/'];

if strncmp(winZ, path, numel(winZ))
    pathType = 'winZ';
elseif strncmp(winL, path, numel(winL))
    pathType = 'winL';
elseif strncmp(unixZ, path, numel(unixZ))
    pathType = 'unixZ';
else
    if ispc
        convertSep = '\';
    else
        convertSep = '/';
    end
    nativepath = path;
    nativepath(nativepath==convertSep)=filesep;
    return;
end

switch pathType
    case 'winZ'
        if ispc
            nativepath = path;
        else
            nativepath = path(numel(winZ)+1:end);
            nativepath(nativepath=='\')='/';
            nativepath = [unixZ nativepath];
        end
    case 'winL'
        if ispc
            nativepath = path;
        else
            nativepath = path(numel(winL)+1:end);
            nativepath(nativepath=='\')='/';
            nativepath = [unixL nativepath];
        end
    case 'unixZ'
        if isunix
            nativepath = path;
        else
            nativepath = path(numel(unixZ)+1:end);
            nativepath(nativepath=='/')='\';
            nativepath = [winZ nativepath];
        end
end