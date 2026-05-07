function stdpath
%%
%%   set global stdpath variables used for data localization
%%   called from startup.m
%%
%%   Sep 2001 - JP

global STDPATH
global HOME

		% !! use backslash at the end of PATH definition !!
        
%HOME = 'D:/';
%[s,HOME] = unix('echo $HOME');
%HOME = strcat(HOME(1:length(HOME)-1), '/')
%STDPATH.pv = 	[HOME 'data/'];

STDPATH.pv = 	['M:/mridata/'];
STDPATH.idl =   ['E:/data_idl/'];
STDPATH.matlab= ['D:/matlab/data/'];
STDPATH.lcm =   ['M:/mpidata/lcm/'];


% remaining VARS NOT yet updated !
%STDPATH.vnmr = 	[HOME 'data/vnmr/'];
%STDPATH.vnp = 	[HOME 'data/vnmr/vnp/'];
%STDPATH.matlab= [HOME 'data/matlab/'];
%STDPATH.idl =   [HOME 'data/idl/'];
%STDPATH.hsvd = 	[HOME 'data/hsvd/'];
%STDPATH.lcm = 	[HOME 'data/lcm/'];
%STDPATH.ascii =	[HOME 'data/ascii/'];
%STDPATH.shape =	[HOME 'vnmrsys/shapelib/'];
%STDPATH.ps = 	[HOME 'data/ps/'];
%STDPATH.hgl = 	'';
%STDPATH.scans =	[HOME 'data/scans/'];

