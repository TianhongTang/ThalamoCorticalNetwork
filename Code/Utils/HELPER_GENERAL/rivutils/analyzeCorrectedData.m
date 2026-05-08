function analyzeCorrectedData(DGZc,MUAc)
% analyzeCorrectedData(DGZc,MUAc); use corrected MUA and DGZ for further analyzis
% FUNCTION: analyzeCorrectedData(DGZc,MUAc)
%           make MUA = MUAc and DGZ = DGZc and go on with analysis
% ARGUMENTS:
%  DGZc, MUAc (from individFileFix.m)
%
% RETURNS:
%	nothing (modifies DGZ and MUA global)
% 	
% 08.10.02 AM

global DGZ DAT MUA 

DGZ = DGZc;
MUA = MUAc;

