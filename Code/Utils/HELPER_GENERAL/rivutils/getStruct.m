function getStruct()
% getStruct(); returns most essential fields for DAT structure
% FUNCTION: getStruct()
%           checks the file date and reads out events accordingly
%           the fields that contain information about the type of 
%           experiment ('free running' vs. 'unambiguous - ambiguous')
%           as well as the valid observation periods:
%           DAT.complete_obs / complete_prevbias_obs / complete_unambiguous_obs
% ARGUMENTS:
%  none
%
% RETURNS:
%	adds fields to global variable DAT
% 	
% 08.10.02 AM

global DGZ

[dir, filename, ext] = fileparts(DGZ.filename);

switch getFileDateNum(filename)
    case getFileDateNum('120802')
          makeDatStruct_prefinalfix;
    case getFileDateNum('130802')
          makeDatStruct_prefinalfix;
    case getFileDateNum('140802')
          makeDatStruct_prefinalfix;
    case getFileDateNum('150802')
          makeDatStruct_prefinalfix;
    case getFileDateNum('160802')
          makeDatStruct_prefinalfix;
    case getFileDateNum('170802')
          makeDatStruct_prefinalfix; 
    case getFileDateNum('180802')
          makeDatStruct_prefinalfix;
    case getFileDateNum('190802')
         makeDatStruct_prefinalfix;
    case getFileDateNum('200802')
          makeDatStruct_prefinalfix;
    case getFileDateNum('210802')
          makeDatStruct_prefinalfix;
    case getFileDateNum('020902')
          makeDatStruct_prefinalfix; 
    case getFileDateNum('030902')
          makeDatStruct_prefinalfix;
    case getFileDateNum('040902')
          makeDatStruct_prefinalfix; 
    case getFileDateNum('090902')
          makeDatStruct_prefinalfix;
    case getFileDateNum('100902')
          makeDatStruct_prefinalfix;
    case getFileDateNum('110902')
          makeDatStruct_prefinalfix;
    case getFileDateNum('120902')
          makeDatStruct_prefinalfix;
    case getFileDateNum('130902')
          makeDatStruct_prefinalfix;
    case getFileDateNum('160902')
          makeDatStruct_prefinalfix; 
    case getFileDateNum('170902')
          makeDatStruct_prefinalfix;
    case getFileDateNum('180902')
          makeDatStruct_prefinalfix; 
    case getFileDateNum('190902')
          makeDatStruct_prefinalfix; 
    case getFileDateNum('210902')
          makeDatStruct_prefinalfix;
    case getFileDateNum('230902')
          makeDatStruct_prefinalfix; 
    case getFileDateNum('270902')
          makeDatStruct_prefinalfix;
    case getFileDateNum('290902')
          makeDatStruct_prefinalfix; 
    case getFileDateNum('300902')
          makeDatStruct_prefinalfix;
    case getFileDateNum('300902')
          makeDatStruct_prefinalfix;
    otherwise
          makeDatStruct;
end  


%------------------------------------------------------
function makeDatStruct_prefinalfix()
global DGZ DAT MUA

% create a structure for further analysis
if ~isempty(MUA)
DAT.nchan                       = size(MUA{1},2);
end
DAT.complete_obs                = findCompleteTrials(DGZ);
DAT.complete_prevbias_obs       = findCompletePrevBiasTrials_pff(DGZ);
DAT.complete_unambiguous_obs    = findCompleteBiasTrials_pff(DGZ);


%------------------------------------------------------
function makeDatStruct()
global DGZ DAT MUA

% create a structure for further analysis
if ~isempty(MUA)
DAT.nchan                       = size(MUA{1},2);
end
DAT.complete_obs                = findCompleteTrials(DGZ);
DAT.complete_prevbias_obs       = findCompletePrevBiasTrials(DGZ);
DAT.complete_unambiguous_obs    = findCompleteBiasTrials(DGZ);
