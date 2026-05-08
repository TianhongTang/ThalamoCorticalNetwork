function stm = getSTMFileInfo(stmfile)
% PURPOSE : This function retrieves stmfile information.
% USAGE :   stm = getSTMFileInfo2(stmfile)

global STM
  
% init outpus
stm = [];

% pickup a file
if nargin < 1
   stmfile = pickfile('Load STM File',pwd,'*.stm');
   if ~length(stmfile), return; end
   fprintf('STMFILE: %s\n',stmfile);
end
stmdir  = getFileDirectory(stmfile);
stmfile = getFileName(stmfile);
if ~length(stmdir), stmdir = pwd; end

% read text
stmfullpath = sprintf('%s/%s',stmdir,stmfile);
stm.fullpath = stmfullpath;
stmline = loadTclFile(stmfullpath);

% process text
for i = 1:length(stmline)
  token = sscanf(stmline{i},'set %s');
  value = sscanf(stmline{i},'set %*s %f');
  if isstr(value)	
    value = str2num(value);
  end
  if isempty(value)
    value = sscanf(stmline{i},'set %*s %s'); 
  end
  if isstr(token) & ~isempty(value)
    p0 = findstr(token,'(');
    p1 = findstr(token,')');
    if ~isempty(p0) & ~isempty(p1)
      token = token([1:p0-1 p0+1:p1-1]);
    end
    if isstr(value)
      eval(sprintf('stm.%s = ''%s'';',token,value));
    else
      eval(sprintf('stm.%s = %f;',token,value));
    end
  end	
end

if nargout == 0, STM = stm; end
