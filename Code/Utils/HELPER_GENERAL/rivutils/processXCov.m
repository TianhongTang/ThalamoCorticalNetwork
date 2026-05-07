function XCV = processXCov(data,sampt,method,win,lags,shuffle,file,A,E,maxpts,triangular,C)

global FILES

procfilename = sprintf('%s/%s/%s', FILES.PROC_XCVPath,FILES.FILERoot,file);

% check whether the file exists or not
XCV = checkFile(procfilename);
if ~isempty(XCV)
  fprintf(' processXcov: Successfully read XCV file %s\n', procfilename);
  return
end;

% check argins
len  = size(data,1);   
npts = floor(len/win);
if nargin < 10,  maxpts = npts;   end
if nargin < 11,  triangular = 1;  end
if nargin < 12,  C = [];          end

% cut data into windows
nwin = min(maxpts,npts);
fprintf(' processXCov: cutting into %d windows .',nwin);
cutdat  = zeros(win,nwin,size(data,2));
for i=1:nwin
  cutdat(:,i,:) = data([(i-1)*win+1:i*win],:);
  %if ~mod(i,100),  fprintf('.');  end
end
data = [];  % no more need
fprintf(' done\n');

% compute xcov/xcorr
fprintf(' processXCov: computing xcov ');
XCV.xcdata     = [];
XCV.shuffle    = [];
XCV.rawdata    = data;
[XCV.xcdata, XCV.shuffle] = computeXCov(cutdat,method,win,lags,shuffle,triangular);
XCV.method     = method;
XCV.win        = win;
XCV.lags       = lags;
XCV.sampt      = sampt;
XCV.electrodes = E;
XCV.areas      = A;
XCV.channels   = C;
fprintf(' done\n');
cutdat = [];

% save data
if ~isempty(XCV.xcdata)
  fprintf(' processXCov: writing to %s',file);
  procdir = sprintf('%s/%s',FILES.PROC_XCVPath,FILES.FILERoot);
  tmp = mkdir(FILES.PROC_XCVPath,FILES.FILERoot);
  save(procfilename,'XCV');
  fprintf(' done\n');
end
