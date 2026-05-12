function mrsDat = mrsBC(mrsDat, optin);
% %  mrsDat = mrsBC(mrsDat, optin)
% %
% %  baseline correction of MRS data
% %
% %  default options:
% %     opt(  
% %           'VERBOSE','1')
% %
% %  return: 
% %
% %  Functions called: 
% %  Tested for: 
% %
% %  Nov 2002 -  Josef Pfeuffer
% %
FCTNAME = 'mrsBC';

DCSTARTFT = 0.95;	% startof relevant data (BCFT - correction)
DCENDFT   = 0.99;	% edn of  relevant data (BCFT - correction)
          
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.BCFT    = 1;     %% corrects 1st time domain data point (DC correction in FT domain)
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 1;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

f_bcft    = dopt.BCFT;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

s_tdat = size(mrsDat.tdat);
td = s_tdat(1);
nr = s_tdat(2);
nelemBcft = fix( (DCENDFT-DCSTARTFT)*td );
bcftLeft  = fix( (1-DCENDFT)*td );
bcftRight = fix( DCENDFT*td );

if f_verbose
    fprintf('--- %s: BCFT limits [%.3f,%.3f]\n', FCTNAME, DCSTARTFT, DCENDFT)
end
for inr = 1:nr      %% 2D loop
    
    dat = mrsDat.tdat(:,inr);

    datFt = fftshift(fft(dat,[],1),1); 
%    bcArr = [datFt(bcftLeft:bcftLeft+nelemBcft-1)' datFt(bcftRight-nelemBcft+1:bcftRight)']';
      % left side only
    bcArr = datFt(bcftLeft:bcftLeft+nelemBcft-1);
    bcMeanReal = sum( real(bcArr) )/ length(bcArr);
    bcMeanImag = sum( imag(bcArr) )/ length(bcArr);
    datFt = complex(real(datFt) - bcMeanReal, imag(datFt) - bcMeanImag);
    dat = ifft(ifftshift(datFt,1),[],1);

    mrsDat.tdat(1,inr) = dat(1);  %% correct first point in FID
end
%------------------------------------------------------------
