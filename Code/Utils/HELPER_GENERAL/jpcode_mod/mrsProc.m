function mrsDat = mrsProc(mrsDat, optin);
% %  mrsDat = mrsProc(mrsDat, optin)
% %
% %  processing module for MRS data
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
FCTNAME = 'mrsProc';

FIGUREID = 5;
          
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.ECC     = 0;
dopt.CONVDTA = 1;
dopt.AVG     = 1;
dopt.PLOT    = 1;
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 1;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

f_ecc     = dopt.ECC;
f_convd2a = dopt.CONVDTA;
f_avg     = dopt.AVG;
f_plot    = dopt.PLOT;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

if f_ecc(1) > 0
    eccFilenum = f_ecc(1);
    if length(f_ecc) > 1
        mrsDat = mrsEcc(mrsDat, eccFilenum, opt('THRES',f_ecc(2),'VERBOSE',f_verbose));
    else
        mrsDat = mrsEcc(mrsDat, eccFilenum, opt('VERBOSE',f_verbose));
    end
end

if f_convd2a > 0
    if f_convd2a == 1
        shiftManual = 0;
    else
        shiftManual = f_convd2a;
    end
    mrsDat = mrsConvdta(mrsDat,opt('SHIFT',shiftManual,'VERBOSE',f_verbose));
end

if f_avg(1) > 1
    groupNum = f_avg(1);
    s_tdat = size(mrsDat.tdat);
    dimsLeft = prod(s_tdat)/s_tdat(1)/groupNum;
    if f_verbose
        fprintf('--- %s: AvgGroup = %d, DimsLeft = %d\n', FCTNAME, groupNum, dimsLeft)
    end

    tdatAvg = reshape(mrsDat.tdat, s_tdat(1), groupNum,  dimsLeft);
    tdatAvg = avg(tdatAvg, 2);
    tdatAvg = reshape(tdatAvg, s_tdat(1), dimsLeft );
    datForPlotOrig = mrsDat.tdat(:,1);
    datForPlotAvg  = tdatAvg(:,1);
    mrsDat.tdat = tdatAvg;          %%mrsDat.tdat(:,1:dimsLeft) = tdatAvg;

    mrsDat.flagAvgDone = mrsDat.flagAvgDone + groupNum;

    if f_plot
		figure(FIGUREID)
		subplot(3,1,1)
		plot(abs(fftshift(fft(datForPlotOrig))))
		subplot(3,1,2)
		plot(abs(fftshift(fft(datForPlotAvg))))
		subplot(3,1,3)
		plot(abs(fftshift(fft(datForPlotAvg))))
		for inr = 1:dimsLeft
            plot(abs(fftshift(fft(tdatAvg(:,inr)))))
            ylim([0 max(abs(fftshift(fft(datForPlotAvg))))])
            %key
		end
    end
end

%------------------------------------------------------------
