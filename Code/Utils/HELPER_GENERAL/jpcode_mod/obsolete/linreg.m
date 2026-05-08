function [yfit,m,a0,bint,r,rint,stats] = linreg(x, y, flagVerbose)
%
% modified version for onedim. REGRESS()
% INPUT X, Y
% OUTPUT yfit, result and plotted fit with data
%
if nargin<3
    flagVerbose = 0; 
end;

X = [ones(length(x),1) x];
[b bint r rint stats] = regress(y, X);
m = b(2);
a0 = b(1);
yfit = a0 + m.*x;

if flagVerbose
    fprintf('slope: %f\nA0:   %f\nR2/Fval/pval: %f/%f/%f\n', m, a0, stats(1:3))
    plot(x,y,'.')
    hold on
    plot(x,yfit)
    hold off
end

return
%-------------------------------------------------------
