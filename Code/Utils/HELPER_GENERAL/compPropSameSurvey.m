% compPropSameSurvey         Compare proportions from the same survey
% 
%     [h,pval,z,ci,pd] = compPropSameSurvey(p,n,C,alpha,tail,flag)
%
%     For example, if you surveyed a group of people and asked them if they 
%     prefer Pepsi or Coke, and you wanted to know whether these proportions 
%     differed within this sample, you should use COMPPROPSAMESURVEY instead 
%     of a t-test, standard z-test or chi2 test, which assume the proportion 
%     of Pepsi preferers and Coke preferers were indepedently sampled from 
%     two different surveys.
%
%     INPUTS
%     p       - vector of probabilities to compare
%     n       - sample size
%
%     OPTIONAL
%     C       - Contrast matrix. Each row represents a contrast and should 
%               have the same number of columns as elements in p. The elements 
%               should be +1 or -1 to indicate which elements of p to compare.
%               Multiple +1s or -1s allow elements to be grouped. 
%               Defaults to all pairwise comparisons of the elements in p.
%     alpha   - significance level (100*alpha), default=0.05
%     tail    - 'both'  p1 ~= p2 (default)
%               'right' p1 > p2
%               'left'  p1 < p2
%     flag    - 's' for simultaneous confidence intervals (Bonferroni)
%               'i' for individual confidence intervals (default)
%     verbose - set false to suppress output to stdout (defaults true)
%
%     OUTPUTS
%     h        - boolean indicating whether null is rejected 
%     pval     - p-value (uncorrected)
%     z        - z-value 
%     ci       - confidence intervals
%     pd       - proportion differences for requested contrasts
%
%     EXAMPLE
%     % Reproduce example from Berry & Hurtado (1994)
%     r = [114 71 175 40]; n = sum(r);
%     compPropSameSurvey(r/n,n,[1 -1 0 0 ; 1 1 -1 1],[],[],'i');
%     compPropSameSurvey(r/n,n,[1 -1 0 0 ; 1 1 -1 1],[],[],'s');
%
%     REFERENCE
%     Scott & Seber (1983) Difference of proportions from the same survey.
%       The American Statistician 37(4): 319-320
%     Berry & Hurtado (1994) Comparing non-independent proportions.
%       Observations: The Technical Journal for SAS Software Users 3(4): 21-27. 
%
%     SEE ALSO
%     compProp2

%     $ Copyright (C) 2011 Brian Lau http://www.subcortex.net/ $
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%     REVISION HISTORY:
%     brian 12.21.11 written
%     brian 01.12.11 added printout to stdout, fixed failure to pass flag
%                    w/ multiple contrasts

function [h,pval,z,ci,pd] = compPropSameSurvey(p,n,C,alpha,tail,flag,verbose)

if (nargin < 7) || isempty(verbose)
   verbose = true;
end
if (nargin < 6) || isempty(flag)
   flag = 'i';
end
if (nargin < 3) || isempty(C)
   % All pairwise contrasts
   np = length(p);
   nk = nchoosek(1:np,2);
   nContrasts = size(nk,1);
   C = zeros(nContrasts,np);
   ind = sub2ind(size(C),(1:nContrasts)',nk(:,1));
   C(ind) = 1;
   ind = sub2ind(size(C),(1:nContrasts)',nk(:,2));
   C(ind) = -1;
else
   % Need to check validity of contrasts
   nContrasts = size(C,1);
end
if (nargin < 4) || isempty(alpha)
   alpha = 0.05;
end
if nargin < 5 || isempty(tail)
    tail = 0;
elseif ischar(tail) && (size(tail,1)==1)
    tail = find(strncmpi(tail,{'left','both','right'},length(tail))) - 2;
end

p = p(:)';
pd = C*p'; % Proportion differences

if nContrasts > 1
   for i = 1:nContrasts
      [h(i,1),pval(i,1),z(i,1)] = compPropSameSurvey(p,n,C(i,:),alpha,tail,flag,false);
   end
else
   p1 = sum(p(C==1));
   p2 = sum(p(C==-1));
   z = (p1-p2) / sqrt((p1+p2)/n);
   
   if tail == 0 % two-tailed test
      pval = 2*normcdf(-abs(z),0,1);
   elseif tail == 1 % right one-tailed test
      pval = normcdf(-z,0,1);
   elseif tail == -1 % left one-tailed test
      pval = normcdf(z,0,1);
   else
      error('compPropSameSurvey:BadTail',...
         'TAIL must be ''both'', ''right'', or ''left'', or 0, 1, or -1.');
   end
   h = pval <= alpha;
end

% Confidence intervals
if (nargout >= 3) || verbose
   % covariance matrix of p (see Berry & Hurtado (1994), pg. 3)
   V = (diag(p) - p'*p)/n;
   se = sqrt(C*V*C');
   switch lower(flag)
      case 's' % simultaneous (Bonferonni-corrected)
         z0 = norminv(1-alpha/(2*nContrasts));
      otherwise % individual
         z0 = norminv(1-alpha/2);
   end
   
   ci(:,1) = pd - z0*diag(se);
   ci(:,2) = pd + z0*diag(se);
end
% c*V*c' using loops
% 
%ncat = length(p);
% for i = 1:ncat
%    for j = 1:ncat
%       if i == j
%          V(i,j) = p(i)*(1-p(i));
%       else
%          V(i,j) = -p(i)*p(j);
%       end
%    end
% end
% V = V/n;

if verbose || (nargout==0)
   fprintf('\n')
   fprintf('%g comparisons from the same survey\n',nContrasts);
   fprintf('n = %g\n',n);
   switch lower(flag)
      case 's'
         str = 'simultaneous (Bonferroni)';
      otherwise
         str = 'individual';
   end
   fprintf('CI type = %s\n',str);
   minlen = length('  contrast ');
   len = max(minlen,length(num2str(C(1,:)))+1);
   fprintf('%s------------------------------------------------------\n',repmat('-',len,1));
   fprintf('  contrast %s|   Diff.     z      P>|z|     [%g%% Conf. Interval]\n',repmat(' ',len-minlen,1),100*(1-alpha));
   fprintf('%s+-----------------------------------------------------\n',repmat('-',len,1));
   for i = 1:nContrasts
      str = num2str(C(i,:));
      pad = repmat(' ',1,max(0,minlen-length(str)-1));
      if (pval(i)<alpha) & (~((ci(i,1)<=0)&(ci(i,2)>=0)))
         fprintf('%s%s | %+1.4f   %+1.2f   %+1.4f*    %+1.4f    %+1.4f*\n',pad,str,pd(i),z(i),pval(i),ci(i,1),ci(i,2));
      elseif (pval(i)<alpha)% & ((ci(i,1)<=0)&(ci(i,1)>=0))
         fprintf('%s%s | %+1.4f   %+1.2f   %+1.4f*    %+1.4f    %+1.4f\n',pad,str,pd(i),z(i),pval(i),ci(i,1),ci(i,2));
      elseif ~((ci(i,1)<=0)&(ci(i,2)>=0))
         fprintf('%s%s | %+1.4f   %+1.2f   %+1.4f     %+1.4f    %+1.4f*\n',pad,str,pd(i),z(i),pval(i),ci(i,1),ci(i,2));
      else
         fprintf('%s%s | %+1.4f   %+1.2f   %+1.4f     %+1.4f    %+1.4f\n',pad,str,pd(i),z(i),pval(i),ci(i,1),ci(i,2));
      end
   end
   fprintf('%s------------------------------------------------------\n',repmat('-',len,1));
end


