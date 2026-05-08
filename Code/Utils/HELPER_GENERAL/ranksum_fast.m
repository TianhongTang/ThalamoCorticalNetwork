function [p] = ranksum_fast(x,y,tail,method)
% [p] = ranksum_fast(x,y,tail,method)
% 2019-03-17 ESBM version of 'ranksum' that is simpler & faster
% - doesn't take tons of time on parameter checking each time it is called
%  (both for this function, and by inlining & simplifying its call to 
%   'tiedrank')
% - only has a single output (the p-value)
%
% tail: 'both', 'left', or 'right' 
% (default: 'both')
%
% method: 'exact' or 'approximate' 
% (default: 'exact' for small samples, 'approximate' otherwise)

% Check that x and y are vectors
if ~isvector(x) || ~isvector(y)
   error(message('stats:ranksum:InvalidData'));
end

% Remove missing data
x = x( ~isnan(x) );
y = y( ~isnan(y) );
nx = numel(x);
ny = numel(y);
ns = min(nx,ny);

if nx == 0 || ny == 0
	error(message('stats:signrank:NotEnoughData'));
end

% Set value for 'tail'
if nargin < 3 || isempty(tail)
    tail = 'both';
end

% Set value for 'method'
if nargin < 4 || isempty(method)
   if (ns < 10)  &&  ((nx+ny) < 20)
      method = 'exact';
   else
      method = 'approximate';
   end
end

% Set computational 'technique'
switch method
	case 'approximate'
		technique = 'normal_approximation';
	case 'exact'
		if (nx+ny) < 10
            technique = 'full_enumeration';
        else
            technique = 'network_algorithm';
        end
    otherwise
        error('unknown method, expected "exact" or "approximate"');
end

%      %      %      %      %      %      %      %      %      %
% Calculations for Rank Sum Test

x = x(:);   % ensure columns
y = y(:);
if nx <= ny
   smsample = x;
   lgsample = y;
   same_order = true;
else
   smsample = y;
   lgsample = x;
   same_order = false;
end

% Compute the rank sum statistic based on the smaller sample

[ranks,tieadj] = tr([smsample; lgsample]);
srank = ranks(1:ns);
w = sum(srank);


switch technique
	case 'full_enumeration'
		allpos = nchoosek(ranks,ns);   % enumerate all possibilities
		sumranks = sum(allpos,2);
		np = size(sumranks, 1);
		
		switch tail
			case 'both'
				plo = sum( sumranks <= w) / np ;
				phi = sum( sumranks >= w) / np ;
				p_tail = min(plo,phi);
				p = min(2*p_tail, 1);   % 2-sided, p>1 means middle is double-counted
				
			case 'right'
				switch same_order
					case true
						p = sum( sumranks >= w) / np ;
					case false
						p = sum( sumranks <= w) / np;
				end
				
			case 'left'
				switch same_order
					case true
						p = sum( sumranks <= w) / np ;
					case false
						p = sum( sumranks >= w) / np;
                end
            otherwise
                error('unknown setting of tail, expected "both", "right", or "left"');
				
		end
		
		%     %     %     %     %     %     %      %      %      %
		
	case 'network_algorithm'
		[p_net, pvals] = exactprob(smsample, lgsample, w);
		
		if any(isnan(p_net)) || any(isnan(pvals))
			warning(message('stats:ranksum:NanResult'));
			p = NaN;
			
		else
			switch tail
				case 'both'   % two-tailed test
					p = min(2*p_net, 1);   % p>1 means the middle is double-counted
					
				case 'right'   % right-tail test
					switch same_order
						case true
							p =  pvals(2) + pvals(3);
						case false
							p = pvals(2) + pvals(1);
					end
					
				case 'left'   % left-tail test
					switch same_order
						case true
							p =  pvals(2) + pvals(1);
						case false
							p = pvals(2) + pvals(3);
                    end
                otherwise
                    error('unknown setting of tail, expected "both", "right", or "left"');
					
			end   % conditional on 'tail'
		
		end
		
		%     %     %     %     %     %     %      %      %      %
				
	case 'normal_approximation'
		wmean = ns*(nx + ny + 1)/2;
		tiescor = 2 * tieadj / ((nx+ny) * (nx+ny-1));
		wvar  = nx*ny*((nx + ny + 1) - tiescor)/12;
		wc = w - wmean;

		% compute z-value, including continuity correction
		switch tail
			case 'both'
				z = (wc - 0.5 * sign(wc))/sqrt(wvar);
				if ~same_order
					z = -z;
				end
				p = 2*normcdf(-abs(z));
				
			case 'right'
				if same_order
					z = (wc - 0.5)/sqrt(wvar);
				else
					z = -(wc + 0.5)/sqrt(wvar);
				end
				
				p = normcdf(-z);
				
			case 'left'
				if same_order
					z = (wc + 0.5)/sqrt(wvar);
				else
					z = -(wc - 0.5)/sqrt(wvar);
				end
				
				p = normcdf(z);
            otherwise
                error('unknown setting of tail, expected "both", "right", or "left"');
		end

end   % conditional on 'technique'

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [p1, pvals] = exactprob(x,y,w)
%EXACTPROB Exact P-values for Wilcoxon Mann Whitney nonparametric test
%   [P1,PVALS]=EXACTPROB(X,Y,W) computes the p-value P for the test
%   statistic W in a Wilcoxon-Mann-Whitney nonparametric test of the
%   hypothesis that X and Y come from distributions with equal medians.

% Create a contingency table with element (i,j) indicating how many
% times u(j) appears in sample i
u = unique([x(:); y(:)]);
t = zeros(2,length(u));
t(1,:) = histc(x,u)';
t(2,:) = histc(y,u)';

% Compute weights for wmw test
colsum = sum(t,1);
tmp = cumsum(colsum);
wts = [0 tmp(1:end-1)] + .5*(1+diff([0 tmp]));

% Compute p-value using network algorithm for contingency tables
[p1, pvals] = statctexact_fast(t,wts,w);


end

% --------------------------------
function [r,tieadj] = tr(x)
%TR Local tiedrank function to compute results for one column
%
% ESBM: modified to always behave as if:
%  tieflag = false
%  bidirectional = 0
%  epsx = zeros(size(x))
% and to assume:
%  x is a column vector (since we assure this as input)
%  there are no NaNs (since we already removed them)

% Sort
[sx, rowidx] = sort(x);
xLen = numel(x);

% Use ranks counting from low end
ranks = (1:xLen)';

tieadj = 0;
if isa(x,'single')
   ranks = single(ranks);
   tieadj = single(tieadj);
end

% Adjust for ties.  Avoid using diff(sx) here in case there are infs.
ties = sx(1:xLen-1) >= sx(2:xLen);
tieloc = [find(ties); xLen+2];
maxTies = numel(tieloc);

tiecount = 1;
while (tiecount < maxTies)
    tiestart = tieloc(tiecount);
    ntied = 2;
    while(tieloc(tiecount+1) == tieloc(tiecount)+1)
        tiecount = tiecount+1;
        ntied = ntied+1;
    end

    tieadj = tieadj + ntied*(ntied-1)*(ntied+1)/2;
    
    % Compute mean of tied ranks
    ranks(tiestart:tiestart+ntied-1) = ...
                  sum(ranks(tiestart:tiestart+ntied-1)) / ntied;
    tiecount = tiecount + 1;
end

% Broadcast the ranks back out
r(rowidx) = ranks;

end