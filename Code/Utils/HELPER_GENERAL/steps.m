function [y,x] = steps(x,y,col,fl,wid)

if nargin < 1
 error('USAGE: steps(y) or steps(x,y,col,fill)');
end

if ~exist('fl') | isempty(fl), fl = 0; end
if ~exist('col') | isempty(col), col = 'r'; end
if ~exist('wid') | isempty(wid), wid = 0.5; end

if nargin == 1
  x = [0 0 repeat(1:length(x),2)];
  y = [0 repeat(x,2) 0];
else
  x = [x(1) repeat(x,2)];
  y = [repeat(y,2) y(end)];
end

if nargout == 0
  if length(col) == 1
    if fl  
      fill(x,y,col,'LineWidth',wid);
    else 
      plot(x,y,col,'LineWidth',wid);
    end
  else
    if fl  
      fill(x,y,col,'LineWidth',wid) 
    else 
      plot(x,y,'Color',col,'LineWidth',wid);
    end
  end
end