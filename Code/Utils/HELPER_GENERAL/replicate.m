function l2 = replicate(l1,n,one_dim)
%
% FUNCTION
% l2 = replicate(l1,n)
%
% ARGS
% l1 = input list
% n  = number of times to replicate     i.e.  12341234123412341234
% 
% see also: repeat


  if ~exist('one_dim') | isempty(one_dim)
       one_dim = 0;
  end 
  if nRealDims(l1) == 1 & ~one_dim
    l1 = l1(:);
    % make a 2-d array
    l2 = zeros(length(l1),n); %allocate
    for i = 1:n
      l2(:,i) = l1;
    end
  else
    l2 = [];
    for i = 1:n
      l2 = [l2 l1];
    end
  end
    
function n = nRealDims(array)
  if ndims(array) > 2
    n = ndims(array);
  else
    sz = size(array);
    if sz(1) == 1 & sz(2) == 1
      n = 0;
    elseif sz(1) == 1 | sz(2) == 1
      n = 1;
    else
      n = 2;
    end
  end


