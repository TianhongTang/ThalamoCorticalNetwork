function cl = doConv(l,k)

  if isempty(k)
    cl = l;
    return;
  end
if ndims(l) == 2
  tmp = size(l);
  if (tmp(1) > tmp(2))     % allow for columnar matrix
    sz = tmp(2);
    l = l';
  else
    sz = tmp(1);
  end
  for i=1:sz
    cl(i,:) = convList(l(i,:),k);
  end
else
  cl = convList(l,k);
end

  
%--------------------------------------------------

function cl = convList(l,k)
  hlen    = floor(length(k)/2);
  llen    = length(l);	
  begmean = mean(l(1:hlen))*ones(1,hlen);
  endmean = mean(l(llen-hlen:llen))*ones(1,hlen);
  ltmp    = [begmean l endmean]; 	
  tmp    = conv(ltmp,k);
  start  = round(length(k)-1);
  finish = start+length(l)-1;
  cl     = tmp([start:finish]); 
  