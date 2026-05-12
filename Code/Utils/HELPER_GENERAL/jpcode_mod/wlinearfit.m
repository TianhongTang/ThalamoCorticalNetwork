function	[ A0, A1, Yout ]=wlinearfit(X,Y,W);
%function	[ output ]=wlinearfit(X,Y,W);
if nargin ~= 3
   error('Usage: 	[ A0, A1, Yout ]=wlinearfit(X,Y,W)')
end

xsize=size(X);
ysize=size(Y);
wsize=size(W);
if ( ndims(X) ~= 2 ) | (ndims(Y) ~= 2) | (ndims(W) ~= 2)...
      | (xsize(1) > 1 & xsize(2) > 1 ) ...
      | ((ysize(1) ~=  1 & ysize(2) ~= 1)) ...
      | ((wsize(1) ~= 1 ) & (wsize(2) ~= 1))
   error('Does not work on multidim arraies')
end
if xsize(1) ~= 1
   X=X';
end
if ysize(1) ~= 1
   Y=Y';
end
if wsize(1) ~= 1
   W=W';
end

sumW = sum(W);
A1 = ( sumW.*sum(X.*Y.*W) - sum(W.*X).*sum(W.*Y) )   / ( sum(W)*sum(W.*X.^2) - ((sum(W.*X)).^2) );
A0 = ( sum(Y) - A1*sum(X) ) / (max(xsize)) ;

Yout = A0 + A1.*X;