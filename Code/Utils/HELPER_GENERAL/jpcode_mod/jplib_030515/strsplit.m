function [strarr, nelems] = strsplit(str, delimiter)
%
%%   splits string into array of tokens
%%   [strarr, nelems] = strsplit(str, delimiter)
%%   [strarr, nelems] = strsplit(str)
%%   
%%   effects: cleans leading and trailing spaces of each element
%%
%%   JP Apr 2000

narg = nargin;
error(nargchk(1,2,narg));

istr = 0;	
rem = deblank(strjust(str,'left'));
while ( ~isempty(rem) )
      istr = istr + 1;
      if (narg == 1)	
         [token,rem] = strtok( rem );
      else
         [token,rem] = strtok( rem, delimiter);
      end
      strarr(istr,:) = {deblank(strjust(token,'left'))};   
      rem = deblank(strjust(rem,'left'));
      %%%fprintf('<%s>\n',strarr(istr,:))
end
nelems = istr;

%--------------------------------------------------------------

