function key(arg1)
% %  
% %
if (nargin >= 1)
    str = num2str(arg1);
else
    str = '';
end
pause on ;
fprintf(strcat(str,'==> ')); 
pause;