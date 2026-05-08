function nbytes=checkfilesize(filename, verbose)

%(from pfvdm) checkfilesize.m Usage: checkfilesize(filename)
%written by Pierre-Francois Van de Moortele

if nargin < 1 
	error('Usage : checkfilesize(filename,verbose)');
else if nargin >= 2
        f_verbose = verbose;
    else
        f_verbose = 0;
    end
end
thedir=dir(filename);
thesize=size(thedir);
if thesize(1) ~= 1 
   if f_verbose
       whos thedir
       thesize
       dir(filename);
       fprintf('problem : size of <%s>\n',filename);
   end
   nbytes = 0;
   return
end
if thedir(thesize(1)).isdir == 1
   if f_verbose
       fprintf('OOps... This is a Directory !\n');
   end
   nbytes = 0;
   return
end
nbytes=thedir(thesize(1)).bytes;
