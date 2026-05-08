function s = opt(varargin)
%
%%   builds struct from arguments
%%   optstruct = opt('OPTION1',var1,'OPTION2',var2, ...)
%%
%%   effects: fielddescriptors are converted to UPPERCASE
%%
%%   JP Apr 2000

tagfield = 'STRUCTNAME';
tagstr = 'optstruct';

len_arg = length(varargin);
s = [];
s = setfield(s,tagfield,tagstr);      % init s - return to call w/o argument
for ilen = 1:2:len_arg-1
   if (~ischar(varargin{ilen}))
      varargin
      error('opt(): descriptor is not a string !');
   end
   fielddescriptor = upper(num2str( varargin{ilen} ));
   s = setfield(s, fielddescriptor, varargin{ilen+1});
end

%--------------------------------------------------------------
