function optout = setopt(optall, optin)
%
%%   updates default optall struct with input optin
%%   optout = setopt(optall, optin)
%%
%%   effects: warning on unkown option field
%%            ignores tagfield
%%
%%   JP Apr 2000

narg = nargin;
error(nargchk(2,2,narg));

tagfield = 'STRUCTNAME';
tagstr = 'optstruct';

optout = optall;
if (isempty(optin))
   return
end
descin = fieldnames(optin);
l_descin = length(descin);

for il=1:l_descin
   if ( strcmp(descin(il), tagfield) &  ...
	strcmp(getfield(optin,descin{il}), tagstr) )
      %% do nothing
      %% fprintf('field <%s> contains <%s>',descin{il}, ...
      %%             getfield(optin,descin{il}));
   elseif (  isfield(optall, descin(il)) )
      optout = setfield(optout, descin{il}, getfield(optin,descin{il}));
   else
      errmsg = sprintf('! option <%s> ignored !',descin{il});
      warning(errmsg);
   end
end % for

%--------------------------------------------------------------
