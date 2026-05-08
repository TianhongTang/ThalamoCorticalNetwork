function COORD = LCMrdCoord(dir, filenum, optin);
% %      COORD = LCMrdCoord(dir, filenum, optin)
% %
% %  --- reads LCModel output file  ---
% % 
% %  INPUT:
% % 	dir
% % 	filenum 
% % 
% %  default options:
% %     opt(  'DIRAPP','',    append 'DIRAPP' to dir 
% %           'TWODIM',0,
% %           'VERBOSE',1)    file info & plot
% %
% %  OUTPUT:
% %       COORDstruct = 
% % 	 	coord.metab
% % 	 	coord.conc  
% %  	 	coord.relconc 
% % 	 	coord.SD    
% % 	 	coord.ndata 
% % 	 	coord.data  
% % 	 	coord.mdata  
% % 	 	coord.header:   header1
% % 	 	coord.nconc 
% % 	 	coord.tconc :   header2
% % 	 	coord.tmisc :   header3
% % 	 	coord.tdiag :   header4
% % 	 	coord.tinput:   header5
% % 
% %  ORDER OF reading:
% %  	header1: general info, date
% % 	header2: concentration table
% %     header3: misc. output table
% % 	data: (NY, NDATABLOCKS)
% % 		1: ppm axis
% % 		2: phased data points
% % 		3: fit to the follow
% % 		4: background values
% % 	mdata: (NY, nconc)
% % 	       metabolite data, index assingment according coord.metab
% % 	header4: diagnostic table
% % 	header5: table of input changes
% % 		
% % 
% %  Tested for: PC (not Linux)
% %
% %  March 2003 -  Josef Pfeuffer
% %
FCTNAME = 'LCMrdCoord';

global STDPATH

STDPATHLOCAL = STDPATH.lcm; 
STDEXT  = '.fid';	% of VN data directory
LCMEXT  = '.lcm';	% of LCM data directory
LCMEXT4  = '.COORD';
FNAMEDEF = 'lcm';	% of LCModel raw data filename
TWODIMAPP = 'S'; 	% TwoDim appendix

NUMHLINES = 4;		% initial header lines
NDATABLOCKS = 4;	% ppm-axis, phased data, fit, background 
NELTDIAG = 10;		% minimum number of lines in tdiag

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.DIRAPP  = '';
dopt.TWODIM  = 0;
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

dirapp    = dopt.DIRAPP;
flagTwoDim  = dopt.TWODIM;      
flagVerbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

if (filenum > 0)     % stdpath handling
   file = sprintf('%s%s/%s', STDPATHLOCAL, num2str(dir), ...
      num2str(filenum));
else              % give full path as arg1
   file = sprintf('%s', num2str(dir));
end

% IF(strpos(dir, STDEXT) GT -1) THEN $
%    fdir = strmid(dir, 0, strpos(dir, STDEXT)) $
% ELSE IF(strpos(dir, LCMEXT) GT -1) THEN $
%    fdir = strmid(dir, 0, strpos(dir, LCMEXT)) $
% ELSE IF(strpos(dir, LCMEXT4) GT -1) THEN $
%    fdir = strmid(dir, 0, strpos(dir, LCMEXT4)) $
% ELSE $
%    fdir = dir
fdir = dir;

		%% --- whole path: +SLICEAPP+islice+LCMEXT+'/'
if narg == 2 
%     j = checkfile(STDPATH.lcm + fdir, /write, IS_DIR=is_dir)
%     if is_dir EQ 1 then $
%        dirstr = "/" $
%    else $
%       dirstr = "_"
    dirstr = '/';
    fdir = strcat(STDPATH.lcm, fdir, dirstr, sprintf('%d',filenum), dirapp);
else
    fdir = strcat(fdir, dirapp);
end
lcmfile = strcat(FNAMEDEF, LCMEXT4);

%%strsl = SLICEAPP + '*' 
strsl = '';
% j = findfile(fdir + strsl + LCMEXT + '/' + lcmfile, count=nslices)
% if (nslices LE 0) then begin
%    print,'Not found <' + fdir + strsl + LCMEXT + '/' + lcmfile + '> !'
%    return, 0
% endif 
% 
% if f_slices(0) EQ 0 then $
nslices = 1;
line = '';
ymax = 0;

numTwoDim = nslices;
if flagTwoDim == 0 
    numTwoDim = 1;
end
	
for itwodim=1:numTwoDim 		% ---- LOOP:  2D data set
    
    if flagTwoDim ~= 0
        strTwoDim = sprintf('%s%d', TWODIMAPP, itwodim);
    else
        strTwoDim = '';
    end
    fdirslice = strcat(fdir, strTwoDim, LCMEXT, '/');
    file = strcat(fdirslice, lcmfile);
    
    coord = LCMinitCoord;
    unit = fopen(file,'r');
    if unit < 0
        error(sprintf('Can not open <%s>', file));
    end
    if flagVerbose ~= 0
        fprintf('reading <%s>\n', file);
    end

			% --- read header1 ---------------------------
    for i=1:NUMHLINES
        coord.header(i,1) = { fgetl(unit) };
        if flagVerbose
            fprintf('%s\n', char(coord.header(i,1)));
        end
    end

  			% --- read header2 (tconc), determine NCONC ----------
    line = fgetl(unit);
    if flagVerbose
        fprintf('%s\n', line);
    end
    coord.nconc = strread(strtok(line),'%f') - 1;
    for i=1:coord.nconc+1
        coord.tconc(i,1) = { fgetl(unit) };
        if flagVerbose
            fprintf('%s\n', char(coord.tconc(i,1)));
        end
    end
    for i=2:coord.nconc+1
%         [a,b,c,d,e] = strread( char(coord.tconc(i,1)),'%f%f%c%f%s');
%         coord.conc(i-1,1) = a;
%         coord.SD(i-1,1) = b;
%         coord.relconc(i-1,1) = d;
%         coord.metab(i-1,1) = e;

		% dum = temporary dummy string
		[a,dum1] = strtok(char(coord.tconc(i,1)),' ');
		[b,dum2] = strtok(dum1,'%');
		[dum3,dum4] = strtok(dum2,' ');
		[dum5,dum6] = strtok(dum4,' ');
		if ~strcmp(dum5,strrep(dum4,' ',''))
            d = dum5;
            e = strtok(dum6,' ');
		else
            if ~isempty(findstr(dum5,'+'))
                [d,dum7] = strtok(dum5,'+');
                [e,dum8] = strtok(dum7,'+');
            elseif ~isempty(findstr(dum5,'-'))
                [d,dum7] = strtok(dum5,'-');
                [e,dum8] = strtok(dum7,'-');
            else
                error(sprintf('%s -> reading metabolite concentrations failed'));
            end
		end
		coord.conc(i-1,1) = str2double(a);
		coord.SD(i-1,1) = str2double(b);
		coord.relconc(i-1,1) = str2double(d);
		coord.metab(i-1,1) = {e};
    end
     
 			% --- read header3 (tmisc) ---------------------------
    line = fgetl(unit);
    if flagVerbose
        fprintf('%s\n', line);
    end
    nmisc = strread(strtok(line),'%f');
    for i=1:nmisc
        coord.tmisc(i,1) = { fgetl(unit) };
        if flagVerbose
            fprintf('%s\n', char(coord.tmisc(i,1)));
        end
    end
 
 			% --- read data blocks ----------------------
    for idata=1:NDATABLOCKS 
        %   --- determine NY (coord.ndata)
        line = fgetl(unit);
        if flagVerbose >= 2
            fprintf('%s\n', line);
        end
        if idata == 1    % --- read NY from line/ init dataarr
            coord.ndata = strread(strtok(line),'%f');
            coord.data = zeros(coord.ndata, NDATABLOCKS);
        end
        [bdata,count] = fscanf(unit, '%f ', coord.ndata);   % !!! ' ' to get rest of line !!!
        if count ~= coord.ndata
            error(sprintf('%s: read format error <%s/%s> read\n', FCTNAME, count, coord.ndata));
        else
            coord.data(:,idata) = bdata(:);
        end
        %numcols = length( strread(line,'%f') );
        %numlines = fix( coord.ndata / numcols);
        %restelements = mod(coord.ndata, numcols);
        
        flagNY = 0;
        while ~flagNY    % read rest of line after last float
            unitpos = ftell(unit);
            line = fgetl(unit);
            NYstr = strtok(line);
            if strcmp(NYstr,'NY') | findstr(line,'Conc') | findstr(line,'line')
                flagNY = 1;
                fseek(unit,unitpos,'bof');
            end
        end    
    end   % idata
    
 			% --- read mdata blocks ----------------------
    line = fgetl(unit);
    if flagVerbose  >= 2
        fprintf('%s\n', line);
    end
    coord.mdata = zeros(coord.ndata, coord.nconc);
    while isempty( findstr(line,'diagnostic') )
        curMetabolite = strtok(line);
        [bdata,count] = fscanf(unit, '%f ', coord.ndata);   % !!! ' ' to get rest of line !!!
        if count ~= coord.ndata
            error(sprintf('%s: read format error <%s/%s> read\n', FCTNAME, count, coord.ndata));
        else
            midx = 0;
            for i=1:length(coord.metab)
                if strcmp(char(coord.metab(i)), curMetabolite)
                    midx = i;
                end
            end
            if midx > 0
                coord.mdata(:,midx) = bdata(:);
            else
                error(sprintf('%s: metabolite assignment index not found <%s> read\n', FCTNAME, curMetabolite));
            end
        end
        
        line = fgetl(unit);
        if flagVerbose
            fprintf('%s\n', line);
        end
    end
    
 			% --- read header4 (tdiag) ---------------------------
    ndiag = strread(strtok(line),'%f');
    ndiag = max([ndiag 1]);    % read at least 1 line 
    for i=1:ndiag
        coord.tdiag(i,1) = { fgetl(unit) };
        if flagVerbose
            fprintf('%s\n', char(coord.tdiag(i,1)));
        end
    end

 			% --- read header5 (tinput) ---------------------------
    line = fgetl(unit);
    if flagVerbose
        fprintf('%s\n', line);
    end
    ninput = strread(strtok(line),'%f');
    for i=1:ninput
        coord.tinput(i,1) = { fgetl(unit) };
        if flagVerbose
            fprintf('%s\n', char(coord.tinput(i,1)));
        end
    end

    fclose(unit);

    if flagVerbose
        d = coord.data;
        m = coord.mdata;
        subplot(2,1,1)
        plot(d(:,1),d(:,2))
        set(gca,'xDir','reverse')
        xlim([min(coord.data(:,1)) max(coord.data(:,1))]);
        ylim([min(coord.data(:,2)) max(coord.data(:,2))]);
        xlabel('ppm')
        ylabel('Spectrum / a.u.')
        title(strcat(' ',char(coord.header(2))))
        hold on
        plot(d(:,1),d(:,3),'k')
        plot(d(:,1),d(:,4),'k-')
        hold off
        subplot(2,1,2)
        plot(d(:,1),d(:,2)-d(:,3),'k')
        set(gca,'xDir','reverse')
        xlim([min(coord.data(:,1)) max(coord.data(:,1))]);
        xlabel('ppm')
        ylabel('Residual')
    end
end 	% for: --- itwodim

if nargout,
	COORD = coord;
end

return
%------------------------------------------------------------
function coord = LCMinitCoord;
% %  coord = LCMinitCoord
% %
% %  LCModel
% %  inits coord struct
% %
% %
FCTNAME = 'LCMinitCoord';

	coord.metab   = {''};
	coord.conc    = 0;
 	coord.relconc = 0;
	coord.SD      = 0;
	coord.ndata   = 0;
	coord.data    = 0;
	coord.mdata   = 0;
	coord.header  = {''};
	coord.nconc   = 0;
	coord.tconc   = {''};
	coord.tmisc   = {''};
	coord.tdiag   = {''};
	coord.tinput  = {''};
% --- end of COORD field definition !
end
% --------------------------------------------------------------
function LCMrdCoordTest(dir, filenum, optin);

global testdat

if nargin < 1
    COORD = LCMrdCoord('G02.gM1', 112)
end
    
% --------------------------------------------------------------
