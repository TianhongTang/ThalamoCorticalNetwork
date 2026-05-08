function	[abcisse, mypowerabs, myzeropowerabs, mypowernorm, myzeropowernorm, FFTDATA]=autopower(MYFILENAME,skipbegin,skipend,TR_second,dimx_or_writewhat,dimy,dimz,nb_reptot)

%(pfvdm)  usage: autopower(filename,skipbegin,skipend,TR_second [,writewhat]) 
%
%values for [writewhat]: (put as many as you want, including the brackets, e.g [1,2,4,6])
% 1: spectrum intensities
% 2: normalized spectrum intensities
% 3: power spectrum 
% 4: normalized power spectrum 
% 5: phase of spectrum
% 6: matlab file with complex FFT product
%
% IF ARGUMENT [writewhat] NOT SUPPLIED, ALL FILES ARE WRITTEN
%
%WARNING: currently, powernorm is multiplied by 10e6 !


%pfvdm [abcisse mypowerabs myzeropowerabs mypowernorm myzeropowernorm FFTDATA]...
% ...=autopower(MYFILENAME,skipbegin,skipend,TR_second [,dimx,dimy,dimz,skipbegin,nb_reptot])
% files written:
%	- power abs (puissance au carre)
%	- power norm (puissance au carre normalisee sur freq nulle)
%	- FFTDATA (phase !)
%	- spectrum abs (amplitude spectre (not puissance))
%	- cpxfftata (complex fourier Transform , from freq nulle to nyquist, MATLAB file)

%global DATA FFTDATA INDEXFREQMAX ABCISSE PWNORM MY_PW_ZERO
%global MYPOWERABS MYZEROPOWERABS MYSPECTRUMABS MYNORMPOWERABS
globalvar_power
global myshortfilename scale_factor
global DATA NB_VOL_KEEP 

DEBUG = 0;
BYTEORDER = 'ieee-be';

   % --- defs
%precind=    1       2        3       4          5       6
precStr = {'';   'int16'; 'int32'; 'float32'; '';      ''};
typeStr = {'BYTE'; 'WORD';'LWORD'; 'REAL'; 'COMPLEX'; 'TBD'};
precind = 4;

fullfilename=sdtname(MYFILENAME);
if nargin == 4
	%readstimulate3(fullfilename,skipbegin,skipend);
elseif nargin == 5
	writewhat=dimx_or_writewhat
	%readstimulate3(fullfilename,skipbegin,skipend);
elseif nargin == 9
	%readstimulate3(fullfilename,dimx,dimy,dimz,nb_reptot,skipbegin,skipend);
else
	disp 'usage: [abcisse MYPOWERABS mypowernorm FFTDATA]='
	disp ' '
	disp '        autopower2(filename,skipbegin,skipend,TR_second [,writewhat]) '
	disp ' '
	disp '	      or, instead of [writewhat]: [,dimx,dimy,dimz,skipbegin,nb_reptot]'
	disp ' '
	disp 'values for [writewhat]: (put as many as you want, including the brackets, e.g [1,2,4,6])'
	disp ' 1: spectrum intensities'
	disp ' 2: normalized spectrum intensities'
	disp ' 3: power spectrum '
	disp ' 4: normalized power spectrum '
	disp ' 5: phase of spectrum'
	disp ' 6: matlab file with complex FFT product'
	disp ' '
	disp 'WARNING: not properly handled in .spr currently:'
	disp ' - 4th dim of fov, interval, extent'
	disp ' - displayRange'
	disp 'WARNING: currently, powernorm is multiplied by 10e6 !'
	return
end

    % adapted by JP
DATA = STIMrdSdt( MYFILENAME );
if skipbegin >= 1 
	DATA(:,:,:,[1:skipbegin])=[];
end
if skipend >= 1
	DATA(:,:,:,[end-skipend+1:end])=[];
end
s_data = size(DATA);
NB_VOL_KEEP = s_data(4)

figure(3)
imagesc(DATA(:,:,1,1))
pause(1)
    % end JP

fprintf('<autopower> job starting on file :%s\n',myshortfilename);

%---analyze writewhat---
absspectrumflag=0;
normspectrumflag=0;
abspowerflag=0;
normpowerflag=0;
phaseflagflag=0;
cpxmatlabflag=0;

if nargin == 5 
	if (intersect(writewhat,1))
		absspectrumflag=1;
	end
	if (intersect(writewhat,2))
		normspectrumflag=1;
	end
	if (intersect(writewhat,3))
		abspowerflag=1;
	end
	if (intersect(writewhat,4))
		normpowerflag=1;
	end
	if (intersect(writewhat,5))
		phaseflagflag=1;
	end
	if (intersect(writewhat,6))
		cpxmatlabflag=1;
	end
elseif nargin == 4
	absspectrumflag=1;
	normspectrumflag=1;
	abspowerflag=1;
	normpowerflag=1;
	phaseflagflag=1;
	cpxmatlabflag=1;
end

%---start program----

fprintf('computing fft...\n')
if DEBUG, tic, end
fftdimtpos;
if DEBUG, toc, end
fprintf('converting abcisse to freq ...\n')
freqpos(INDEXFREQMAX-1,TR_second);
fprintf('computing powerabs...\n')
if DEBUG, tic, end
MYSPECTRUMABS=abs(FFTDATA);
MYPOWERABS=MYSPECTRUMABS.^2;
%MYPOWERABS=abs(FFTDATA).^2;
if DEBUG, toc, end
if normpowerflag
	fprintf('computing powernorm...\n')
	if DEBUG, tic, end
	powernorm;
	if DEBUG, toc, end
end
if normspectrumflag
	fprintf('computing spectrumnorm...\n')
	if DEBUG, tic, end
	spectrumnorm;
	if DEBUG, toc, end
end
fprintf('\nEnd of computing for %s \n',fullfilename);

%IGNORE DC COMPONENT (already done in powernorm, not in abs)
%scaling for display in stimulate
%--------------------------------
fprintf('...just reshaping some matrix..\n');
MYZEROPOWERABS=MYPOWERABS(:,:,:,1);
%MYPOWERABS(:,:,:,1)=[];

if normpowerflag
	scale_factor=1e6;
	PWNORM=PWNORM*scale_factor;
	MY_PW_ZERO=MY_PW_ZERO*scale_factor;
end


%writing files
%---------------
[pa,na,ex,ve]=fileparts(fullfilename);
myshortfilename=fullfile(pa,na);
nativeheader=sprname(MYFILENAME);
%fidnativeheader=fopen(nativeheader,'r');

%--some of the parameters from native header
%------------------------------------------
fidnativeheader=fopen(nativeheader,'r');
done=0;

fovstr='';
originstr='';
extentstr='';
intervalstr='';
myorient='';

while( done == 0 )
 line = fgetl(fidnativeheader);
 if (line == -1)
   done = 1;
  else
   [attr, val] = strtok(line,':');
   if isequal(attr,'fov'), fovstr = strtok(val,':'); end
   if isequal(attr,'origin'), originstr = strtok(val,':'); end
   if isequal(attr,'extent'), extentstr = strtok(val,':'); end
   if isequal(attr,'interval'), intervalstr = strtok(val,':'); end
   if isequal(attr,'sdtOrient'), myorient = strtok(val,':'); end
  end
end
myfov=str2num(fovstr);
myorigin=str2num(originstr);
myextent=str2num(extentstr);
myinterval=str2num(intervalstr);
if 0
	fprintf('output fov: %g %g %g\n',myfov(1),myfov(2),myfov(3));
	fprintf('output origin: %g %g %g\n',myorigin(1),myorigin(2),myorigin(3));
	fprintf('output extent: %g %g %g\n',myextent(1),myextent(2),myextent(3));
	fprintf('output interval: %g %g %g\n',myinterval(1),myinterval(2),myinterval(3));
	fprintf('output orient: %s\n',myorient);
end
fclose(fidnativeheader);

	%temporary shorcut
	%-----------------
	tempo_orient_string='';

	%---power abs files names
if  abspowerflag
%'abs' is confusing with 'power', let's take it out of the filename
mypowerabsfile=[myshortfilename,'pw.sdt'];
mypowerabsheader=[myshortfilename,'pw.spr'];

%'abs' is confusing with 'power', let's take it out of the filename
mypowerabsfile_zero=[myshortfilename,'pw_zero.sdt'];
mypowerabsheader_zero=[myshortfilename,'pw_zero.spr'];

	%----powerabs_--
dimpowerabs=size(MYPOWERABS);
ndimpowerabs=ndims(MYPOWERABS);
mydimx=dimpowerabs(1);
mydimy=dimpowerabs(2);
if ndimpowerabs==3 
	mydimz=1;
	mydimt=dimpowerabs(3);
else
	mydimz=dimpowerabs(3);
	mydimt=dimpowerabs(4);
end
fprintf('\n..writing %s\n',mypowerabsfile);
fidabsf=fopen(mypowerabsfile,'w',BYTEORDER);
fwrite(fidabsf, MYPOWERABS, char(precStr(precind)));
minPOWERABS=min(MYPOWERABS(:));
maxPOWERABS=max(MYPOWERABS(:));
fclose(fidabsf);
fidabsheader=fopen(mypowerabsheader,'w');
fidnativeheader=fopen(nativeheader,'r');
while ~feof(fidnativeheader)
	ss=fgetl(fidnativeheader);
	[ attr, val ]=strtok(ss,':');
	if ss(1,1:4) == 'dim:'
		fprintf(fidabsheader,'dim: %d %d %d %d\n',mydimx,mydimy,mydimz,mydimt);
		fprintf('dim: %d %d %d %d\n',mydimx,mydimy,mydimz,mydimt);
	elseif ss(1,1:8) == 'fidName:'
		fprintf(fidabsheader,'fidName: %s\n',mypowerabsfile);
		fprintf('fidName: %s\n',mypowerabsfile);
	elseif isequal(attr,'fov')
		fprintf(fidabsheader,'fov: %g %g %g %g\n',myfov(1),myfov(2),myfov(3),mydimt);
		fprintf('fov: %g %g %g %g\n',myfov(1),myfov(2),myfov(3),mydimt);
	elseif isequal(attr,'extent')
		fprintf(fidabsheader,'extent: %g %g %g %g\n',myextent(1),myextent(2),myextent(3),mydimt-1);
		fprintf('extent: %g %g %g %g\n',myextent(1),myextent(2),myextent(3),mydimt-1);
	elseif isequal(attr,'displayRange')
		fprintf(fidabsheader,'displayRange: %g %g\n',minPOWERABS,maxPOWERABS);
		fprintf('displayRange: %g %g\n',minPOWERABS,maxPOWERABS);
	else
		fprintf(fidabsheader,'%s\n',ss);
		fprintf('%s\n',ss);
	end
	%if ss(1,1:9) == 'sdtOrient'
	%	tempo_orient_string=ss;		
	%	[mybegin myend]=strtok(tempo_orient_string,':');
	%	tempo_orient_string=myend(2:end);
	%end
end
fclose(fidabsheader);
fclose(fidnativeheader);

fprintf('\n..writing %s\n',mypowerabsfile_zero);
fidabsf_zero=fopen(mypowerabsfile_zero,'w',BYTEORDER);
fwrite(fidabsf_zero,MYZEROPOWERABS,char(precStr(precind)));
fclose(fidabsf_zero);
fidabsheader_zero=fopen(mypowerabsheader_zero,'w');
fidnativeheader=fopen(nativeheader,'r');
while ~feof(fidnativeheader)
        ss=fgetl(fidnativeheader);
	[ attr, val ] = strtok(ss,':');
	if ss(1,1:4) == 'dim:'
		fprintf(fidabsheader_zero,'dim: %d %d %d\n',mydimx,mydimy,mydimz);
		fprintf('dim: %d %d %d\n',mydimx,mydimy,mydimz);
	elseif ss(1,1:8) == 'fidName:'
		fprintf(fidabsheader_zero,'fidName: %s\n',mypowerabsfile_zero);
		fprintf('fidName: %s\n',mypowerabsfile_zero);
	elseif isequal(attr,'numDim')
		fprintf(fidabsheader_zero,'numDim: %d\n',3);
		fprintf('numDim: %d\n',3);
	elseif isequal(attr,'fov')
		fprintf(fidabsheader_zero,'fov: %g %g %g\n',myfov(1),myfov(2),myfov(3));
		fprintf('fov: %g %g %g\n',myfov(1),myfov(2),myfov(3));
	elseif isequal(attr,'extent')
		fprintf(fidabsheader_zero,'extent: %g %g %g\n',myextent(1),myextent(2),myextent(3));
		fprintf('extent: %g %g %g %g\n',myextent(1),myextent(2),myextent(3));
	elseif isequal(attr,'origin')
		fprintf(fidabsheader_zero,'origin: %g %g %g\n',myorigin(1),myorigin(2),myorigin(3));
		fprintf('origin: %g %g %g\n',myorigin(1),myorigin(2),myorigin(3));
	elseif isequal(attr,'interval')
		fprintf(fidabsheader_zero,'interval: %g %g %g\n',myinterval(1),myinterval(2),myinterval(3));
		fprintf('interval: %g %g %g\n',myinterval(1),myinterval(2),myinterval(3));
	else
		fprintf(fidabsheader_zero,'%s\n',ss);
		fprintf('%s\n',ss);
	end
end
fclose(fidabsheader_zero);
fclose(fidnativeheader);

end %end flag


	%---power norm files names
if normpowerflag
mypowernormfile=[myshortfilename,'pwnorm.sdt'];
mypowernormheader=[myshortfilename,'pwnorm.spr'];

mypowernormfile_zero=[myshortfilename,'pwnorm_zero.sdt'];
mypowernormheader_zero=[myshortfilename,'pwnorm_zero.spr'];

nativeheader=sprname(MYFILENAME);
fidnativeheader=fopen(nativeheader,'r');

	%----powernorm_--
dimpowernorm=size(PWNORM );
ndimpowernorm=ndims(PWNORM );
mydimx=dimpowernorm(1);
mydimy=dimpowernorm(2);
if ndimpowernorm==3 
	mydimz=1;
	mydimt=dimpowernorm(3);
else
	mydimz=dimpowernorm(3);
	mydimt=dimpowernorm(4);
end
fprintf('\n..writing %s\n',mypowernormfile);
fidnormf=fopen(mypowernormfile,'w',BYTEORDER);
fwrite(fidnormf,PWNORM,char(precStr(precind)));
fclose(fidnormf);
fidnormheader=fopen(mypowernormheader,'w');
while ~feof(fidnativeheader)
	ss=fgetl(fidnativeheader);
	[ attr, val ] = strtok(ss,':');
	if ss(1,1:4) == 'dim:'
		fprintf(fidnormheader,'dim: %d %d %d %d\n',mydimx,mydimy,mydimz,mydimt);
		fprintf('dim: %d %d %d %d\n',mydimx,mydimy,mydimz,mydimt);
	elseif ss(1,1:8) == 'fidName:'
		fprintf(fidnormheader,'fidName: %s\n',mypowernormfile);
		fprintf('fidName: %s\n',mypowernormfile);
	elseif isequal(attr,'extent')
		fprintf(fidnormheader,'extent: %g %g %g %g\n',myextent(1),myextent(2),myextent(3),mydimt-1);
		fprintf('extent: %g %g %g %g\n',myextent(1),myextent(2),myextent(3),mydimt-1);

	elseif isequal(attr,'fov')
		fprintf(fidnormheader,'fov: %g %g %g %g\n',myfov(1),myfov(2),myfov(3),mydimt);
		fprintf('fov: %g %g %g %g\n',myfov(1),myfov(2),myfov(3),mydimt);

	else
		fprintf(fidnormheader,'%s\n',ss);
		fprintf('%s\n',ss);
	end
end
fclose(fidnormheader);
fclose(fidnativeheader);
end % end flag 

	%---cpx fft product as matlab file----
if cpxmatlabflag
cpxfftdatafile=[myshortfilename,'cpxfftdata'];
%save cpxfftdatafile FFTDATA;
save(cpxfftdatafile,'FFTDATA');
fprintf('\n...saved complex file in matlab format: %s \n',cpxfftdatafile);
end % end flag

	%---FFTDATA complex file names
if phaseflagflag
fftdatafile=[myshortfilename,'fftdata.sdt'];
fftdataheader=[myshortfilename,'fftdata.spr'];

fftdatafile_zero=[myshortfilename,'fftdata_zero.sdt'];
fftdataheader_zero=[myshortfilename,'fftdata_zero.spr'];

	%----FFTDATA
dimfftdata=size(FFTDATA);
ndimfftdata=ndims(FFTDATA);
mydimx=dimfftdata(1);
mydimy=dimfftdata(2);
if ndimfftdata==3 
	mydimz=1;
	mydimt=dimfftdata(3);
else
	mydimz=dimfftdata(3);
	mydimt=dimfftdata(4);
end
fprintf('\n..writing %s\n',fftdatafile);
fidfftdata=fopen(fftdatafile,'w',BYTEORDER);
%BAD fwrite(fidfftdata,FFTDATA,'complex');
%good for complex fwrite(fidfftdata,[real(FFTDATA(:).');imag(FFTDATA(:).')],'float');
%currently let's write the phase
fwrite(fidfftdata,angle(FFTDATA),char(precStr(precind)));
fclose(fidfftdata);
fidfftdataheader=fopen(fftdataheader,'w');
fidnativeheader=fopen(nativeheader,'r');
while ~feof(fidnativeheader)
	ss=fgetl(fidnativeheader);
	[ attr, val ] = strtok(ss,':');
	if ss(1,1:4) == 'dim:'
		fprintf(fidfftdataheader,'dim: %d %d %d %d\n',mydimx,mydimy,mydimz,mydimt);
		fprintf('dim: %d %d %d %d\n',mydimx,mydimy,mydimz,mydimt);
	elseif ss(1,1:8) == 'fidName:'
		fprintf(fidfftdataheader,'fidName: %s\n',fftdatafile);
		fprintf('fidName: %s\n',fftdatafile);
	elseif ss(1,1:9) == 'dataType:'
		fprintf(fidfftdataheader,'dataType: %s\n', char(typeStr(precind)));
		fprintf('dataType: %s\n', char(typeStr(precind)));
	elseif isequal(attr,'fov')
		fprintf(fidfftdataheader,'fov: %g %g %g %g\n',myfov(1),myfov(2),myfov(3),mydimt);
		fprintf('fov: %g %g %g %g\n',myfov(1),myfov(2),myfov(3),mydimt);
	elseif isequal(attr,'extent')
		fprintf(fidfftdataheader,'extent: %g %g %g %g\n',myextent(1),myextent(2),myextent(3),mydimt-1);
		fprintf('extent: %g %g %g %g\n',myextent(1),myextent(2),myextent(3),mydimt-1);
	else
		fprintf(fidfftdataheader,'%s\n',ss);
		fprintf('%s\n',ss);
	end
end
fclose(fidfftdataheader);
fclose(fidnativeheader);

end % end flag
	%---spectrum intensities
if absspectrumflag

	%---file names---
	myspectrumabsfile=[myshortfilename,'spectrumabs.sdt'];
	myspectrumabsheader=[myshortfilename,'spectrumabs.spr'];

		%----powerabs_--
	dimspectrumabs=size(MYSPECTRUMABS);
	ndimspectrumabs=ndims(MYSPECTRUMABS);
	mydimx=dimspectrumabs(1);
	mydimy=dimspectrumabs(2);
	if ndimspectrumabs==3 
		mydimz=1;
		mydimt=dimspectrumabs(3);
	else
		mydimz=dimspectrumabs(3);
		mydimt=dimspectrumabs(4);
	end
	fprintf('\n..writing %s\n',myspectrumabsfile);
	fidabsf=fopen(myspectrumabsfile,'w',BYTEORDER);
	fwrite(fidabsf,MYSPECTRUMABS,char(precStr(precind)));
	fclose(fidabsf);
	fidabsheader=fopen(myspectrumabsheader,'w');
	fidnativeheader=fopen(nativeheader,'r');
	while ~feof(fidnativeheader)
		ss=fgetl(fidnativeheader);
		[ attr, val ] = strtok(ss,':');
		if ss(1,1:4) == 'dim:'
			fprintf(fidabsheader,'dim: %d %d %d %d\n',mydimx,mydimy,mydimz,mydimt);
			fprintf('dim: %d %d %d %d\n',mydimx,mydimy,mydimz,mydimt);
		elseif ss(1,1:8) == 'fidName:'
			fprintf(fidabsheader,'fidName: %s\n',myspectrumabsfile);
			fprintf('fidName: %s\n',myspectrumabsfile);
		elseif isequal(attr,'fov')
			fprintf(fidabsheader,'fov: %g %g %g %g\n',myfov(1),myfov(2),myfov(3),mydimt);
			fprintf('fov: %g %g %g %g\n',myfov(1),myfov(2),myfov(3),mydimt);
		elseif isequal(attr,'extent')
			fprintf(fidabsheader,'extent: %g %g %g %g\n',myextent(1),myextent(2),myextent(3),mydimt-1);
			fprintf('extent: %g %g %g %g\n',myextent(1),myextent(2),myextent(3),mydimt-1);
		else
			fprintf(fidabsheader,'%s\n',ss);
			fprintf('%s\n',ss);
		end
		%if ss(1,1:9) == 'sdtOrient'
		%	tempo_orient_string=ss;		
		%	[mybegin myend]=strtok(tempo_orient_string,':');
		%	tempo_orient_string=myend(2:end);
		%end
	end
	fclose(fidabsheader);
	fclose(fidnativeheader);
end % end flag

%---normalized spectrum intensities
if normspectrumflag

	%---file names---
	mynormspectrumfilename=[myshortfilename,'normspectrumabs.sdt'];
	mynormspectrumheader=[myshortfilename,'normspectrumabs.spr'];

		%----powerabs_--
	dimnormspectrum=size(MYSPECTRUMNORM);
	ndimnormspectrum=ndims(MYSPECTRUMNORM);
	mydimx=dimnormspectrum(1);
	mydimy=dimnormspectrum(2);
	if ndimnormspectrum==3 
		mydimz=1;
		mydimt=dimnormspectrum(3);
	else
		mydimz=dimnormspectrum(3);
		mydimt=dimnormspectrum(4);
	end
	fprintf('\n..writing %s\n',mynormspectrumfilename);
	fidabsf=fopen(mynormspectrumfilename,'w',BYTEORDER);
	fwrite(fidabsf,MYSPECTRUMNORM,char(precStr(precind)));
	fclose(fidabsf);
	fidabsheader=fopen(mynormspectrumheader,'w');
	fidnativeheader=fopen(nativeheader,'r');
	while ~feof(fidnativeheader)
		ss=fgetl(fidnativeheader);
		[ attr, val ] = strtok(ss,':');
		if ss(1,1:4) == 'dim:'
			fprintf(fidabsheader,'dim: %d %d %d %d\n',mydimx,mydimy,mydimz,mydimt);
			fprintf('dim: %d %d %d %d\n',mydimx,mydimy,mydimz,mydimt);
		elseif ss(1,1:8) == 'fidName:'
			fprintf(fidabsheader,'fidName: %s\n',mynormspectrumfilename);
			fprintf('fidName: %s\n',mynormspectrumfilename);
		elseif isequal(attr,'fov')
			fprintf(fidabsheader,'fov: %g %g %g %g\n',myfov(1),myfov(2),myfov(3),mydimt);
			fprintf('fov: %g %g %g %g\n',myfov(1),myfov(2),myfov(3),mydimt);
		elseif isequal(attr,'extent')
			fprintf(fidabsheader,'extent: %g %g %g %g\n',myextent(1),myextent(2),myextent(3),mydimt-1);
			fprintf('extent: %g %g %g %g\n',myextent(1),myextent(2),myextent(3),mydimt-1);
		else
			fprintf(fidabsheader,'%s\n',ss);
			fprintf('%s\n',ss);
		end
		%if ss(1,1:9) == 'sdtOrient'
		%	tempo_orient_string=ss;		
		%	[mybegin myend]=strtok(tempo_orient_string,':');
		%	tempo_orient_string=myend(2:end);
		%end
	end
	fclose(fidabsheader);
	fclose(fidnativeheader);
end % end flag


fprintf('\n  <autopower> done for: %s\n',myshortfilename);
