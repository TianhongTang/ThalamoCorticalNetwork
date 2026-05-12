;;; ---------------------------------------------------------------------
PRO lcm, dir, FILE=file, DISPLAY=f_display, MACHINE=f_machine, $
		VERBOSE = f_verbose

;----- !!! should be compiled for update: .r lcm; compile,'lcm' !!! ---------

;+
; JP Jan 1998
;
; --- batch execution of LCModel  ---
;     revised version, avoids use of rlogin/rsh
;
; INPUT:
;	dir: location of lcm files
;	reading LCM.BAT
;
; OUTPUT:
;	LCModel output files
;
; FLAGS:
;	VERBOSE: file info
;
; RUNNING FROM UNIX:
;	compile: WAVE>.run lcm.pro; cd,'~/wave/lib'; compile,'lcm'
;	run    : pfeuffer@crusher 88% wave -r "lcm,'string'"
;	
;-

USAGE = "usage: LCM, dir, FILE=[batch file list]"

COMMON paths, stdpath

LCMBATFILE = 'lcmodel.bat'
LCMODEL = 'lcm'		; running script
FNAMEDEF = 'lcm'	; of LCModel raw data filename
LCMEXT3  = '.CONTROL'
LCMRUNEXT = '.out'
MACHINEDEF = ['ashley','ashley','ashley','ashley','ashley','ashley','ashley','amy','amy','ashley','sealth']	; strarray of possible machines
CHECKTIME = 1		; wait [sec] for repeated checking

if NOT keyword_set(f_machine) then begin
   machines = MACHINEDEF
endif else begin
   s_mach = size(f_machine)
   if s_mach( s_mach(0)+1 ) EQ 7 then $		; type 7 EQ string
      machines = f_machine $
   else $
      machines = MACHINEDEF
endelse
if NOT KEYWORD_SET(f_verbose) then $
   f_verbose = 1
if NOT KEYWORD_SET(f_display) then $
   f_display = 0

stdpath		; update stdpath, necessary for shell execution

if n_params() GE 1 then begin
   fdirlist = dir 
   sfile = ''
endif else begin 
   if KEYWORD_SET(file) then $
      batfile = file $
   else $
      batfile = STDPATH.lcm + LCMBATFILE
   sfile = ' ('+batfile+') '
		; ---- determine number of lines
   spawn, 'cat ' + batfile + ' | wc -l', spawnresult
   fdirlist = strarr( spawnresult(0) )
		; ---- read BATFILE
   openr, unit, batfile, /get_lun
   readf, unit, fdirlist
   free_lun, unit
endelse
		; ---- check fdirlist
dirarr = strarr(n_elements(fdirlist))
idirmax = 0l	
for idir=0l,n_elements(fdirlist)-1 do begin
   fdir = strsplit( strcompress(strtrim( fdirlist(idir), 2)), ' ')
   if strmid(fdir(0), strlen(fdir(0))-1, 1) EQ '/' then $
	 fdirtmp = fdir(0) $
      else $
         fdirtmp = fdir(0) + '/'
   if idir LE 240 then begin
      if checkfile( fdirtmp + FNAMEDEF + LCMEXT3, /read) then begin
         dirarr(idirmax) = fdirtmp
         idirmax = idirmax + 1l
      endif else $
         print,'!! <' + fdir(0) + '> excluded from fdirlist !!'
   endif else begin
      dirarr(idirmax) = fdirtmp
      idirmax = idirmax + 1l
   endelse
endfor
if idirmax LT 1 then begin
   print,'! nothing to do !'
   return
endif else $
   dirarr = dirarr(0:(idirmax-1))

	; ---- determine DISPLAY and MACHINE
spawn,'echo $DISPLAY', curdisplay
curdisplay = curdisplay(0)
curdisplay = strmid(curdisplay, 0, strpos(curdisplay, ':0'))
spawn,'uname -a', curmachine
curmachine = strsplit(strcompress(curmachine(0)),' ')
curmachine = curmachine(1)

ijobmax = n_elements(machines)
jobpid = lonarr(ijobmax)		; PIDs of running jobs
jobarr = replicate( -1l, ijobmax )	; LT 0 : job free
					; GE 0 : job working on dirarr(idir)

IF f_verbose GE 1 THEN $
      print,'----->> LCModel daemon - batch execution '+sfile+'<<-----'

	; ---- work through the dirarr / fdirlist
idir = 0l
while (idir LT idirmax) do begin
		; ---- start new jobs
   freejobind = where( jobarr LT 0, c_freej)
   ic_freej=0l
;;   for ic_freej=0,c_freej-1 do begin
   while (ic_freej LT c_freej) AND (idir LT idirmax) do begin
      freejob = freejobind(ic_freej)
      runmachine = machines(freejob)
      if f_display then $
	 if curdisplay EQ '' then $
            if runmachine NE curmachine then $
               rundisplay = curmachine + ':0.0' $
            else $
               rundisplay = ':0.0' $
         else $
            rundisplay = curdisplay + ':0.0' $
      else $
   	 rundisplay = ''
      rundir = dirarr(idir)

		; ---- run LCMODEL on rundir file
      tmpstr = "rsh " + runmachine + " 'setenv DISPLAY " + rundisplay + $
	   "; cd " + rundir + "; " + LCMODEL + " " + FNAMEDEF + "' &"
      ;print,tmpstr
      spawn, tmpstr, spawnres
	 		; ---- determine PID of rsh call
      if n_elements(spawnres) EQ 1 then begin
	 tmpstr1 = strsplit( strcompress(strtrim( spawnres(0), 2)), ' ')
         if n_elements(tmpstr1) GE 2 then begin
            runpid = long(tmpstr1(1))
         endif else $
            runpid = 0
      endif else $
         runpid = 0
      if runpid EQ 0 then begin
         print,'--- !! can not determine RUNPID from <'+spawnres+'>'
      endif

      IF f_verbose NE 0 THEN $
         print, string( freejob, format="('job ',I0,' started in  <')" ) $
		+ rundir +"> @" + runmachine $
		+ string(runpid, format="('[',I0,']')")
      if (runpid NE 0) then begin
         jobpid(freejob) = runpid		; store running PID
         jobarr(freejob) = idir		; store idir in jobarr
      endif	;; else skip checking
      idir = idir + 1l
      ic_freej = ic_freej + 1l
   endwhile   ; ic_freej 

   IF f_verbose GE 2 THEN $
      print,jobarr,jobpid

		; ---- create list of running jobs 
   tmpstr = "ps -e | grep rsh"
   spawn,tmpstr, spawnres
	;;info,spawnres
	;;print,'<'+spawnres+'>'
   pids = lonarr(n_elements(spawnres))
   for ipids=0,n_elements(spawnres)-1 do begin
				; eg. '7185 pts/12   0:00 rsh' 
      tmpstr = strsplit( strcompress(strtrim( spawnres(ipids), 2)), ' ')
      pids(ipids) = long(tmpstr(0))
   endfor
	;;print,'------------------'
	;;print,pids

		; ---- check whether running jobs are ready 
   for ijob=0,ijobmax-1 do $
      if jobarr(ijob) GE 0 then begin	; ---- running job
	 curdir = dirarr(jobarr(ijob))
	 curpid = jobpid(ijob)
         pids_index = where( pids EQ curpid, count)
         if count LE 0 then begin
	    jobarr(ijob) = -1l		; release running job info
      	    IF f_verbose GE 2 THEN $
      	       print, string( ijob, format="('job ',I0,' released in <')" ) $
		+ curdir +"> @" + machines(ijob) $
		+ string(curpid, format="('[',I0,']')")
  	 endif
      endif

   if idir LT idirmax then $
      wait, CHECKTIME
endwhile

IF f_verbose GE 1 THEN $
      print,'----->> LCModel daemon finished <<-----'

end
;;; ---------------------------------------------------------------------

;vn2lcm, 'jp71216', 101, DIRAPP='P1', SLICES=1, LSHIFT=1, $
;	BASIS='94T20_20', TRAMP=1.0, VOLUME= 2.03, DEGZER=0, DEGPPM=0, $
;	LCM=0, SHOWRAW=0, verbose=2, MACHINE=1, PPMST=4.5, PPMEND=0.5


;end
