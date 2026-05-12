;;; ---------------------------------------------------------------------
PRO vn2lcm, dir, filenum, DIRAPP=dirapp, SLICES=f_slices, LSHIFT=f_ls, $
	CT = f_ct, WDWEM = f_wdwem, BASIS = basisfile, $
	SHOWRAW = f_showraw,  VERBOSE = f_verbose, NODISPLAY = f_nodisplay, $
	MACHINE=f_machine, XSIZE=xsize, $
	TRAMP = tramp, VOLUME = volume, DEGZER = degzer, DEGPPM = degppm, $
	PPMST = ppmst, PPMEND = ppmend, $
	FWHMBA = fwhmba, RFWHM = rfwhm, SHIFMN = shifmn, PPMSHF = ppmshf, $
	NEACH = neach, CHOMIT = chomit, CHKEEP = chkeep, VITRO = f_vitro, $
	DKNTMN = dkntmn, SDDEGP = sddegp, NOXXT2 = f_noxxt2

;+
; JP Jan 1998
;
; --- reads VNMR data  ---
;
; INPUT:
;	dir
;	num = filenum (optional: then using STDPATH.VNMR)
;	COMMON PROCPAR, procpar, headerfile, headerblock
;
; OUTPUT:
;	written LCModel file and run
;
; FLAGS:	default	
;	DIRAPP 	''	append 'DIRAPP' to dir
;	SLICES	0	process all slices of a 2D data set
;	CT	0	explicit input of number of scans
;	LSHIFT	0	data left shift
;	WDWEM	0	data windowing: exponential multiplication in [Hz]
;	BASIS	! must be specified !	FILBAS - LCModel basis set 
;	SHOWRAW	0	PLOTRAW data from LCModel, 
;	VERBOSE 1	verbose on output
;	NODISPLAY 0	no display on X window
;
;	MACHINE	0	use default or specify eg. ['amy','amy','sealth']
;	XSIZE	2048	num complex data points: truncate/fill data 
;
; LCModel input parameter: (see manual for use)
;		default		example
;	TRAMP	1.0
;	VOLUME	1.0
;	DEGZER	0.
;	DEGPPM	0.
;	PPMST	999.
;	PPMEND	-999.
;	FWHMBA	0.009
;	RFWHM	2.5
;	SHIFMN	[-0.2,-0.1]
;	PPMSHF	0.
;	NEACH	0.
;	CHOMIT	''		['Ser','Asp']	
;	CHKEEP	''		['NAAG','Ala']	
;	VITRO	0.
;	DKNTMN	0.
;	SDDEGP	0.		0: no PK 1st / >0: set prob. / <0: LCM default 
; special:
;	NOXXT2	0.		NE 0: NO fixing of linewidth of some metab.
;
;-

USAGE = "usage: VN2LCM, dir, filenum"

COMMON paths, stdpath
COMMON PROCPAR, procpar, headerfile, headerblock

VERSION = 'VN2LCM() - JP Jan 1998'

PLOTRAW = 'plotraw'	; commands
LCMODEL = 'lcm'		;

STDEXT = '.fid'		; VNMR  data directory
LCMEXT  = '.lcm'	; of LCM data directory
LCMEXT1  = '.RAW'
LCMEXT2  = '.PLOTIN'
LCMEXT3  = '.CONTROL'
FNAMEDEF = 'lcm'	; of LCModel raw data filename
SLICEAPP = 'S'		; slice appendix
MACHINEDEF = ['ashley','ashley','ashley','ashley','ashley','ashley','ashley','amy','amy','ashley','sealth']	; strarray of possible machines
LCMBATFILE = 'lcmodel.bat'

h = { , $			; ---- RAW parameters
	comm   : strarr(3), $	; arbitrary length
	id     : '', $		; ID string
	tramp  : 0., $		;
	volume : 0., $		;

	hzpppm : 0., $		; ---- PLOTIN parameters
	nunfil : 0l, $		;
	deltat : 0., $		;
	filraw : '', $		;
	filps  : '', $		;
	ppmst  : 999., $	;
	ppmend : -999., $	;
	degzer : 0., $		;	
	degppm : 0., $		;

	title  : '', $		; ---- CONTROL parameters
	filbas : '', $		;
	filps2 : '', $		;
	filcoo : '', $		;
	vitro  : 0, $		;
	dkntmn : 0., $		;
	sddegp : 0., $		;
	fwhmba : 0., $	;
	rfwhm  : 2.5, $	;
	shifmn : [-0.2,-0.1], $	;
	ppmshf : 0., $		;
 	neach  : 0l, $		;
 	chomit : strarr(50), $	;
 	chkeep : strarr(50), $	;
	noxxt2 : 0 $		;
}

if n_params() LT 1 then begin
   print, USAGE
   return
endif 
if NOT KEYWORD_SET(dirapp) then $
   dirapp = ''
if NOT KEYWORD_SET(f_slices) then $
   f_slices = 0
if NOT keyword_set(f_verbose) then $
   f_verbose = 1
if NOT keyword_set(f_nodisplay) then $
   f_nodisplay = 0
if NOT keyword_set(f_ls) then $
   f_ls = 0
if NOT keyword_set(xsize) then $
   xsize = 2048				; complex data points
if NOT keyword_set(f_machine) then begin
   machines = MACHINEDEF
   f_machine = 0
endif else begin
   s_mach = size(f_machine)
   if s_mach( s_mach(0)+1 ) EQ 7 then $		; type 7 EQ string
      machines = f_machine $
   else $
      machines = MACHINEDEF
endelse
if NOT keyword_set(f_ct) then $
   f_ct = 0
if NOT keyword_set(basisfile) then $
   basisfile = 'unknown'
if NOT keyword_set(fwhmba) then $
   fwhmba = 0.009		; on 9.4: (1.5+2) hz/400
if NOT keyword_set(rfwhm) then $
   rfwhm = 2.5			; 
if NOT keyword_set(f_vitro) then $
   f_vitro = 0		
if NOT keyword_set(dkntmn) then $
   dkntmn = 0		
if NOT keyword_set(sddegp) then $
   sddegp = 0.		
if NOT keyword_set(shifmn) then $
   shifmn = [-0.2,-0.1]
if n_elements(shifmn) NE 2 then begin
   print,'! n_elements(SHIFMN) ne 2 -> default setting used !'
   shifmn = [-0.2,-0.1]
endif
if NOT keyword_set(ppmshf) then $
   ppmshf = 0
if NOT keyword_set(neach) then $
   neach = 0
if NOT keyword_set(chomit) then $
   chomit = ''
if NOT keyword_set(chkeep) then $
   chkeep = ''
if NOT keyword_set(f_noxxt2) then $
   f_noxxt2 = 0		
IF(strpos(dir, STDEXT) GT -1) THEN $
   fdir = strmid(dir, 0, strpos(dir, STDEXT)) $
ELSE $
   fdir = dir

	; ---- determine vndir
IF (n_params() EQ 2) THEN BEGIN
   if dirapp NE '' then $
      vndir = STDPATH.vnp + fdir + "_" + mystring(filenum) $
		+ dirapp + STDEXT $
   else $
      vndir = STDPATH.vnmr + fdir + "_" + mystring(filenum) $
		+ STDEXT
ENDIF ELSE $
   vndir = fdir + dirapp + STDEXT

	; ---- determine lcmdir
	;      whole path: +SLICEAPP+islice+LCMEXT+'/'
IF (n_params() EQ 2) THEN begin
   j = checkfile(STDPATH.lcm + fdir, /write, IS_DIR=is_dir)
   if is_dir NE 1 then begin		; create directory
      spawn,'mkdir '+STDPATH.lcm + fdir, result
      if result(0) NE '' then print,'MKDIR :'+result
   endif
   j = checkfile(STDPATH.lcm + fdir, /write, IS_DIR=is_dir)
   if is_dir EQ 1 then $
      dirstr = "/" $
   else $
      dirstr = "_"
   lcmdir = STDPATH.lcm + fdir + dirstr + mystring(filenum) + dirapp 
endif ELSE $
   lcmdir = fdir + dirapp


IF (n_params() EQ 2) THEN $
   data = readvnmr( dir, filenum, DIRAPP=dirapp, ct=ctcount ) $
else $
   data = readvnmr( dir, DIRAPP=dirapp, ct=ctcount )

if KEYWORD_SET(f_wdwem) then $
   wdwem = double(f_wdwem)/procpar.sw $
else $
   wdwem = 0

;;;dataorg = data
;;;bc, dataorg, LSHIFT = f_ls
bc, data, LSHIFT = f_ls, WDWEM=wdwem, /BCFT	; conversion to float

		; ---- calibrate according to ctcount
if f_ct NE 0 then begin
   data = data / float(f_ct) 
   print, string(f_ct, format='("--- CT = ",G0.0,"  Calibration")')
endif else begin
   j = where(ctcount LT 1, c_ct)
   if c_ct EQ 0 then begin
      for ict=0,n_elements(ctcount)-1 do $
          data(*,ict) = data(*,ict) / float(ctcount(ict))
      print,  "--- CT = " + strjoin(string(ctcount,format='(G0.0)'))  $
		+ "  Calibration"
   endif else begin
      print, 'CTCOUNT is zero !'
      print, ctcount
      print, '!!! use CT=ctcount to calibrate !'
      stop
   endelse
endelse

		; ---- truncate/fill according XSIZE
s_td = size(data)
datalenin = s_td(1)
datalen = long(xsize)*2
datain = reform( data, s_td(1), s_td(s_td(0)+2)/ s_td(1))
data = fltarr( datalen, s_td(s_td(0)+2)/ s_td(1) )
if datalenin LE datalen then $
   data(0:(datalenin-1),*) = datain $
else $
   data = datain(0:(datalen-1),*)
if datalenin NE datalen then $
   print,string(xsize,format='("-> data truncated/filled up to XSIZE = ",I0)')

				; ---- RAW parameters
h.comm(0) = 'written by ' + VERSION	    ; check h.comm: strarr(?) !!
h.id = vndir
if KEYWORD_SET(tramp) then $
   h.tramp = tramp $
else $
   h.tramp = 1.0
if KEYWORD_SET(volume) then $
   h.volume = volume $
else $
   h.volume = 1.0
				; ---- PLOTIN parameters
h.hzpppm = procpar.sfrq
h.nunfil = datalen/2
h.deltat = procpar.at/(procpar.np/2)
h.filraw = FNAMEDEF + '.RAW'
h.filps  = FNAMEDEF + '_RAW.PS'
if KEYWORD_SET(ppmst) then $
   h.ppmst = ppmst 
if KEYWORD_SET(ppmend) then $
   h.ppmend = ppmend 
if KEYWORD_SET(degzer) then $
   h.degzer = degzer 
if KEYWORD_SET(degppm) then $
   h.degppm = degppm 
				; ---- CONTROL parameters
str_title = ' '
for i=0,n_elements(procpar.text)-1 do begin
   if procpar.text(i) NE '' then $
      str_title = str_title + procpar.text(i) + ' '
endfor
h.title   = strmid( str_title, 0, 120)		; max length is 120 chars
spawn,'echo ' + STDPATH.lcm, result
h.filbas  = result(0) + basisfile + '.BASIS'
h.filps2  = FNAMEDEF + '.PS'
h.filcoo  = FNAMEDEF + '.COORD'
h.vitro   = f_vitro
h.dkntmn  = dkntmn
h.sddegp  = float(sddegp)
h.fwhmba  = fwhmba
h.rfwhm   = rfwhm
h.shifmn  = shifmn
h.ppmshf  = ppmshf
h.neach   = neach
if chomit(0) NE '' then $
   for iomit=0,n_elements(chomit)-1 do $
      h.chomit(iomit)  = chomit(iomit)
if chkeep(0) NE '' then $
   for ikeep=0,n_elements(chkeep)-1 do $
      h.chkeep(ikeep)  = chkeep(ikeep)
h.noxxt2  = f_noxxt2

IF (n_params() EQ 2) THEN $
   writelcmraw, data, dir, filenum, DIRAPP=dirapp, SLICES=f_slices, $
	HEADER=h, VERBOSE = f_verbose $
else $
   writelcmraw, data, dir, DIRAPP=dirapp, SLICES=f_slices, $
	HEADER=h, VERBOSE = f_verbose

if f_slices(0) NE 0 then $
   strsl = SLICEAPP + '*' $
else $
   strsl = ''
j = findfile(lcmdir + strsl + LCMEXT + '/' + FNAMEDEF + LCMEXT1, count=nslices)
if (nslices LE 0) then begin
   print,'Not found <' + lcmdir+strsl+LCMEXT+'/'+FNAMEDEF+LCMEXT1 + '> !'
   return
endif 
fdirlist = strarr(nslices)

	; ---- determine DISPLAY and MACHINE
spawn,'echo $DISPLAY', curdisplay
curdisplay = curdisplay(0)
curdisplay = strmid(curdisplay, 0, strpos(curdisplay, ':0'))
spawn,'uname -a', curmachine
curmachine = strsplit(strcompress(curmachine(0)),' ')
curmachine = curmachine(1)


for islice=0,nslices-1 do begin

   if f_slices(0) NE 0 then $
      strsl = SLICEAPP + string(islice, format='(I0)') $
   else $
      strsl = ''
   fdirslice = lcmdir + strsl + LCMEXT + '/' 
   fdirlist(islice) = fdirslice

   runmachine = machines( long(abs(islice)) MOD n_elements(machines) )
   if curdisplay EQ '' then $
      if runmachine NE curmachine then $
         rundisplay = curmachine + ':0.0' $
      else $
         rundisplay = ':0.0' $
   else $
      rundisplay = curdisplay + ':0.0'

   if KEYWORD_SET(f_showraw) then begin
      IF f_verbose NE 0 THEN $
         print, 'writing <' + fdirslice + h.filps + '> (' + runmachine + ')'

      tmpstr = "rsh " + runmachine + " 'setenv DISPLAY " + rundisplay + $
   	   "; cd " + fdirslice + "; " + PLOTRAW + " " + FNAMEDEF + "' &"
      spawn, tmpstr, result
   endif

   if KEYWORD_SET(f_lcm) then begin
      IF f_verbose NE 0 THEN BEGIN
         print, 'LCModel started on   ' + runmachine + '   <' + $
	   fdirslice + FNAMEDEF + '.CONTROL>'
         print, '-> output is <' + h.filps2 + '><' + h.filcoo + '>' 
      ENDIF
      tmpstr = "rsh " + runmachine + " 'setenv DISPLAY " + rundisplay + $
	   "; cd " + fdirslice + "; " + LCMODEL + " " + FNAMEDEF + "' &"
      spawn, tmpstr, result
      print,result
   endif

endfor 	; --- islice

		; ---- write lcmodel batch file
ch = ''
read, '? overwrite batchfile <'+STDPATH.lcm+LCMBATFILE+'> ?  (yes/no)[no] ',ch
if strmid(strupcase(strtrim(ch, 1)), 0, 3) EQ 'YES' then $
   openw, unit, STDPATH.lcm + LCMBATFILE, /get_lun $
else $
   openw, unit, STDPATH.lcm + LCMBATFILE, /get_lun, /APPEND
for islice=0,nslices-1 do $
   printf, unit, fdirlist(islice)
free_lun, unit

	; ----- running of LCModel ----
ch = ''
read, '? start LCModel ?  (yes/no)[yes] ',ch
if strmid(strupcase(strtrim(ch, 1)), 0, 3) NE 'NO' then $
   if f_nodisplay NE 0 then $ 
      lcm, fdirlist, DISPLAY=0, MACHINE=f_machine $
   else $
      lcm, fdirlist, DISPLAY=1, MACHINE=f_machine

end
;;; ---------------------------------------------------------------------

;vn2lcm, '~/data/vnmr/jp71104/11' , '~/LCModel/jptest', /SHOWRAW
;vn2lcm, '~/data/vnmr/jk/jk04287.naa.s1', '~/LCModel/4T_NAA', $
;	TRAMP=2, VOLUME= 0.5, DEGZER=-63, DEGPPM=0, PPMST=4.5, PPMEND=0.1, $
;	/SHOWRAW, verbose=2

;vn2lcm, '~/data/vnmr/jk/kk_090296/jk09026.s1' , '~/LCModel/jk09026.s1', $
;	TRAMP=1.0, VOLUME= 1.0, DEGZER=-30, DEGPPM=0, PPMST=4.5, PPMEND=0.1, $
;	LCM=1, /SHOWRAW, verbose=2

;end
