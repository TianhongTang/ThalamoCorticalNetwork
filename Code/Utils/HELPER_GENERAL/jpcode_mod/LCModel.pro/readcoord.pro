FUNCTION readcoord, dir, filenum, DIRAPP=dirapp, SLICES=f_slices, $
	VERBOSE=f_verbose

;+
; JP Jan 1998
;
; --- reads data in COORD format (output from LCModel)  ---
;
; INPUT:
;	dir
;	filenum 	(optional: then using STDPATH.LCM)
;
; OUTPUT:
;      lcmcoord = { , $
;	metab : tcmetabolite
;	conc  : tcconc
; 	relconc:tcrelconc
;	SD    :	tcSD
;	data  : data
;	mdata : mdata
;	header: header1
;	tconc : header2
;	tmisc : header3
;	tdiag : header4
;	tinput: header5
;      }
;   endif 
;
; FLAGS:
;	DIRAPP : append 'DIRAPP' to dir 
;	VERBOSE: file info
;
; ORDER OF reading:
; 	header1: general info, date
;	header2: concentration table
;	header3: misc. output table
;	data: dblarray(NY, NDATABLOCKS+nconc-NCONCRED)
;		0: ppm axis
;		1: phased data points
;		2: fit to the follow
;		3: background values
;	       >3: metabolite data
;	header4: diagnostic table
;	header5: table of input changes
;		
;-

USAGE = "usage: READCOORD, data, dir, filenum"
		   
COMMON paths, stdpath

if n_params() LT 1 then begin
   print, USAGE
   stop
endif
if NOT KEYWORD_SET(dirapp) then $
   dirapp = ''
if NOT KEYWORD_SET(f_slices) then $
   f_slices = 0 $
else if n_elements(f_slices) LT 1 then $
   f_slices = [fslices(0),0]
if NOT KEYWORD_SET(f_verbose) then $
    f_verbose = 0

STDEXT = '.fid'		; of data directory
LCMEXT = '.lcm'
LCMEXT4 = '.COORD'
FNAMEDEF = 'lcm'	; of LCModel raw data filename
SLICEAPP = 'S'		; slice appendix

NUMHLINES = 4		; initial header lines
NDATABLOCKS = 4		; ppm-axis, phased data, fit, background 
NELTDIAG = 10		; minimum number of lines in tdiag

;;if f_slices(0) NE 0 then $
;;   strsl = SLICEAPP + string(f_slices(1), format='(I0)') $
;;else $
;;   strsl = ''

IF(strpos(dir, STDEXT) GT -1) THEN $
   fdir = strmid(dir, 0, strpos(dir, STDEXT)) $
ELSE IF(strpos(dir, LCMEXT) GT -1) THEN $
   fdir = strmid(dir, 0, strpos(dir, LCMEXT)) $
ELSE IF(strpos(dir, LCMEXT4) GT -1) THEN $
   fdir = strmid(dir, 0, strpos(dir, LCMEXT4)) $
ELSE $
   fdir = dir

		; --- whole path: +SLICEAPP+islice+LCMEXT+'/'
IF (n_params() EQ 2) THEN begin
   j = checkfile(STDPATH.lcm + fdir, /write, IS_DIR=is_dir)
   if is_dir EQ 1 then $
      dirstr = "/" $
   else $
      dirstr = "_"
   fdir = STDPATH.lcm + fdir + dirstr + mystring(filenum) + dirapp 
endif ELSE $
   fdir = fdir + dirapp 
lcmfile = FNAMEDEF + LCMEXT4

if f_slices(0) NE 0 then $
   strsl = SLICEAPP + '*' $
else $
   strsl = ''
j = findfile(fdir + strsl + LCMEXT + '/' + lcmfile, count=nslices)
if (nslices LE 0) then begin
   print,'Not found <' + fdir + strsl + LCMEXT + '/' + lcmfile + '> !'
   return, 0
endif 

if f_slices(0) EQ 0 then $
   nslices = 1
	
line = ''
ymax = 0
for islice=0,nslices-1 do begin		; ---- 2D data set

   if f_slices(0) NE 0 then $
      strsl = SLICEAPP + string(islice, format='(I0)') $
   else $
      strsl = ''
   fdirslice = fdir + strsl + LCMEXT + '/'
   file = fdirslice + lcmfile


   print, 'reading <'+file+'>'
   openr, unit, file, /get_lun

			; --- read header1 ---------------------------
   header1 = strarr( NUMHLINES )
   readf, unit, header1
   IF f_verbose GE 1 THEN $
      for i=0,n_elements(header1)-1 do $   
         print, header1(i)

			; --- read header2, determine NCONC ----------
   readf, unit, line
   IF f_verbose GE 1 THEN print, line
   nconc = long(line) - 1

   tcmetabolite = strarr(nconc)
   tcconc = fltarr(nconc)
   tcrelconc = fltarr(nconc)
   tcSD = lonarr(nconc)

   header2 = strarr(nconc+1)
   readf, unit, header2
   for i=0,n_elements(header2)-1 do begin 
      tmpstr =  strsplit(strcompress(strtrim(header2(i), 1)),' ')
      if i GE 1 then begin
   	 tcmetabolite(i-1) = tmpstr(3)
   	 tcconc(i-1) = float(tmpstr(0))
   	 tcrelconc(i-1) = float(tmpstr(2))
   	 tcSD(i-1) = long(tmpstr(1))	 
      endif
      IF f_verbose GE 1 THEN $
         print, header2(i)
   endfor

			; --- read header3 ---------------------------
   readf, unit, line
   IF f_verbose GE 1 THEN print, line
   header3 = strarr( long(line))
   readf, unit, header3
   IF f_verbose GE 1 THEN $
      for i=0,n_elements(header3)-1 do $   
         print, header3(i)

			; --- read data blocks ----------------------
   FOR idata=0,NDATABLOCKS-1 DO BEGIN

 			; --- determine numNY
      readf, unit, line
      IF f_verbose GE 2 THEN print, line
      IF idata EQ 0 THEN BEGIN	; --- read NY from line/ init dataarr
         numNY = long(line)	
	 data = dblarr(numNY, NDATABLOCKS)
      ENDIF
 			; --- determine numcols, numlines, restelements
      point_lun, -unit, lunpos 	
      readf, unit, line
      point_lun, unit, lunpos
      numcols = n_elements( strsplit(strcompress(strtrim(line,2)),' ') )
      numlines = long( numNY / numcols)
      restelements = numNY MOD numcols
      bdata = dblarr(numcols, numlines)
      readf, unit, bdata
      if restelements GT 0 then begin
         rdata = dblarr(restelements)
	 readf, unit, rdata
	 data(*,idata) = [reform(bdata, numcols*numlines), rdata]
      endif else $
	 data(*,idata) = reform(bdata, numcols*numlines)
	 
   ENDFOR   ; idata

			; --- read mdata blocks -----------------
   readf, unit, line
   mdata = dblarr(numNY, nconc)

   WHILE strpos(strlowcase(line), 'diagnostic') LT 0 DO BEGIN

      IF f_verbose GE 2 THEN print, line
      tmpstr = strsplit(strcompress(strtrim(line,2)),' ')
      metabolite = tmpstr(0)

 	; --- numcols/numlines/restelements, bdata/rdata as ABOVE
      readf, unit, bdata
      if restelements GT 0 then begin
	 readf, unit, rdata
	 tmpdata = [reform(bdata, numcols*numlines), rdata]
      endif else $
	 tmpdata = reform(bdata, numcols*numlines)
	 
      mind = where(tcmetabolite EQ metabolite, count)
      if count GT 0 then $
         mdata(*, mind(0)) = tmpdata

      readf, unit, line

   ENDWHILE
			; --- read header4 ---------------------------
   IF f_verbose GE 1 THEN print, line
   header4 = strarr( max([1l, long(line)]) ) 	; minimum 1 line
   readf, unit, header4
   IF f_verbose GE 1 THEN $
      for i=0,n_elements(header4)-1 do $   
         print, header4(i)

			; --- read header5 ---------------------------
   readf, unit, line
   IF f_verbose GE 1 THEN print, line
   header5 = strarr( long(line) )
   readf, unit, header5
   IF f_verbose GE 1 THEN $
      for i=0,n_elements(header5)-1 do $   
         print, header5(i)

   free_lun, unit

		; ---- definition of LCMCOORD struct

   if islice EQ 0 then begin
      nel_tdiag = max( [n_elements(header4), NELTDIAG] )
      lcmcoord = { , $
	metab:  ( strarr(nconc, nslices)), $
	conc: 	( fltarr(nconc, nslices)), $
 	relconc:( fltarr(nconc, nslices)), $
	SD:	( lonarr(nconc, nslices)), $
	data  : ( dblarr(numNY, NDATABLOCKS, nslices)), $
	mdata : ( dblarr(numNY, nconc, nslices)), $
	header: ( strarr(n_elements(header1), nslices)), $
	tconc : ( strarr(n_elements(header2), nslices)), $
	tmisc : ( strarr(n_elements(header3), nslices)), $
	tdiag : ( strarr(nel_tdiag, nslices)), $
	tinput: ( strarr(n_elements(header5), nslices)) $
      }
   endif 

   lcmcoord.metab(*,islice)   = tcmetabolite
   lcmcoord.conc(*,islice)    = tcconc
   lcmcoord.relconc(*,islice) = tcrelconc
   lcmcoord.SD(*,islice)      = tcSD
   lcmcoord.data(*,*,islice)  = data
   lcmcoord.mdata(*,*,islice) = mdata
   lcmcoord.header(*,islice)  = header1
   lcmcoord.tconc(*,islice)   = header2
   lcmcoord.tmisc(*,islice)   = header3
;;;lcmcoord.tdiag(*,islice)   = header4
   lcmcoord.tinput(*,islice)  = header5

		; ---- tdiag may vary in number of lines
   for itdiag=0,n_elements(header4)-1 do $
       lcmcoord.tdiag(itdiag,islice) = header4(itdiag)


   IF f_verbose GE 1 THEN BEGIN
      ymax = max([max(data(*,1), min=ymin),ymax])
      !p.multi = [0,1,2]
      plot, data(*,0), data(*,1), xstyle=1, xminor=10, ystyle=1, $
	 xrange=[data(0,0),data(numNY-1,0)], $
	 yrange=[ymin,ymax], psym=3, $
	 xtitle = 'ppm', ytitle = 'Spectrum / a.u.'
      oplot, data(*,0), data(*,2)
      oplot, data(*,0), data(*,3), psym=3 ;;, linestyle=1
      plot, data(*,0), data(*,1) - data(*,2), xstyle=1, xminor=10, $
 	 xrange=[data(0,0),data(numNY-1,0)], $
	 xtitle = 'ppm', ytitle = 'Residuum'
      !p.multi = 0
   ENDIF

   IF f_verbose GE 2 AND nslices GT 1 THEN $
      read, '-> RET to continue ', line

endfor 	; --- islice

IF f_verbose GE 2 THEN $
   info,/struct,lcmcoord

return, lcmcoord

END
;;; ---------------------------------------------------------------------

;result = readcoord('~/LCModel/96.version/Test.COORD', /verbose)

;j = readcoord('/home/raid20/pfeuffer/data/lcm/jp71216_101P1S0.lcm', verbose=2)
;j = readcoord('jp71216', 101, DIRAPP='P1', /SLICES, verbose=1)
; j = readcoord('jp80115', 101, DIRAPP='P1', /SLICES, verbose=0)

;end