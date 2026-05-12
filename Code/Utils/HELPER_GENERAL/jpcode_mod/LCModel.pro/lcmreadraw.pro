PRO lcmreadraw, filename, data, header, VERBOSE=f_verbose

COMMON PATHS, STDPATH

STDEXT = '.RAW'
EOH = '$END'
MAXLINES = 80		; max number of text lines in header
MAXDATA = 8*1024l	; max number of data lines to read

IF NOT KEYWORD_SET(f_verbose) THEN $
   f_verbose = 1

file = filename
IF (strpos(file, STDEXT) EQ -1) THEN $
    file = file + STDEXT

IF f_verbose THEN $
   print, 'reading <'+file+'>'
openr, unit, file, /get_lun

iline = 0
line = ''
header = strarr(MAXLINES)
f_eoh = 0		; flag: END OF HEADER
WHILE (NOT eof(unit)) AND (f_eoh EQ 0) DO BEGIN
   readf, unit, line
   if (iline LT MAXLINES) then begin
      header(iline) = line
      iheader = iline
   endif else if (iline EQ MAXLINES) then $
      print, '!!! <'+file+'> header: MAXLINES exeeded !!!'

   iline = iline + 1
   if (strpos(line, EOH) GE 0) then $
      f_eoh = 1
ENDWHILE
header = header(0:iheader)


iline = 0
WHILE (NOT eof(unit)) DO BEGIN
 			; --- determine numcols
   if iline EQ 0 then begin
      point_lun, -unit, lunpos
      readf, unit, line
      point_lun, unit, lunpos
      numcols = n_elements( strsplit(strcompress(strtrim(line,2)),' ') )
      rdat = fltarr(numcols)
      data = fltarr(MAXDATA*numcols)
   endif 
   readf, unit, rdat
   if iline LT MAXDATA then begin
      data( iline*numcols:(iline*numcols+numcols-1) ) = rdat
      idata = iline
   endif else if (iline EQ MAXDATA) then $
      print, '!!! <'+file+'> data: MAXDATA exeeded !!!'

   iline = iline + 1
ENDWHILE
data = data( 0:(idata*numcols+numcols-1) )

free_lun, unit

IF f_verbose THEN BEGIN
   print, '--- ndata = ' + $
	string(n_elements(data), format = '(I0)' )
   for i=0,n_elements(header)-1 do $   
      print, header(i)
   !p.multi = [0,2,2]
      plot, data
      fdata = ft(data)
      plot, float(fdata)
      plot, abs(fdata)
      plot, imaginary(fdata)
  !p.multi = 0
ENDIF

END
;;; ---------------------------------------------------------------------

;filename = '~/LCModel/Glu_3'
;filename = '~/LCModel/Tau_3'
;filename = '~/LCModel/Test'
;lcmreadraw, filename, data, header
;info, data, header

;end