PRO writelcmraw, data, dir, filenum, SLICES=f_slices, HEADER=header, $
	DIRAPP=dirapp, VERBOSE=f_verbose

;+
; JP Jan 1998
;
; --- write files in LCModel format  ---
;     called by VN2LCM.PRO
;
; INPUT:
;	data
;	dir
;	num = filenum (optional: then using STDPATH.VNMR)
; 	header
;
; OUTPUT:
;	LCModel files are written
;
; FLAGS:
;	DIRAPP : append 'DIRAPP' to dir 
;	VERBOSE: file info
;	
;-

USAGE = "usage: WRITELCMRAW, data, dir, filenum, header"

COMMON paths, stdpath

STDEXT = '.fid'		; of VN data directory
LCMEXT  = '.lcm'	; of LCM data directory
FNAMEDEF = 'lcm'	; of LCModel raw data filename
SLICEAPP = 'S'		; slice appendix

			; definitions for LCM files
LCMEXT1  = '.RAW'
NLBEGIN = '$NMID'
NLEND   = '$END'
FMTDAT  = '(8E13.5)'		; to define FORTRAN format for RAW data
NCOLS   = 8

LCMEXT2  = '.PLOTIN'
NLBEGIN2 = '$PLTRAW'
NLEND2   = '$END'

LCMEXT3  = '.CONTROL'
NLBEGIN3 = '$LCMODL'
NLEND3   = '$END'


if NOT KEYWORD_SET(header) then begin
   print,'HEADER must be supplied !'
   stop
endif
if NOT KEYWORD_SET(dirapp) then $
   dirapp = ''
if NOT KEYWORD_SET(f_slices) then $
   f_slices = 0
IF NOT KEYWORD_SET(f_verbose) THEN $
   f_verbose = 0

IF(strpos(dir, STDEXT) GT -1) THEN $
   fdir = strmid(dir, 0, strpos(dir, STDEXT)) $
ELSE IF(strpos(dir, LCMEXT1) GT -1) THEN $
   fdir = strmid(dir, 0, strpos(dir, LCMEXT1)) $
ELSE IF(strpos(dir, LCMEXT2) GT -1) THEN $
   fdir = strmid(dir, 0, strpos(dir, LCMEXT2)) $
ELSE IF(strpos(dir, LCMEXT3) GT -1) THEN $
   fdir = strmid(dir, 0, strpos(dir, LCMEXT3)) $
ELSE IF(strpos(dir, LCMEXT) GT -1) THEN $
   fdir = strmid(dir, 0, strpos(dir, LCMEXT)) $
ELSE $
   fdir = dir

		; --- whole path: +SLICEAPP+islice+LCMEXT+'/'
IF (n_params() EQ 3) THEN begin
   j = checkfile(STDPATH.lcm + fdir, /write, IS_DIR=is_dir)
   if is_dir EQ 1 then $
      dirstr = "/" $
   else $
      dirstr = "_"
   fdir = STDPATH.lcm + fdir + dirstr + mystring(filenum) + dirapp 
endif ELSE $
   fdir = fdir + dirapp 

if f_slices(0) NE 0 then $
   strsl = SLICEAPP + '0' $
else $
   strsl = ''
j = findfile(fdir + strsl + LCMEXT + '/' + FNAMEDEF + '*', count=c_fdir)
if (c_fdir GT 0) then begin
   ch = ''
   read, '!! overwrite <'+fdir+strsl+LCMEXT+'> ?  (yes/no)[no] ', ch
   if strmid(strupcase(strtrim(ch, 1)), 0, 3) NE 'YES' then begin
      return
   endif
endif 

IF f_verbose GT 1 THEN $
   info, /struct, header

s_data = size(data)
if s_data(0) GE 2 then begin		
   nslices = s_data( s_data(0)+2 ) / s_data(1)
endif else $
   nslices = 1

ldata = reform(data, s_data(1), nslices )	; --- lcm data
if f_slices(0) EQ 0 then $
   nslices = 1
strtitle = header.title
	
for islice=0,nslices-1 do begin		; ---- 2D data set

   if f_slices(0) NE 0 then $
      strsl = SLICEAPP + string(islice, format='(I0)') $
   else $
      strsl = ''
   fdirslice = fdir + strsl + LCMEXT + '/'

   j = findfile(fdirslice + FNAMEDEF + '*', count=c_fdir)
   if (c_fdir LE 0) then begin
      spawn, 'mkdir ' + fdirslice, result
      if result NE '' then $
         print, result
   endif
		; ----  ; max length is 120 chars
   header.title = strmid( fdirslice + strtitle, 0, 120)  
   isdata = ldata(*, islice)

			; ---- writing RAW file
file = fdirslice + FNAMEDEF + LCMEXT1
IF f_verbose NE 0 THEN $
   print, 'writing <'+file+'>'
openw, unit, file, /get_lun

for i = 0, n_elements(header.comm)-1 do $
   printf, unit, ' ' + header.comm(i)
printf, unit, ' ' + NLBEGIN
str_id = strcompress(strupcase(header.id),/remove_all)
if strlen(str_id) GT 16 then $
   str_id = strmid(str_id, strlen(str_id)-16, 16)
printf, unit, " ID='" + str_id + "'"
printf, unit, ' FMTDAT=' + "'" + FMTDAT + "'"
printf, unit, ' TRAMP=' + string(header.tramp,format='(E13.5)')
printf, unit, ' VOLUME=' + string(header.volume,format='(E13.5)')
printf, unit, ' ' + NLEND

nrows = long( n_elements(isdata) / NCOLS )
for irow= 0,nrows-1 do $
   printf, unit, string( isdata( (irow*NCOLS):((irow+1)*NCOLS-1)) , $
	format=FMTDAT )

free_lun, unit

ndatalost = n_elements(isdata) MOD NCOLS
if ndatalost NE 0 then $
   print, '--- writelcmraw: ' + string(ndatalost) + ' data points lost !'


			; ---- writing PLOTIN file
file = fdirslice + FNAMEDEF + LCMEXT2
IF f_verbose NE 0 THEN $
   print, 'writing <'+file+'>'
openw, unit, file, /get_lun

printf, unit, ' ' + NLBEGIN2
printf, unit, ' HZPPPM=' + string(header.hzpppm,format='(G0.0)')
printf, unit, ' NUNFIL=' + string(header.nunfil,format='(I0)')
printf, unit, ' DELTAT=' + string(header.deltat,format='(G0.0)')
printf, unit, " FILRAW='" + header.filraw + "'"
printf, unit, " FILPS='" + header.filps + "'"
printf, unit, ' PPMST=' + string(header.ppmst,format='(G0.0)')
printf, unit, ' PPMEND=' + string(header.ppmend,format='(G0.0)')
printf, unit, ' DEGPPM=' + string(header.degppm,format='(G0.0)')
printf, unit, ' DEGZER=' + string(header.degzer,format='(G0.0)')
printf, unit, ' ' + NLEND2

free_lun, unit

			; ---- writing CONTROL file
file = fdirslice + FNAMEDEF + LCMEXT3
IF f_verbose NE 0 THEN $
   print, 'writing <'+file+'>'
openw, unit, file, /get_lun

if strpos(header.filbas, '94TC13') GE 0 then begin

	; -------------------------- 94TC13 basis handling
printf, unit, ' ' + NLBEGIN3
printf, unit, " TITLE='" + header.title + "'"
printf, unit, " OWNER='Center for MR Research, University of Minnesota Medical School'"
printf, unit, " KEY=   940737318, 818517978, 595502896, 274233856, 731795398, 83917292, "		;;; old keys
printf, unit, " KEY(7) =856120280, 674658387,  5741772, 633041300, 190794767,446867949 "		;;; old keys, ..., amy
printf, unit, " KEY(13)=814679863, 264741136, 792921115,856690523, 519810613,  88318785, 134386697"	;;;lefty,ely,ashley(?),poplar,mink,mudro,cortex

printf, unit, " PGNORM='US'"
printf, unit, " IETCOU=3"
printf, unit, " NREFPK=1"
printf, unit, " PPMREF(1,2)=2.35"	; Glu C4 reference
printf, unit, " NUSE1=6"
printf, unit, " CHUSE1='Glu4','Gln4','Glu3','Gln3','Glu2','Gln2','Lac','NAA','Cr'"

printf, unit, " FILPS='" + header.filps2 + "'"
printf, unit, " FILCOO='" + header.filcoo + "'"
printf, unit, " CHCOMB(9)='Gln3+Glu3'"	
printf, unit, " CHCOMB(8)='Gln2+Glu2'"	
printf, unit, " CHCOMB(3)='Cr+PCr'"	
printf, unit, " CHCOMB(2)='GPC+PCho'"
printf, unit, ' NCOMBI=9'		
printf, unit, " NAMREL='Glu4'"	
printf, unit, ' LPRINT=6'
printf, unit, " FILPRI='lcm.PRI'"
printf, unit, ' LCOORD=9'
ind = where(header.chkeep NE '', count)
if count GT 0 then begin
   chkeep = header.chkeep(ind)
   printf, unit, ' NKEEP=' + string(n_elements(chkeep),format='(I0)')
   for ikeep=0,n_elements(chkeep)-1 do $
      printf, unit, string(ikeep+1, format='(" CHKEEP(",I0,")=")') $
		+ "'" + chkeep(ikeep) + "'"
endif
ind = where(header.chomit NE '', count)
if count GT 0 then begin
   chomit = header.chomit(ind)
   printf, unit, ' NOMIT=' + string(n_elements(chomit),format='(I0)')
   for iomit=0,n_elements(chomit)-1 do $
      printf, unit, string(iomit+1, format='(" CHOMIT(",I0,")=")') $
		+ "'" + chomit(iomit) + "'"
endif
printf, unit, " NAMEAC(5)='GABA'"
printf, unit, " NAMEAC(4)='GSH'"
printf, unit, " NAMEAC(3)='Ins'"
printf, unit, " NAMEAC(2)='Glc'"
printf, unit, " NAMEAC(1)='Mac'"
printf, unit, ' NEACH=' + string(header.neach,format='(I0)')
if header.sddegp GE 0 then $
   printf, unit, ' SDDEGP=' + string(header.sddegp,format='(G0.0)')
if header.degzer NE 0 then begin
   printf, unit, ' DEGZER=' + string(header.degzer,format='(G0.0)')
   printf, unit, ' SDDEGZ=1'
endif
if header.ppmshf NE 0 then $
   printf, unit, ' PPMSHF=' + string(header.ppmshf,format='(G0.0)')
printf, unit, ' FWHMBA=' + string(header.fwhmba,format='(G0.0)')
printf, unit, ' RFWHM=' + string(header.rfwhm,format='(G0.0)')
if header.dkntmn NE 0 then $
   printf, unit, ' DKNTMN=' + string(header.dkntmn,format='(G0.0)') $
else begin
   printf, unit, " ECCDON=.TRUE."
   if header.vitro NE 0 then $
      printf, unit, ' VITRO=.TRUE.'
endelse
printf, unit, ' PPMEND=' + string(header.ppmend,format='(G0.0)')
printf, unit, ' PPMST=' + string(header.ppmst,format='(G0.0)')
printf, unit, " FILRAW='" + header.filraw + "'"
printf, unit, " FILBAS='" + header.filbas + "'"
printf, unit, ' DELTAT=' + string(header.deltat,format='(G0.0)')
printf, unit, ' NUNFIL=' + string(header.nunfil,format='(I0)')
printf, unit, ' HZPPPM=' + string(header.hzpppm,format='(G0.0)')
printf, unit, ' ' + NLEND3

	; -------------------- end of 94TC13 basis handling
endif else begin

	; -------------------------- common basis handling
printf, unit, ' ' + NLBEGIN3
printf, unit, " TITLE='" + header.title + "'"
printf, unit, " OWNER='Center for MR Research, University of Minnesota Medical School'"
printf, unit, " KEY=   940737318, 818517978, 595502896, 274233856, 731795398, 83917292, "		;;; old keys
printf, unit, " KEY(7) =856120280, 674658387,  5741772, 633041300, 190794767,446867949 "		;;; old keys, ..., amy
printf, unit, " KEY(13)=814679863, 264741136, 792921115,856690523, 519810613,  88318785, 134386697"	;;;lefty,ely,ashley(?),poplar,mink,mudro,cortex

printf, unit, " PGNORM='US'"
printf, unit, " FILPS='" + header.filps2 + "'"
printf, unit, " FILCOO='" + header.filcoo + "'"
if header.noxxt2 EQ 0 then begin
   printf, unit, " ALSDT2(5)=0.01"
   printf, unit, " ALSDT2(4)=0.01"
   printf, unit, " ALSDT2(3)=0.01"
   printf, unit, " ALSDT2(2)=0.01"
   printf, unit, " ALSDT2(1)=0.01"
   printf, unit, " ALEXT2(5)=0.5"
   printf, unit, " ALEXT2(4)=0.5"
   printf, unit, " ALEXT2(3)=0.5"
   printf, unit, " ALEXT2(2)=0.5"
   printf, unit, " ALEXT2(1)=0.5"
   printf, unit, " CHSDT2(5)='Ala'"
   printf, unit, " CHSDT2(4)='Gln'"
   printf, unit, " CHSDT2(3)='Lac'"
   printf, unit, " CHSDT2(2)='Tau'"
   printf, unit, " CHSDT2(1)='Ins'"
   printf, unit, " CHEXT2(5)='Ala'"
   printf, unit, " CHEXT2(4)='Gln'"
   printf, unit, " CHEXT2(3)='Lac'"
   printf, unit, " CHEXT2(2)='Tau'"
   printf, unit, " CHEXT2(1)='Ins'"
   printf, unit, " NSDT2=5"
   printf, unit, " NEXT2=5"
endif
printf, unit, " CHCOMB(3)='Cr+PCr'"	
printf, unit, " CHCOMB(2)='GPC+PCho'"
printf, unit, ' NCOMBI=8'		
printf, unit, " NAMREL='Cr+PCr'"	
printf, unit, ' LPRINT=6'
printf, unit, " FILPRI='lcm.PRI'"
printf, unit, ' LCOORD=9'
ind = where(header.chkeep NE '', count)
if count GT 0 then begin
   chkeep = header.chkeep(ind)
   printf, unit, ' NKEEP=' + string(n_elements(chkeep),format='(I0)')
   for ikeep=0,n_elements(chkeep)-1 do $
      printf, unit, string(ikeep+1, format='(" CHKEEP(",I0,")=")') $
		+ "'" + chkeep(ikeep) + "'"
endif
ind = where(header.chomit NE '', count)
if count GT 0 then begin
   chomit = header.chomit(ind)
   printf, unit, ' NOMIT=' + string(n_elements(chomit),format='(I0)')
   for iomit=0,n_elements(chomit)-1 do $
      printf, unit, string(iomit+1, format='(" CHOMIT(",I0,")=")') $
		+ "'" + chomit(iomit) + "'"
endif
printf, unit, " NAMEAC(5)='GABA'"
printf, unit, " NAMEAC(4)='GSH'"
printf, unit, " NAMEAC(3)='Ins'"
printf, unit, " NAMEAC(2)='Glc'"
printf, unit, " NAMEAC(1)='Mac'"
printf, unit, ' NEACH=' + string(header.neach,format='(I0)')
if header.sddegp GE 0 then $
   printf, unit, ' SDDEGP=' + string(header.sddegp,format='(G0.0)')
printf, unit, ' SHIFMN=' + string(header.shifmn(0),header.shifmn(1), $
			format='(G0.0,",",G0.0)')
if header.ppmshf NE 0 then $
   printf, unit, ' PPMSHF=' + string(header.ppmshf,format='(G0.0)')
printf, unit, ' FWHMBA=' + string(header.fwhmba,format='(G0.0)')
printf, unit, ' RFWHM=' + string(header.rfwhm,format='(G0.0)')
if header.dkntmn NE 0 then $
   printf, unit, ' DKNTMN=' + string(header.dkntmn,format='(G0.0)') $
else begin
   printf, unit, " ECCDON=.TRUE."
   if header.vitro NE 0 then $
      printf, unit, ' VITRO=.TRUE.'
endelse
printf, unit, ' PPMEND=' + string(header.ppmend,format='(G0.0)')
printf, unit, ' PPMST=' + string(header.ppmst,format='(G0.0)')
printf, unit, " FILRAW='" + header.filraw + "'"
printf, unit, " FILBAS='" + header.filbas + "'"
printf, unit, ' DELTAT=' + string(header.deltat,format='(G0.0)')
printf, unit, ' NUNFIL=' + string(header.nunfil,format='(I0)')
printf, unit, ' HZPPPM=' + string(header.hzpppm,format='(G0.0)')
printf, unit, ' ' + NLEND3

endelse	; -------------------- end of common basis handling

free_lun, unit

endfor 	; --- islice

END
;;; ---------------------------------------------------------------------
