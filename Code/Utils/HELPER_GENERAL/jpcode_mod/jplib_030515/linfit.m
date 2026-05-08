;$Id: linfit.pro,v 1.6 1997/01/15 03:11:50 ali Exp $
;
; Copyright (c) 1994-1997, Research Systems, Inc.  All rights reserved.
;       Unauthorized reproduction prohibited.
; Adapted for PV-WAVE by JP - 07/1998
;+
; NAME:
;       LINFIT
;
; PURPOSE:
;       This function fits the paired data {X(i), Y(i)} to the linear model,
;       y = A + Bx, by minimizing the chi-square error statistic. The result
;       is a two-element vector containing the model parameters,[A,B].
;
; CATEGORY:
;       Statistics.
;
; CALLING SEQUENCE:
;       Result = LINFIT(X, Y)
;
; INPUTS:
;       X:    An n-element vector of type integer, float or double.
;
;       Y:    An n-element vector of type integer, float or double.
;
; KEYWORD PARAMETERS:
;   CHISQ:    Use this keyword to specify a named variable which returns the
;             chi-square error  as the sum of squared errors between
;             Y(i) and A + BX(i). If individual standard deviations are 
;             supplied, then the chi-square error statistic is computed as
;             the sum of squared errors divided by the standard deviations.
;
;  ;;DOUBLE:    If set to a non-zero value, computations are done in double 
;  ;;           precision arithmetic.
;
;    PROB:    Use this keyword to specify a named variable which returns the
;             probability that the computed fit would have a value of CHISQR 
;             or greater. If PROB is greater than 0.1, the model parameters 
;             are "believable". If PROB is less than 0.1, the accuracy of the
;             model parameters is questionable.
;
;    SDEV:    An n-element vector of type integer, float or double that 
;             specifies the individual standard deviations for {X(i), Y(i)}.
;
;   SIGMA:    Use this keyword to specify a named variable which returns a 
;             two-element vector of probable uncertainties for the model par-
;             ameters, [SIG_A,SIG_B].
;
;  NOPLOT:    omit plot and printout of results
;
; EXAMPLE:
;       Define two n-element vectors of paired data.
;         x = [-3.20, 4.49, -1.66, 0.64, -2.43, -0.89, -0.12, 1.41, $
;               2.95, 2.18,  3.72, 5.26]
;         y = [-7.14, -1.30, -4.26, -1.90, -6.19, -3.98, -2.87, -1.66, $
;              -0.78, -2.61,  0.31,  1.74]
;       Define a vector of standard deviations with a constant value of 0.85
;         sdev = replicate(0.85, n_elements(x))
;       Compute the model parameters, A and B.
;         result = linfit(x, y, chisq = chisq, prob = prob, sdev = sdev)
;       The result should be the two-element vector:
;         [-3.44596, 0.867329]
;       The keyword parameters should be returned as:
;         chisq = 11.4998, prob = 0.319925
;
; REFERENCE:
;       Numerical Recipes, The Art of Scientific Computing (Second Edition)
;       Cambridge University Press
;       ISBN 0-521-43108-5
;
; MODIFICATION HISTORY:
;       Written by:  GGS, RSI, September 1994
;                    LINFIT is based on the routines: fit.c, gammq.c, gser.c,
;                    and gcf.c described in section 15.2 of Numerical Recipes,
;                    The Art of Scientific Computing (Second Edition), and is
;                    used by permission.
;         Modified:  SVP, RSI, June 1996
;		     Changed SIG_AB to SIGMA to be consistant with the other
;		     fitting functions. Changed CHISQR to CHISQ in the docs
;                    for the same reason. Note that the chisqr and the SIG_AB
;		     keywords are left for backwards compatibility.
;         Modified:  GGS, RSI, October 1996
;                    Modified keyword checking and use of double precision.
;                    Added DOUBLE keyword.
;-

FUNCTION LinFit, x, y, chisqr = chisqr, prob = prob, NOPLOT=f_noplot, $
                       sdev = sdev, sig_ab = sig_ab, sigma = sigma

  ON_ERROR, 2

  if NOT keyword_set(f_noplot) then $
	f_noplot = 0

  TypeX = SIZE(X)
  TypeY = SIZE(Y)
  nX = TypeX(TypeX(0)+2)
  nY = TypeY(TypeY(0)+2)

  if nX ne nY then $
    MESSAGE, "X and Y must be vectors of equal length."

  ;If the DOUBLE keyword is not set then the internal precision and
  ;result are identical to the type of input.
  ;;if N_ELEMENTS(Double) eq 0 then $
  ;;  Double = (TypeX[TypeX[0]+1] eq 5 or TypeY[TypeY[0]+1] eq 5)

  nsdev = n_elements(sdev)

  if nsdev eq nX then begin ;Standard deviations are supplied.
    wt = 1.0 / sdev^2
    ss = TOTAL(wt)
    sx = TOTAL(wt * x)
    sy = TOTAL(wt * y)
    t =  (x - sx/ss) / sdev
    st2 = TOTAL(t^2)
    b = TOTAL(t * y / sdev)
  endif else if nsdev eq 0 then begin
    ss = nX + 0.0
    sx = TOTAL(x)
    sy = TOTAL(y)
    t = x - sx/ss
    st2 = TOTAL(t^2)
    b = TOTAL(t * y)
  endif else $
    MESSAGE, "sdev and x must be vectors of equal length."

  st2 = FLOAT(st2) & b = FLOAT(b)

  b = b / st2
  a = (sy - sx * b) / ss
  sdeva = SQRT((1.0 + sx * sx / (ss * st2)) / ss)
  sdevb = SQRT(1.0 / st2)

  if nsdev ne 0 then begin
    chisqr = TOTAL( ((y - a - b * x) / sdev)^2 )
    chisqr = FLOAT(chisqr)
    prob = 1 - IGAMMA(0.5*(nX-2), 0.5*chisqr)
  endif else begin
    chisqr = TOTAL( (y - a - b * x)^2 )
    chisqr = FLOAT(chisqr)
    prob = chisqr * 0 + 1 ;Make prob same type as chisqr.
    sdevdat = SQRT(chisqr / (nX-2))
    sdeva = sdeva * sdevdat
    sdevb = sdevb * sdevdat
  endelse

  sig_ab = [sdeva, sdevb]
  sigma = sig_ab

  if f_noplot LE 0 then begin
    if f_noplot EQ 0 then begin
       print,string(b,sdevb, format="('slope: ',G16.8,' +/- ',G16.8)" )
       print,string(a,sdeva, format="('A0   : ',G16.8,' +/- ',G16.8)" )
    endif    

    pmax = max([a + b*x,y], min=pmin)	;pmin,pmax - 8.5,9.5
    plot,x,a + b*x, xtickformat='(I0)',yrange=[pmin,pmax]  ;;,xrange=[0,500]
    oplot,x,y,psym=4,symsize=0.5
  endif

  RETURN, FLOAT([a, b])

END
