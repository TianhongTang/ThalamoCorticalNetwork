function F = funGauss(P,x)

if length(P) < 3
    error('funGauss needs >= 3 arguments!')
end

if P(3) ~= 0
    z = ( x - P(2) )/ P(3);
    ez = exp(-z.*z/2.);	%GAUSSIAN PART
else
    z = 100;
    ez = 0;
end

if length(P) == 3
    F  = P(1).* ez;
elseif length(P) == 4
    F  = P(1).* ez + P(4);
elseif length(P) == 5
    F  = P(1).* ez + P(4) + P(5)*x;
elseif length(P) == 6
    F  = P(1).* ez + P(4) + P(5)*x  + P(6)*x^2;
else
    error('funGauss needs [3..6] arguments!')
end

% PRO	FUNCT,X,A,F,PDER
% ; NAME:
% ;	GAUSS_FUNCT
% ;
% ; PURPOSE:
% ;	EVALUATE THE SUM OF A GAUSSIAN AND A 2ND ORDER POLYNOMIAL
% ;	AND OPTIONALLY RETURN THE VALUE OF IT'S PARTIAL DERIVATIVES.
% ;	NORMALLY, THIS FUNCTION IS USED BY CURVEFIT TO FIT THE
% ;	SUM OF A LINE AND A VARYING BACKGROUND TO ACTUAL DATA.
% ;
% ; CATEGORY:
% ;	E2 - CURVE AND SURFACE FITTING.
% ; CALLING SEQUENCE:
% ;	FUNCT,X,A,F,PDER
% ; INPUTS:
% ;	X = VALUES OF INDEPENDENT VARIABLE.
% ;	A = PARAMETERS OF EQUATION DESCRIBED BELOW.
% ; OUTPUTS:
% ;	F = VALUE OF FUNCTION AT EACH X(I).
% ;
% ; OPTIONAL OUTPUT PARAMETERS:
% ;	PDER = (N_ELEMENTS(X),6) ARRAY CONTAINING THE
% ;		PARTIAL DERIVATIVES.  P(I,J) = DERIVATIVE
% ;		AT ITH POINT W/RESPECT TO JTH PARAMETER.
% ; COMMON BLOCKS:
% ;	NONE.
% ; SIDE EFFECTS:
% ;	NONE.
% ; RESTRICTIONS:
% ;	NONE.
% ; PROCEDURE:
% ;	F = A(0)*EXP(-Z^2/2) + A(3) + A(4)*X + A(5)*X^2
% ;	Z = (X-A(1))/A(2)
% ;	Elements beyond A(2) are optional.
% ; MODIFICATION HISTORY:
% ;	WRITTEN, DMS, RSI, SEPT, 1982.
% ;	Modified, DMS, Oct 1990.  Avoids divide by 0 if A(2) is 0.
% ;	Added to Gauss_fit, when the variable function name to
% ;		Curve_fit was implemented.  DMS, Nov, 1990.
% ;
% 	n = n_elements(a)
% 	ON_ERROR,2                      ;Return to caller if an error occurs
% 	if a(2) ne 0.0 then begin
% 	    Z = (X-A(1))/A(2) 	;GET Z
% 	    EZ = EXP(-Z^2/2.)	;GAUSSIAN PART
% 	endif else begin
% 	    z = 100.
% 	    ez = 0.0
% 	endelse
% 
% 	case n of
% 3: 	F = A(0)*EZ
% 4: 	F = A(0)*EZ + A(3)
% 5: 	F = A(0)*EZ + A(3) + A(4)*X
% 6: 	F = A(0)*EZ + A(3) + A(4)*X + A(5)*X^2 ;FUNCTIONS.
% 	ENDCASE
% 
% 	IF N_PARAMS(0) LE 3 THEN RETURN ;NEED PARTIAL?
% ;
% 	PDER = FLTARR(N_ELEMENTS(X),n) ;YES, MAKE ARRAY.
% 	PDER(*,0) = EZ		;COMPUTE PARTIALS
% 	if a(2) ne 0. then PDER(*,1) = A(0) * EZ * Z/A(2)
% 	PDER(*,2) = PDER(*,1) * Z
% 	if n gt 3 then PDER(*,3) = 1.
% 	if n gt 4 then PDER(*,4) = X
% 	if n gt 5 then PDER(*,5) = X^2
% 	RETURN
% END
