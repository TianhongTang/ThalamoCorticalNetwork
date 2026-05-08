function f = funGamma(a,x)
% f = funGauss(a,x):
% --- Cohen 1997 -- 
%     Gamma variate impulse model G(t, gama_a, gam_b, gam_c)
%   gam_a = 0.008;
%   gam_b = 8.60;
%   gam_c = 0.547;
%
%   JP Apr 2002

if length(a) < 3
   error('funGamma: length(a) < 3 !')
end

f_verbose = 0;

gam_a = a(1);
gam_b = a(2);
gam_c = a(3);

f = gam_a .* (x.^gam_b).* exp(-x/gam_c);
f = f/(sum(f)*(x(2)-x(1)));

if ( f_verbose )
    disp( sprintf('Gamma model: %f %f',gam_b, gam_c));
    disp( sprintf('Gamma model: time-to-peak [s]= %f', gam_b*gam_c) );
    disp( sprintf('Gamma model: FWHM [s]~ %f', 2.35 * sqrt(gam_b) * gam_c) );
    plot(time,f);
end
