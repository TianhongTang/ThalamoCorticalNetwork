function PVconvRF(path, filein, fileout, rfversion)

if nargin > 0
    PVconvRFfun(path, filein, fileout, rfversion)
else

	brukerPath = '//wks5/guest/mpi/PV_install/RF_bruker/';
	%PVconvRF(brukerPath,'bp32.pls','bp32', 1);
	%PVconvRF(brukerPath,'gauss256.pls','gauss256', 1);
	%PVconvRF(brukerPath,'hb_lurie_ex.pls','hb_lurie_ex', 1);
	%PVconvRF(brukerPath,'wd_mao_4_256.pls','wd_mao_4_256', 1);
	%PVconvRF(brukerPath,'wd_mao_6_256.pls','wd_mao_6_256', 1);
	
	vnmrPath = '//wks5/guest/mpi/PV_install/RF_vnmr/';
	%PVconvRF(vnmrPath,'M8.RF','M8', 2);
	%PVconvRF(vnmrPath,'M88.RF','M88', 2);
	%PVconvRF(vnmrPath,'M8sym.RF','M8sym', 2);
	%PVconvRF(vnmrPath,'P10.RF','P10', 2);
	%PVconvRF(vnmrPath,'P10tr.RF','P10tr', 2);
	%PVconvRF(vnmrPath,'P10tr1.RF','P10tr1', 2);
	%PVconvRF(vnmrPath,'P10tr2.RF','P10tr2', 2);
	%PVconvRF(vnmrPath,'P11.RF','P11', 2);
	%PVconvRF(vnmrPath,'at60ph0.RF','at60ph0', 2);
	%PVconvRF(vnmrPath,'bir-54.RF','bir-54', 2);
	%PVconvRF(vnmrPath,'hs10.RF','hs10', 2);
	%PVconvRF(vnmrPath,'hs20.RF','hs20', 2);
	%PVconvRF(vnmrPath,'invpat.2.RF','invpat.2', 2);
	%PVconvRF(vnmrPath,'invpat.25.RF','invpat.25', 2);
	% PVconvRF(vnmrPath,'at60.n29.RF','at60.n29', 2);
	% PVconvRF(vnmrPath,'at60ph0.RF','at60ph0', 2);
	% PVconvRF(vnmrPath,'invpat.10.RF','invpat.10', 2);
	% PVconvRF(vnmrPath,'invpat.15.RF','invpat.15', 2);
	% PVconvRF(vnmrPath,'invpat.20.RF','invpat.20', 2);
	% PVconvRF(vnmrPath,'invpat.25.RF','invpat.25', 2);
	%PVconvRF(vnmrPath,'scewd.RF','scewd', 2);
	%PVconvRF(vnmrPath,'scewdrev.RF','scewdrev', 2);
       % FASTESTMAP pulses
	% PVconvRF(vnmrPath,'at60.n29.RF','at60.n29', 2);
	% PVconvRF(vnmrPath,'bir4.45.RF','bir4.45', 2);
	% PVconvRF(vnmrPath,'invpat.10.RF','invpat.10', 2);
	% PVconvRF(vnmrPath,'invpat.2.RF','invpat.2', 2);
	% PVconvRF(vnmrPath,'invpat2.10.RF','invpat2.10', 2);
	% PVconvRF(vnmrPath,'invpat2.2.RF','invpat2.2', 2);
	
end
return
%--------------------------------------------------------------
function PVconvRFfun(path, filein, fileout, rfversion)
%
% rfversion     1:Bruker, 2: Varian 
%

[RFamp, RFphase] = PVrdRF([path filein], rfversion);
RFnpts = length(RFamp);
RFarr = zeros(2,RFnpts);
RFarr(1,:) = RFamp';
RFarr(2,:) = RFphase';

shapeHeader = {
	'##TITLE= /usr/local/mpi/PV_install/RF_bruker/';
	'##JCAMP-DX= 5.00 Bruker JCAMP library';
	'##DATA TYPE= Shape Data';
	'##ORIGIN= Max-Planck Institute for Biological Cybernetics';
	'##OWNER= <josef>';
	'##DATE= 03/01/25';
	'##TIME= 10:43:01';
	'##MINX= 0.000000e+00';
	'##MAXX= 0.000000e+00';
	'##MINY= 0.000000e+00';
	'##MAXY= 0.000000e+00';
	'##$SHAPE_EXMODE= Excitation';   % 'Excitation,Inversion,Refocussing'
	'##$SHAPE_TOTROT= 0';            % flipangle used to calibrate SHAPE_BWFAC
	'##$SHAPE_BWFAC= 0';             % BW FWHM [Hz*s]
	'##$SHAPE_INTEGFAC= 0';          % area compared to rectangular pulse
	'##$SHAPE_REPHFAC= 0';           % dephasing property: [0-100%] for Exc, [-100-100%] for Inv,Rfc
	'##$SHAPE_TYPE= conventional';   % 'conventional, adiabatic'        
	'##$SHAPE_MODE= 0';              % only 0 allowed
	['##NPOINTS= ' num2str(RFnpts)];
	'##XYPOINTS= (XY..XY)';
};
shapeFooter = '##END=';

fp = fopen([path fileout], 'w');
for iarr=1:size(shapeHeader)
   istr = char( shapeHeader(iarr));
   fprintf(fp,'%s\n', istr);
end
fprintf(fp,'%g, %g\n', RFarr);
fprintf(fp,'%s\n', shapeFooter);
fclose(fp);

disp( sprintf('--- RF pulse <%s> written in JCAMP format\n', [path fileout]) );
key

return
%--------------------------------------------------------------
function [RFamp, RFphase] = PVrdRF(file, rfversion)

disp( sprintf('--- RF pulse <%s> read\n', file) );

RFdat = load(file);
whos RFdat

if rfversion == 1           % Bruker format
   RFamp = RFdat(:,1);
   RFphase = RFdat(:,2);
else                        % Varian format
   RFamp = RFdat(:,2);
   RFphase = RFdat(:,1);
end

figure(1);
colormap(gray);
plot(RFamp);
figure(3);
colormap(gray);
plot(RFphase);

    % scale RFamp to 100%
RFamp = RFamp/max(RFamp)*100;

return
%--------------------------------------------------------------
