function jp210a
%
global STDPATH
STDPATH.pv = 'M:/mridata/'     % wks5

M02gs1

%--------------------------------------------------------------
function M02gs1
% demo of DORK correction: dynamic offresonance in k-space
% uses navigator FID line or internal navigator echo (center k-space) to determine
%      **global** offresonance changes in the image.
%      these changes make a pixel shift in PE direction
% top: DORK frequency extracted from data
% middle: Center of mass before/after DORK correction
%         the 10 Hz offresonance shift makes about 0.2 pixel shift
% bottom: using the same pixels as mask: BOLD time course before/after
%       MAY BE better to show relative changes after normalization!!

global STDPATH
f_create = 0
filename = 'M02.gs1_28_DORKdemo.mat'

if f_create 
	CenterofMassPE = load('M:/mridata_wks8/M02.gs1/28_massCtr.asc');
	CenterofMassPE = reshape(CenterofMassPE(:,2), 484, 2);
	BOLDtimecourse = reshape( load('M:/mridata_wks8/M02.gs1/28_avgTC.asc'), 484, 2);
	
	STDPATH.pv = 'M:/mridata_wks8/';
	deltaFreq = dorkTrace('M02.gs1', 28, opt('NR',[53 0],'REFNUM',0,'FIDFILE','fid.orig','VERBOSE',0) );
	deltaFreq = squeeze( avg(deltaFreq.deltaFreq, 1));
	
	key('SAVE ? ')
	save(filename,'CenterofMassPE','BOLDtimecourse','deltaFreq');
	fprintf('\n<%s> saved\n', filename)
else
	load(filename)

    subplot(3,2,1)
	plot(deltaFreq)
	ylim([-3 11])
	ylabel('Dynamic Offresonance Freq (DORK) / Hz')
	subplot(3,2,3)
	plot(CenterofMassPE(:,1))
	ylim([15.5 15.8])
	ylabel('Center of Mass / Voxeldim')
	subplot(3,2,4)
	plot(CenterofMassPE(:,2))
	ylim([15.5 15.8])
	ylabel('Center of Mass / Voxeldim')
	subplot(3,2,5)
	plot(BOLDtimecourse(:,1))
	ylim([1250 1460])
	ylabel('BOLD signal change / a.u.')
	subplot(3,2,6)
	plot(BOLDtimecourse(:,2))
	ylim([1250 1460])
	ylabel('BOLD signal change / a.u.')
end

%--------------------------------------------------------------
