function B0map = procPhaseMap(dir,infile1,infile2,outfile)
%     procPhaseMap('E02.lb1',[77 2],[76 2],[76 3])
%
%
%  calc B0 map of FLASH images with different TE
%  infile must be COMPLEX
%
%STDPATH.pv = '//wks5/josef/mridata/' 
%procPhaseMap('E02.lb1',[77 2],[76 2],[76 3])
%
%
global acqp reco STDPATH

f_write = 1;
f_verbose = 1;
FITTHRES = 5;    %% SD of noise
noiseInd = [1 1 10 10];

img1 = PVrd2dseq(dir, infile1(1), opt('RECO',infile1(2)));
TE1 = acqp.PVM_EchoTime/1000;
img2 = PVrd2dseq(dir, infile2(1), opt('RECO',infile2(2)));
TE2 = acqp.PVM_EchoTime/1000;

dimg = ( double(img1) ).*( conj(double(img2)) );
dfreq = angle(dimg) / ((TE1 - TE2)*2*pi);  % [Hz]
B0map = dfreq;

% noiseDat = reshape( double(idat(noiseInd(1):noiseInd(3),noiseInd(2):noiseInd(4),1,1)), ...
%     (noiseInd(3)-noiseInd(1)+1)*(noiseInd(4)-noiseInd(2)+1),1 );
% noisemean = mean(noiseDat);
% noisestd  = std(noiseDat);
% imageThres = noisemean + FITTHRES*noisestd;
%     imagesc(double(idat(:,:,refImage,1)) > imageThres)

s_dfreq = size(dfreq);
if f_verbose
    imagesc(dfreq(:,:,ceil(s_dfreq(3)/2))');
end

if f_write
    img3 = int16(dfreq*10);   % accuracy 0.1 Hz
    STIMwrSdt(img3,sprintf('%s%s/%d_dfreq',STDPATH.pv,dir,outfile(1)),opt('ROTATE',0))

    if f_write > 1
        info = PVrd2dseq(dir, outfile(1), opt('RECO',outfile(2),'GETINFO',1));
        savename = sprintf('%s.dfreq',info.file);
        fid = fopen(savename, 'w', info.byteorder);
        if (f_verbose)
            fprintf('writing <%s>\n', savename);
        end
        fcount = fwrite(fid, img3, info.precision); 
        fclose(fid);
        flen = info.nx*info.ny*info.nslices*info.nr;
        if (fcount < flen)
            error(sprintf('%s: only %d / %d written', FCTNAME, fcount, flen))
        end
    end
end

%------------------------------------------------------------
% b0map()		/* recon first two images, calc field, fit plane */
% {
%   int i, j, k;
%   float img1[MXNIM][MXNIM][2], img2[MXNIM][MXNIM][2]; 
%   float dmi, dmr, dmr1, dmi1, dmr2, dmi2, wt;
%   float x, y;
%   double sxx, syy, sxy, sx, sy, sxw, syw, sw, snn, den;
%   float scale, fscale;
%   char str[80];
% 
%   printf("Calculating Bo map .. "); fflush(stdout);
%   scale = 1.0/(mapdt*.000001);		/* freq in radians/s  */
%   fscale = scale/(2*M_PI);
% 
% /*  make the complex images */
% 
%   foff = wx = wy = 0;
%   set_traj();		/* reset to no correction  */
%   load_ft(0);		/* read data and grid  */
%   ft2d();		/*  do fft */
%   for (i=0; i<npix; i++)
%     for (j=0; j<npix; j++)  {
%       img1[i][j][0] = img[i][j][0];
%       img1[i][j][1] = img[i][j][1];
%     }
%   load_ft(0);		/* read data and grid  */
%   ft2d();		/*  do fft */
%   for (i=0; i<npix; i++)
%     for (j=0; j<npix; j++)  {
%       dmr2 = img[i][j][0];
%       dmi2 = img[i][j][1];
%       dmr1 = img1[i][j][0];
%       dmi1 = img1[i][j][1];
%       dmr = dmr1*dmr2 + dmi1*dmi2;
%       dmi = dmi1*dmr2 - dmr1*dmi2;
%       img2[i][j][0] = hypot(dmi2, dmr2);
%       img2[i][j][1] = (dmr != 0)? atan2(dmi, dmr) : (dmi >= 0)? M_PI : -M_PI;
%       img[i][j][0] = fscale*img2[i][j][1];
%     }
%   if(bmapflag==2)  {	/* give P.xxxxx.B0, also save the mag images */
%     sprintf(str, "%s.B0", pfile);
%     write_image(str, isl+1, REAL+NOSUBSC);
%     for (i=0; i<npix; i++)
%       for (j=0; j<npix; j++)  
%         img[i][j][0] = img2[i][j][0];
%     write_image(pfile, isl+1, REAL+NOSUBSC);
%   }
%   else 
%     write_image("B0", isl+1, REAL+NOSUBSC);
% 
% }
% 
