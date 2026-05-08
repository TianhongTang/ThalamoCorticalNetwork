function MUA = processMUAFile(adffile,muafile)
global FILES

MUA = subProcess(adffile,muafile);  % adfw file

%----------------------------------------------
function MUA = subProcess(adffile,muafile)
[nchan,nobs,sampt] = adfw_info(adffile);

wb = waitbar(0, sprintf('Processing Multi Data'));
valobs = 1;
for obs=1:nobs
    waitbar((obs/nobs),wb);
    for j=1:nchan
        %[tmp, resampt] =  readMultiObs(adffile,j,obs-1);
        % changed on 29-01-03 AM
        % [seems like the original processMUAFile was
        % deleted in the utils directory (???)
        % since that one did not have the bug !!! -
        % see Melanies Computer ...]
        [tmp, resampt] =  readMultiObs(adffile,j,obs);
        MUA.data{valobs}(:,j) = tmp;     
        lastlen = length(tmp);
    end
    valobs = valobs+1;
end
MUA.resampt = resampt;	
close(wb);

%
% Here just save the thing for future use
save(muafile,'MUA');




