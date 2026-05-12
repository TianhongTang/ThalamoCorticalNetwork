function [morletLFP] = convMorletWaveform(inputLFP, morletParameters)

%% Setup parameters & LFP data
% Ephys Parameters
samplingFreq = morletParameters.samplingFreq;
frequencies = morletParameters.frequencies;
cycles = morletParameters.cycle;

%% Extract beta activity for each trial.
% For each beta frequency
for freqIdx = 1:length(frequencies)
    
    % Clear data for loop
    clear sine_wave gaussian_win normalization_factor wavelet halfwaveletsize...
        n_conv fft_w fft_e ift wavelet_conv_data
    
    % Frequency and cycle parameters
    f = frequencies(freqIdx);
    s = cycles/(2*pi*f);
%     fprintf('Extracting %iHz beta activity from LFP...  \n', f)
    
    % Wavelet setup
    time = -1:1/samplingFreq:1;
    sine_wave = exp(1i*2*pi*f.*time);
    gaussian_win = exp(-time.^2./(2*s^2));
    normalization_factor = 1 / (s * sqrt(2* pi));
    wavelet = normalization_factor .* sine_wave .* gaussian_win;
    halfwaveletsize = ceil(length(wavelet)/2); % half of the wavelet size
    n_conv = length(wavelet) + size(inputLFP,2) - 1; % compute Gaussian
    
    % Convolve data and perform fast fourier transform
    fft_w = fft(wavelet,n_conv);
        
    % For each trial
    for trl = 1:size(inputLFP,1)
        % Clear loop variabls
        clear fft_e ift
        
        % Perform fourier transform
        fft_e = fft(inputLFP(trl,:),n_conv);
        ift   = ifft(fft_e.*fft_w,n_conv);
        morletLFP(trl,:,freqIdx) = abs(ift(halfwaveletsize:end-halfwaveletsize+1)).^2;
    end
    
end

end
