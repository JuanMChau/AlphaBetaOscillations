function [outS,f,t,outPSD] = multiTrialSpectrogram(data,window, ...
    noverlap,nfft,fs,returnRaw)
%MULTITRIALSPECTROGRAM Summary of this function goes here
%   Detailed explanation goes here

% It is assumed that data is recorded in the format: trial x channel x
% time, data will be returned in the format: trial x channel x frequency x
% time

data = permute(data,[2 3 1]);

for i = 1:size(data,3)
    % Display iteration for easier understanding of what's going on
    disp(['Calculating step ' num2str(i) '/' num2str(size(data,3))]);
    
     % Calculate spectrogram/psd
    [s,f,t,psd] = lib.reusable.multiChannelSpectrogram(data(:,:,i), ...
            window,noverlap,nfft,fs);
    
    % Allocate memory space if needed
    if (i==1)        
        outS = zeros([size(data,3) size(s)]);
        outPSD = zeros([size(data,3) size(psd)]);
    end
    
    % Fit new data into preallocated slots
    outS(i,:,:,:) = s;
    outPSD(i,:,:,:) = psd; 
end

% Remove one dimension if needed: channel x time x frequency
if (~returnRaw)
    outS = squeeze(mean(outS,1));
    outPSD = squeeze(mean(outPSD,1));
end

end

