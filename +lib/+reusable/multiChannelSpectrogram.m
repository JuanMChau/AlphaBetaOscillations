function [outS,f,t,outPSD] = multiChannelSpectrogram(data,window, ...
    noverlap,nfft,fs)
%MULTICHANNELSPECTROGRAM Summary of this function goes here
%   Detailed explanation goes here

% It is assumed that data comes in the format: channel x time, data will be
% returned in the format: channel x frequency x time

for i = 1:size(data,1)
    % Calculate spectrogram/psd
    [s,f,t,psd] = spectrogram(data(i,:),window,noverlap,nfft,'yaxis',fs);
    
    % Allocate memory space if needed
    if (i==1)
        outS = zeros([size(data,1) size(s)]);
        outPSD = zeros([size(data,1) size(psd)]);
    end
    
    % Fit new data into preallocated slots
    outS(i,:,:) = s;
    outPSD(i,:,:) = psd;
end

end

