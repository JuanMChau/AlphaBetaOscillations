function [outData,outT] = baselineCorrection(data,t,toffset, ...
    tbaseline,logtransform,type)
%BASELINECORRECTION Summary of this function goes here
%   Detailed explanation goes here

outT = t-toffset;

if (logtransform)
    myspec = 10*log10(data);
else
    myspec = data;
end

meanData = mean(myspec(:,:,outT>=tbaseline(1) & outT<=tbaseline(2)),3);
stdData = std(myspec(:,:,outT>=tbaseline(1) & outT<=tbaseline(2)),[],3);

switch(type)
    case 0
        outData = myspec;
    case 1
        outData = 100*(myspec-meanData)./meanData;
    case 2
        outData = (myspec-meanData)./stdData;
    case 3
        outData = myspec-meanData;
end

end

