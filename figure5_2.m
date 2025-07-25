%% Reset

clear all;
close all;
clc;

rng(1);

if (ispc)
    folderSep = '\';
else
    folderSep = '/';
end

%% setup

% savepath figures
path_figures = 'figures/boxplots/laplacian/';

%list subjects. Exclude subjects rejecting during processing
sublist = [01 02 03 04 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 31 32 33];
nsubs   = length(sublist);

%% Read PSD data

% load channel layout for later
loadedLayout = load('config/layout.mat');

% channels used for plotting
channelsForFigures = {loadedLayout.layout, ...
    {'AF3','AF4','AFz','F1','F2','F3','F4','Fz'}, ...
    {'C1','C2','C3','C4','Cz'}, ...
    {'O1','O2','Oz','PO3','PO4','POz'}};

path_spectrogram_psd_avg = 'spectrograms/psd/avg/laplacian/all/';

% high reward
filename = sprintf('%s/allSubjects_high.mat',path_spectrogram_psd_avg);
load(filename);
highReward = avgPSD;
% low reward
filename = sprintf('%s/allSubjects_low.mat',path_spectrogram_psd_avg);
load(filename);
lowReward = avgPSD;

% diff psd
subjectAvgDiff = lowReward.psd-highReward.psd;

%% Calculate reaction time differences

% all subjects' reaction times
subLowRTs = {};
subHighRTs = {};

% average PSD for each subject
subjectAvgPSD = [];

path_spectrogram_config = 'spectrograms/config/';

for isub=1:nsubs    
    % load sub config file
    filename = sprintf([path_spectrogram_config 's%d/configpsd.mat'], ...
        sublist(isub));
    if ~exist(filename,'file')
        error(sprintf('%s does not exist!',filename));
    end
    load(filename);
    
    % get trigger indices
    conditionIndex = lib.specific.getConditionIndex( ... 
        subjectData.trigidx,subjectData.trialData); 

    % make variables for RTs
    subLowRTs{isub} = subjectData.RTs(subjectData.RTs.* ...
        (conditionIndex<=8)~=0);
    subHighRTs{isub} = subjectData.RTs(subjectData.RTs.* ...
        (conditionIndex>8)~=0);
end

% average subject reaction time difference between low and high reward
meanDiffRT = zeros(1,length(subHighRTs));
% median subject reaction time difference between low and high reward
medianDiffRT = zeros(1,length(subHighRTs));
% fill in the values
for i = 1:length(subHighRTs)
    meanDiffRT(i) = mean(subLowRTs{i})-mean(subHighRTs{i});
    medianDiffRT(i) = median(subLowRTs{i})-median(subHighRTs{i});
end

% median split on subjects
globalMedianRT = medianDiffRT<median(medianDiffRT);

%% Grab data in desired time window

close all;

bothFreqs = [8 12; 20 30];
avgWindow = [0.5 2.5];
freqNames = {'alpha','beta'};
roiNames = {'All','Frontal','Central','Posterior'};
splitNames = {'HighSplit','LowSplit'};

allRewards = [];
allFreqs = [];
allValues = [];
allSubjects = repmat(1:size(lowReward.psd,1),1, ...
    length(channelsForFigures)*length(freqNames)*2)';
allRegions = [];
allGroups = repmat(globalMedianRT,1, ...
    length(channelsForFigures)*length(freqNames)*2)';

for j = 1:length(channelsForFigures)
    psdTime = (lowReward.t>=avgWindow(1)) & (lowReward.t<=avgWindow(2));
    psdAlpha = (lowReward.f>=bothFreqs(1,1)) & (lowReward.f<=bothFreqs(1,2));
    psdBeta = (lowReward.f>=bothFreqs(2,1)) & (lowReward.f<=bothFreqs(2,2));

    [~,selectedChannels,~] = intersect(loadedLayout.layout, ...
        channelsForFigures{j});

    lowRewardAlpha = mean(mean(mean(lowReward.psd(:, ...
        selectedChannels,psdAlpha,psdTime),4),3),2);
    lowRewardBeta = mean(mean(mean(lowReward.psd(:, ...
        selectedChannels,psdBeta,psdTime),4),3),2);
    highRewardAlpha = mean(mean(mean(highReward.psd(:, ...
        selectedChannels,psdAlpha,psdTime),4),3),2);
    highRewardBeta = mean(mean(mean(highReward.psd(:, ...
        selectedChannels,psdBeta,psdTime),4),3),2);

    allValues = [allValues;lowRewardAlpha;lowRewardBeta; ...
        highRewardAlpha;highRewardBeta];

    allFreqs = [allFreqs;zeros(size(lowRewardAlpha,1),1); ...
        ones(size(lowRewardAlpha,1),1);zeros(size(lowRewardAlpha,1),1); ...
        ones(size(lowRewardAlpha,1),1)];
    allRewards = [allRewards;zeros(size(lowRewardAlpha,1),1); ...
        zeros(size(lowRewardAlpha,1),1);ones(size(lowRewardAlpha,1),1); ...
        ones(size(lowRewardAlpha,1),1)];
    allRegions = [allRegions;j*ones(4*size(lowRewardAlpha,1),1)];
end

dataTable = table;
dataTable.Subject = allSubjects;
dataTable.ROI = allRegions;
dataTable.Reward = allRewards;
dataTable.FrequencyBand = allFreqs;
dataTable.AveragePower = allValues;
dataTable.SubjectSplit = allGroups;

% check for folder existence
if (exist(path_figures,'file')==0)
    mkdir(path_figures);
end

for i = 1:length(roiNames)
    for j = unique(allGroups)'
        selectedTable = dataTable(dataTable.ROI==i & ...
            dataTable.SubjectSplit==j,:);
        
        figure; hold on;
        myBox = boxchart(selectedTable.FrequencyBand, ...
            selectedTable.AveragePower,'GroupByColor', ...
            selectedTable.Reward,'MarkerStyle','.');

        for k = unique(selectedTable.FrequencyBand)'
            freqTableLow = selectedTable(selectedTable.FrequencyBand==k & ...
                selectedTable.Reward==0,:).AveragePower;
            freqTableHigh = selectedTable(selectedTable.FrequencyBand==k & ...
                selectedTable.Reward==1,:).AveragePower;
            
            [~,p,~,stats] = ttest(freqTableLow,freqTableHigh);

            addSignificanceBar(p,[k 2.5],2.2);
        end

        legend({'Low reward','High reward','',''},'Location','southeast', ...
            'FontSize',18); legend boxoff
        xticks(unique(dataTable.FrequencyBand));
        set(gca,'XTickLabel',{'Alpha','Beta'});
        ylim([-7 3]);
        
        figureName = [splitNames{j+1} roiNames{i}];
        figurePath = [path_figures figureName];
        print(gcf,'-dpng',figurePath);
        print(gcf,'-dpdf','-painters',figurePath);

        close all;
    end
end

%% Additional functions

function addSignificanceBar(pvalue,textPosition,barYPosition)
    significanceText = '';
    if (pvalue<0.0001)
        significanceText = '****';
    elseif (pvalue<0.001)
        significanceText = '***';
    elseif (pvalue<0.01)
        significanceText = '**';
    elseif (pvalue<0.05)
        significanceText = '*';
    else
        significanceText = 'ns';
        return;
    end

    text(textPosition(1),textPosition(2),significanceText, ...
        "FontSize",18,"HorizontalAlignment","center");
    plot([textPosition(1)-1/3 textPosition(1)+1/3], ...
        barYPosition*[1 1],'k','LineWidth',1.5);
end