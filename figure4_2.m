% Author: Juan M. Chau
% Last uploaded: 01.03.2023

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

% TODO: edit this path to your fieldtrip folder
useFieldtrip = '../../fieldtrip/';
lib.fieldtripStartup(useFieldtrip);

%% Initial variables

% apply log transform to psd
logTransform = 1;

% set value of baseline correction used (0, 1, 2, 3)
% 0: no baseline
% 1: percent of variance -> (data-mean)/mean
% 2: zscore -> (data-mean)/std
% 3: center around 0 -> data-mean
freqBaseline = 3;

% load channel layout for later
loadedLayout = load('config/layout.mat');

% groups of channels to consider for regression
[~,selectedChannels,~] = intersect(loadedLayout.layout,loadedLayout.layout);

% list of possible triggers to consider
allTriggers = {9:16, 1:8};

% cluster correction permutations
ccperms = 1e4;

% cluster correction significance level
ccpvalue = 0.05;

% savepath spectrogram
path_spectrogram_psd_avg = 'spectrograms/psd/avg/laplacian/all/';

% savepath figures
path_figures = 'figures/handmade_from_psd/laplacian/all/';

% custom colormap for graphics
myColormap = lib.graphics.brewermap.brewermap([],'RdBu');
myColormap = myColormap(end:-1:1,:,:);
close all;

% graphics scaling limits
xlims = [0.5 1.0];
ylims = [8 12];
zlims = [-1 1];

% graphics scaling limits (tstat)
tstatlims = [-1 1]*4;

%% Load the data
filename = [path_spectrogram_psd_avg 'allSubjects_low.mat'];
if ~exist(filename,'file')
    error(sprintf('%s does not exist!',filename));
end
load(filename);

avgPSDLow = avgPSD;

filename = [path_spectrogram_psd_avg 'allSubjects_high.mat'];
if ~exist(filename,'file')
    error(sprintf('%s does not exist!',filename));
end
load(filename);

avgPSDHigh = avgPSD;

%% Display on fieldtrip

% check for figure folder existence
if (exist(path_figures,'file')==0)
    mkdir(path_figures);
end

tfData = [];
tfData.label = loadedLayout.layout(selectedChannels);
tfData.dimord = 'subj_chan_freq_time';
tfData.freq = avgPSD.f;
tfData.time = avgPSD.t;
tfData.powspctrm = avgPSDHigh.psd - avgPSDLow.psd;

%% Just run this bit to save the figure as it is

% channels used for plotting
channelsForFigures = {{'AF3','AF4','AFz','F1','F2','F3','F4','Fz'}, ...
    {'C1','C2','C3','C4','Cz'}, ...
    {'O1','O2','Oz','PO3','PO4','POz'}};

% trigger times
triggerTimes = [0.0 1.2;
                1.2 1.8;
                1.8 4.5];

% area names for the channel groups
areaNames = {'frontal','central','posterior'};

labelPositionsX = {[2.2],[0.85,1.7],[1.2,3.2]};
labelPositionsY = {[20],[29.5,6],[18.5,21]};

for i = 1:length(channelsForFigures)
    % select desired channels
    [~,selectedChannels,~] = intersect(loadedLayout.layout, ...
        channelsForFigures{i});
    selectedChannels = sort(selectedChannels);
    
    %% work with the psd (all of it)
    % calculate cluster correction for all psd image
    [p,praw] = lib.fromToolbox.stat.ClusterCorrection2(squeeze(mean(tfData.powspctrm(:, ...
        selectedChannels,:,:),2)),ccperms,ccpvalue);
    relevantP = ~isnan(squeeze(p)) & squeeze(p)<0.05;

    myFig = figure;
    imageToDisplay = squeeze(mean(mean(tfData.powspctrm(:, ...
        selectedChannels,:,:),1),2));
    ft_plot_matrix(avgPSD.t,avgPSD.f,imageToDisplay, ...
        'clim',zlims,'tag','cip','highlightstyle','outline', ...
        'highlight',relevantP);
    colormap(myColormap); hold on; axis xy; colorbar;

    % figure beautification
    for j = 1:size(triggerTimes,1)
        plot(triggerTimes(j,1)*[1 1],[avgPSD.f(1) avgPSD.f(end)], ...
            '--','LineWidth',2,'Color','#888888');
    end

    clusterP = [];
    clusterLabels = zeros(size(relevantP));

    labels = bwlabel(relevantP);
    if (sum(labels(:))~=0)
        allLabels = unique(labels(:));
        allLabels = allLabels(2:end)';
        for j = allLabels
            testCluster = labels==j;
            testP = mean(p(logical(testCluster)));
            if (testP<0.05)
                clusterP = [clusterP testP];
                clusterLabels(logical(testCluster)) = j;
            end
        end
    end

    for j = 1:length(clusterP)
        text(labelPositionsX{i}(j),labelPositionsY{i}(j), ...
            lib.text.getPValueText(clusterP(j)), ...
            'HorizontalAlignment','center','FontSize',16);
    end

%     myFig.WindowState = 'maximized';

    % save the figure
    selectedAreaName = areaNames{i};
    figureName = sprintf('highMinusLow_%s',selectedAreaName);
    figurePath = [path_figures figureName];
    print(gcf,'-dpng',figurePath);
    print(gcf,'-dpdf','-painters',figurePath);
    close(myFig);
end
