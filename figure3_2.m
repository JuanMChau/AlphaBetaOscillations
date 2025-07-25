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

% set value of baseline correction used (0, 1, 2, 3)
% 0: no baseline
% 1: percent of variance -> (data-mean)/mean
% 2: zscore -> (data-mean)/std
% 3: center around 0 -> data-mean
freqBaseline = 3;

% list subjects. Exclude subjects rejecting during processing
sublist = [1:4 6:27 30:33];
% sublist = [01 02 03 04 06 07];
% sublist = 7;
nsubs = length(sublist);

% savepath spectrogram
path_regression = 'regression/data_rsa_all/';

% savepath figures
path_figures = 'figures/regression_rsa/all/';

% image scaling limits
minValue = -0.1;
maxValue = 0.1;

% channel groups
channelGroups = 4;
channelGroupLabels = {'all','posterior','central','frontal'};

% regression coefficients
coefficientSelectors = sort(repmat(2:6,1,channelGroups));

% files to process
windowSelectors = [1:8 repmat(9:12,1,3)];

% cluster correction permutations
ccperms = 1e4;

% cluster correction significance level
ccpvalue = 0.05;

% custom colormap for graphics
myColormap = lib.graphics.brewermap.brewermap([],'RdBu');
myColormap = myColormap(end:-1:1,:,:);
close all;

% figure titles
figureTitles = {'Reward coding', ...
                'Task coding', ...
                'Task-relevant feature coding', ...
                'Task-irrelevant feature coding', ...
                'Motor coding'};

% relevant time windows for analysis
allTimeWindows = [0 4.5; 1.2 4.5; 1.8 4.5; 1.8 4.5; 1.8 4.5];

% plot label positions
labelPositionsX = {[0.6,2.6],[2.75],[1.3,2.7,2.75],[0.75,2.5], ...
                   [],[],[],[], ...
                   [2.2],[2.5],[2.1],[], ...
                   [],[],[],[], ...
                   [2.25,2.4,3.3],[2.4],[2.2,3.1,3.3],[]};
labelPositionsY = {[17.5,15.5],[15.5],[17.5,16.5,7],[16,15.5], ...
                   [],[],[],[], ...
                   [10.5],[11.5],[5.5],[], ...
                   [],[],[],[], ...
                   [12.5,26.5,27.5],[7.5],[11.5,4.5,29],[]};

%% Variables to modify

% global selector (ONLY MODIFY THIS ONE FROM 1 TO 20)
globalSelector = 2;

% select time window and regression coefficient
coefficientSelector = coefficientSelectors(globalSelector);
windowSelector = windowSelectors(globalSelector);
timeWindow = allTimeWindows(coefficientSelectors(globalSelector)-1,:);

% select position for cluster labels
plotLabelsX = labelPositionsX{globalSelector};
plotLabelsY = labelPositionsY{globalSelector};

%% Load data, then extract trials
for m = freqBaseline
    % averages matrix to run t-test on
    subjectAverages = [];
    
    for isub=1:nsubs
    
        % load sub psd file
        filename = sprintf([path_regression ... 
            's%d/regression_b%d_w%d.mat'],sublist(isub),m,windowSelector);
        if ~exist(filename,'file')
            error(sprintf('%s does not exist!',filename));
        end
        load(filename);
                
        subjectAverages(isub,:,:,:) = regressionData.coeffs;
    end
    
    %% regular 2d cluster correction

    % calculate cluster correction for all image
    selectedTime = (regressionData.t>=timeWindow(1)&...
        regressionData.t<=timeWindow(2));
    dataToTest = subjectAverages(:,coefficientSelector,:,selectedTime);
    [p,praw] = lib.fromToolbox.stat.ClusterCorrection2(dataToTest,ccperms,ccpvalue);
    relevantP = ~isnan(squeeze(p)) & squeeze(p)<0.05;
    columnsToAdd = zeros(size(relevantP,1), ...
        length(regressionData.t)-sum(selectedTime));
    relevantP = [columnsToAdd relevantP];
    p = [nan(size(columnsToAdd)) squeeze(p)];
    
    % create figure to display
    myFig = figure;
    ft_plot_matrix(regressionData.t,regressionData.f, ...
        squeeze(mean(subjectAverages(:,coefficientSelector,:,:),1)), ...
        'clim',[minValue maxValue],'tag','cip','highlightstyle', ...
        'outline','highlight',relevantP);
    colormap(myColormap);
    hold on; axis xy; colorbar;
%     title(figureTitles{coefficientSelector-1});

    clusterP = [];
    clusterLabels = zeros(size(relevantP));

    labels = bwlabel(relevantP);
    if (sum(labels(:))~=0)
        allLabels = unique(labels(:));
        allLabels = allLabels(2:end)';
        for j = allLabels
            testCluster = labels==j;
            testP = mean(p(logical(testCluster)));
            if (testP<ccpvalue)
                clusterP = [clusterP testP];
                clusterLabels(logical(testCluster)) = j;
            end
        end
    end

    plot(0*[1 1],[regressionData.f(1) regressionData.f(end)],'--', ...
        'LineWidth',2,'Color','#888888');
    plot(1.2*[1 1],[regressionData.f(1) regressionData.f(end)],'--', ...
        'LineWidth',2,'Color','#888888');
    plot(1.8*[1 1],[regressionData.f(1) regressionData.f(end)],'--', ...
        'LineWidth',2,'Color','#888888');
    
    for j = 1:length(clusterP)
        text(plotLabelsX(j),plotLabelsY(j), ...
            lib.text.getPValueText(clusterP(j)),'HorizontalAlignment',...
            'center','FontSize',18);
    end
    
    % check for folder existence
    if (exist(path_figures,'file')==0)
        mkdir(path_figures);
    end
    
    % save the figure
    channelValue = mod(globalSelector,channelGroups);
    if (channelValue==0)
        channelValue = channelGroups;
    end
    figureName = figureTitles{coefficientSelector-1};
    figureName = figureName(1:(regexp(figureName,' ','once')-1));
    figureName = sprintf('%s_%s',figureName, ...
        channelGroupLabels{channelValue});
    figurePath = [path_figures figureName];
    print(gcf,'-dpng',figurePath);
    print(gcf,'-dpdf','-painters',figurePath);
end
