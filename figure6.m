% Author: Juan M. Chau
% Last uploaded: 01.22.2024

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

%% setup

% savepath figures
path_figures = 'figures/rankCorrelation/';

% cluster correction permutations
ccperms = 1e4;

% cluster correction significance level
ccpvalue = 0.05;
ccalpha = 0.025;

% custom colormap for graphics
myColormap = lib.graphics.brewermap.brewermap([],'RdBu');
myColormap = myColormap(end:-1:1,:,:);
close all;

%list subjects. Exclude subjects rejecting during processing
sublist = [01 02 03 04 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 31 32 33];
nsubs   = length(sublist);

figureTitles = {'All channels','Frontal','Central','Posterior'};

%% Read PSD data

% load channel layout for later
loadedLayout = load('config/layout.mat');

% channels used for plotting
channelsForFigures = {loadedLayout.layout, ...
    {'AF3','AF4','AFz','F1','F2','F3','F4','Fz'}, ...
    {'C1','C2','C3','C4','Cz'}, ...
    {'O1','O2','Oz','PO3','PO4','POz'}};

typeOfPSD = {'laplacian','original'};
psdSelector = 1;
path_spectrogram_psd_avg = ['spectrograms/psd/avg/' ...
    typeOfPSD{psdSelector} '/all/'];
path_figures = [path_figures typeOfPSD{psdSelector} '/'];

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

%% Perform correlation permutations

allClusters = {};
allLabels = {};
allR = {};
allP = {};
allDataForCC = {}; 

for i = 1:length(channelsForFigures)
    [~,selectedChannels,~] = intersect(loadedLayout.layout, ...
        channelsForFigures{i});
    allDataForCC{i} = subjectAvgDiff(:,selectedChannels,:,:);

    % Calculate permutations
    [cluster labels r p] = permutationCorrelation(allDataForCC{i}, ...
        medianDiffRT,ccperms,ccalpha,ccpvalue);

    allClusters{i} = cluster;
    allLabels{i} = labels;
    allR{i} = r;
    allP{i} = p;
end

%% Plot results

% check for folder existence
if (exist(path_figures,'file')==0)
    mkdir(path_figures);
end

zmin = -0.7;
zmax = 0.7;

labelPositionsX = {[1.25,1.25],[1.5],[],[1]};
labelPositionsY = {[16,32],[28],[],[13.5]};

for i = 1:length(channelsForFigures)
    % Display differences
    myFig = figure;
    ft_plot_matrix(highReward.t,highReward.f,allR{i},'clim', ...
        [zmin zmax],'tag','cip','highlightstyle', ...
        'outline','highlight',allClusters{i});
    colormap(myColormap);
    hold on; axis xy; colorbar;

    plot(0*[1 1],[highReward.f(1) highReward.f(end)],'--', ...
        'LineWidth',2,'Color','#888888');
    plot(1.2*[1 1],[highReward.f(1) highReward.f(end)],'--', ...
        'LineWidth',2,'Color','#888888');
    plot(1.8*[1 1],[highReward.f(1) highReward.f(end)],'--', ...
        'LineWidth',2,'Color','#888888');
    
    for j = 1:length(allP{i})
        text(labelPositionsX{i}(j),labelPositionsY{i}(j), ...
            lib.text.getPValueText(allP{i}(j)), ...
            'HorizontalAlignment','center','FontSize',16);
    end    
    
    figureName = figureTitles{i};
    figureName = sprintf('%s_normal',figureName);
    figurePath = [path_figures figureName];
    print(gcf,'-dpng',figurePath);
    print(gcf,'-dpdf','-painters',figurePath);
    close(myFig);

end

%% Functions (testing atm)

function [baseCluster,clusterLabels,baseR,clusterP] = ...
        permutationCorrelation(tfData,RTs,nPerm,alphaThreshold,...
        clusterThreshold)
    tfData = squeeze(mean(tfData,2));    
    allClusterMass = zeros(1,nPerm-1);
    
    parfor i = 1:(nPerm-1)
        if (mod(i-1,10)==0)
            disp(['Handling permutation ' num2str(i)]);
        end
        shuffleValues = randperm(length(RTs));
        [~,allClusterMass(i)] = getCluster(tfData, ...
            RTs(shuffleValues),alphaThreshold);
    end
    
    baseR = zeros(size(tfData,2:3));
    baseP = zeros(size(tfData,2:3));

    for i = 1:size(tfData,2)
        for j = 1:size(tfData,3)
            [baseR(i,j),baseP(i,j)] = corr(tfData(:,i,j),RTs', ...
                'Type','Spearman');
        end
    end

    clusterP = [];

    baseCluster = zeros(size(tfData,2:3));
    clusterLabels = zeros(size(tfData,2:3));
    labels = bwlabel(baseP<clusterThreshold);
    if (sum(labels(:))~=0)
        allLabels = unique(labels(:));
        allLabels = allLabels(2:end)';
        for i = allLabels
            testCluster = labels==i;
            testMass = sum(abs(baseR(logical(testCluster))));
            testP = (1+sum(allClusterMass>testMass))/nPerm;
            if (testP<alphaThreshold)
                baseCluster(logical(testCluster)) = 1;
                clusterP = [clusterP testP];
                clusterLabels(logical(testCluster)) = i;
            end
        end
    end
end

function [cluster,clusterMass] = getCluster(tfData,RTs,alphaThreshold)
    rMatrix = zeros(size(tfData,2:3));
    pMatrix = zeros(size(tfData,2:3));
    cluster = zeros(size(tfData,2:3));
    clusterMass = 0;

    for i = 1:size(tfData,2)
        for j = 1:size(tfData,3)
            [rMatrix(i,j),pMatrix(i,j)] = corr(tfData(:,i,j),RTs', ...
                'Type','Spearman');
        end
    end
    labels = pMatrix<alphaThreshold;
    CC = bwconncomp(labels);
    if CC.NumObjects ~= 0
        RP = regionprops(CC);
        [~,maxIndex] = max([RP.Area]);
        cluster(CC.PixelIdxList{maxIndex}) = 1;
        clusterMass = sum(abs(rMatrix(logical(cluster))));
    end
end

function addClusterText(clusters,xrange,yrange,pvalues)
    stats = regionprops(clusters,'BoundingBox');
    [rows,cols] = size(clusters);

    clusterNames = unique(clusters(clusters(:)~=0));
    
    for i = 1:length(clusterNames)
        bbox = stats(clusterNames(i)).BoundingBox;
        yTop = floor(bbox(2));
        yBottom = floor(bbox(2)+bbox(4));
        
        if(yBottom>(rows/2))
            textPosition = [floor(bbox(1)+bbox(3)/2),yTop-1];
        else
            textPosition = [floor(bbox(1)+bbox(3)/2),yBottom+1];
        end
        
        text(xrange(textPosition(1)),yrange(textPosition(2)),...
            lib.text.getPValueText(pvalues(i)),'HorizontalAlignment',...
            'center');
    end
end