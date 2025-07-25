% Author: Juan M. Chau
% Last uploaded: 01.11.2022

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

%% Initial variables

% apply log transform to psd
logTransform = 1;

% set value of baseline correction used (0, 1, 2, 3)
% 0: no baseline
% 1: percent of variance -> (data-mean)/mean
% 2: zscore -> (data-mean)/std
% 3: center around 0 -> data-mean
freqBaseline = 3;

% set number of subsample iterations to match number of high and low reward
% trials
numSubSamp = 1;

% list subjects. Exclude subjects rejecting during processing
sublist = [1:4 6:27 30:33];
% sublist = [01 02 03 04 06 07];
% sublist = 1;
nsubs = length(sublist);

% load channel layout for later
loadedLayout = load('config/layout.mat');

% extracted time window for baseline
% timeWindowBaseline = [-0.200 -0.050; 1.000 1.150; 1.600 1.750; 1.600 1.750];
timeWindowBaseline = [-0.600 -0.300;
                      0.950 1.150;
                      1.550 1.750];

% time window to consider for regression
% fullTimeWindow = [-0.200 4.500; 1.000 4.500; 1.600 4.500; 1.600 4.500];
fullTimeWindow = [-0.600 4.500;
                  0.950 4.500;
                  1.550 4.500];
              
% groups of channels to consider for regression
channelGroups = {loadedLayout.layout, ...
    {'AF3','AF4','AFz','F1','F2','F3','F4','Fz'}, ...
    {'C1','C2','C3','C4','Cz'}, ...
    {'O1','O2','Oz','PO3','PO4','POz'}};

% savepath spectrogram
path_spectrogram_config = 'spectrograms/config/';
path_spectrogram_psd = 'spectrograms/psd/';

% savepath regression
path_regression_data = 'regression/data_rsa_all/';

% build the selector arrays
windowSelector = [];
channelSelector = [];
for i = 1:length(fullTimeWindow)
    for j = 1:length(channelGroups)
        % baseline and time window selectors
        windowSelector = [windowSelector i];
        % channel selector
        channelSelector = [channelSelector j];
    end    
end

%% Variables to modify

% TODO: ONLY MODIFY THIS ONE FROM 1 TO 12 (global selector)
globalSelector = 12;

% apply all selections
baselineWindow = timeWindowBaseline(windowSelector(globalSelector),:);
timeWindow = fullTimeWindow(windowSelector(globalSelector),:);
[~,selectedChannels,~] = intersect(loadedLayout.layout, ...
    channelGroups{channelSelector(globalSelector)});
selectedChannels = sort(selectedChannels);

% selected triggers
selectedTriggers = 1:16;

%% Create model representational dissimilarity matrices (RDMs)

% base RDM matrix (2 by 2)
baseRDM = [0 1; 1 0];

% reward coding matrix
X1 = squareform(round(imresize(baseRDM,8)))';

% task coding matrix
X2 = squareform(repmat(round(imresize(baseRDM,4)),2,2))';

% relevant/irrelevant features coding matrix
corner1 = repmat(baseRDM,2,2);
corner2 = round(imresize(baseRDM,2));
X3 = squareform(repmat([corner1 ones(4); ones(4) corner2],2,2))'; % relevant
X4 = squareform(repmat([corner2 ones(4); ones(4) corner1],2,2))'; % irrelevant

% motor coding matrix
corner3 = round(imresize(repmat(baseRDM,1,2),[4 4]));
X5 = squareform(repmat([corner1 corner3'; corner3 corner2],2,2))';

X = zscore([X1 X2 X3 X4 X5]);

%% Load data, then extract trials
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
        
    % Baselining process
    for i = 1:subjectData.trialDivision:subjectData.length
        % check for last value
        if ((i+subjectData.trialDivision-1)>=subjectData.length)
            selectedValues = i:subjectData.length;
        else
            selectedValues = i:(i+subjectData.trialDivision-1);
        end
        
        % load sub psd file
        filename = sprintf('%s/s%d/laplacian/s%d_%d_%d.mat', ...
            path_spectrogram_psd,sublist(isub),sublist(isub), ...
            selectedValues(1),selectedValues(end));
        if ~exist(filename,'file')
            error(sprintf('%s does not exist!',filename));
        end
        load(filename);
        
        % select desired channels
        tfDataPSD.psd = tfDataPSD.psd(:,selectedChannels,:,:);
                
        % select time window
        timeIndices = tfDataPSD.t>=timeWindow(1) & ...
            tfDataPSD.t<=timeWindow(2);
        tfDataPSD.psd = tfDataPSD.psd(:,:,:,timeIndices);
        tfDataPSD.t = tfDataPSD.t(timeIndices);
        
        % reconstruct all data
        for m = freqBaseline
            eval(['varExists = ~exist(''bigPSD' num2str(m) ...
                ''',''var'');']);
            if (varExists)
                eval(['bigPSD' num2str(m) ...
                    ' = zeros([subjectData.length' ...
                    ' size(tfDataPSD.psd,2:4)]);']);
            end
        end
                
        fprintf('\nBaselining trials for subject %d: %d-%d ', ...
            sublist(isub),selectedValues(1),selectedValues(end));
        tic;
        % now fill big trial variables
        for j = selectedValues            
            % correct trial indices
            myindex = mod(j-1,length(selectedValues))+1;
            
            for m = freqBaseline                
                % perform s baseline correction
%                 s = squeeze(tfDataS.s(myindex,:,:,:));
%                 [s,t] = lib.reusable.baselineCorrection(s,tfDataS.t, ...
%                     0,timeWindowBaseline,logTransform,baselineType);
                
                % perform psd baseline correction
                psd = squeeze(tfDataPSD.psd(myindex,:,:,:));
                [psd,t] = lib.reusable.baselineCorrection(psd, ...
                    tfDataPSD.t,0,baselineWindow,logTransform,m);
                
                eval(['bigPSD' num2str(m) '(j,:,:,:) = psd;']);
            end
        end
        fprintf('took %f seconds',toc);
    end
        
    % Process by baseline method
    for m = freqBaseline 
        % select the PSD matrix for that baseline method
        eval(['currentPSD = bigPSD' num2str(m) ';']);
        
        % calculate the original covariance matrix
        covPSD = zeros(size(currentPSD,2),size(currentPSD,2), ...
            length(tfDataPSD.f),length(tfDataPSD.t));
        for f = 1:length(tfDataPSD.f) % frequency
            for t = 1:length(tfDataPSD.t) % time
                covPSD(:,:,f,t) = lib.fromToolbox.decoding.covdiag(currentPSD(:,:,f,t));
            end
        end
        
        % calculate distances
        fprintf('\nProcessing data');

        % group permutation by condition
        conditionPSD = zeros([length(unique(conditionIndex)) ...
            size(currentPSD,2:4)]);
        for p = unique(conditionIndex)'
            indices = (conditionIndex==p);
            conditionPSD(p,:,:,:) = ...
                squeeze(mean(currentPSD(indices,:,:,:),1));
        end
        
        % remove unwanted triggers
        conditionPSD = conditionPSD(selectedTriggers,:,:,:);

        % just to make dealing with the dimensions easier        
        conditionPSD = permute(conditionPSD,[2 1 3 4]);            

        % create outputs over all conditions
        RDMSubSampIn = zeros(size(X,1),length(tfDataPSD.f), ...
            length(tfDataPSD.t));

        % fill this matrix
        for f = 1:length(tfDataPSD.f) % frequency
            for t = 1:length(tfDataPSD.t) % time
                make2D = conditionPSD(:,:,f,t)';
                stimDistances = pdist(make2D,'mahalanobis', ...
                    covPSD(:,:,f,t));
                stimDistances = zscore(stimDistances);
                RDMSubSampIn(:,f,t) = stimDistances;
            end
        end

        % perform linear regression
        RDMSubSampOut = lib.fromToolbox.stat.massGLM(X,RDMSubSampIn,1);
        
        % check for folder existence
        if (exist([path_regression_data 's' ...
            num2str(sublist(isub))],'file')==0)
            mkdir(path_regression_data,['s' num2str(sublist(isub))]);
        end
        % save the regression file
        regressionData = struct('coeffs',RDMSubSampOut, ...
            't',tfDataPSD.t,'f',tfDataPSD.f);
        filename = sprintf([path_regression_data ...
            's%d/regression_b%d_w%d.mat'],sublist(isub),m,globalSelector);
        save(filename,'regressionData');
                
        clear covPSD currentPSD RDMSubSampOut;
    end
    
    clear -regexp ^bigPSD;
end
