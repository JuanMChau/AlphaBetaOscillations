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

%set number of iterations for correlation between reward prospect and behavioural performance
npermcorr = 10000;

%set behavioural variable to correlation. 1=RT difference (high-low reward). 2=accuracy difference (high-low reward).
behmeas = 1;

%set paths
mydir  = pwd;
folders   = strfind(mydir,folderSep);
path = mydir(1:(folders(end-1)-1));
fs = filesep;
datapath  = 'dropTheDataHere/data_neural/processed';
datapath_behav = 'dropTheDataHere/data_behavioural/experiment_matrix_format';
datapath_triggers = 'dropTheDataHere/data_neural/triggers';

%list subjects. Exclude subjects rejecting during processing
sublist = [01 02 03 04 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 31 32 33];
nsubs   = length(sublist);

% custom colormap for graphics
myColormap = lib.graphics.brewermap.brewermap([],'RdBu');
myColormap = myColormap(end:-1:1,:,:);
close all;

% load participant demographics
participantData = readtable('dropTheDataHere/data_behavioural/experiment_processed/AccCounts_all.csv');
participantAge = participantData.Age;
participantSex = participantData.Sex;

%% perform within-subject RSA analyses
for isub=1:nsubs
    
    %load sub EEG file
    filename = sprintf('%s/TaskRep_%02d_EEG_trials.mat',datapath,sublist(isub));
    if ~exist(filename,'file')
        error(sprintf('%s does not exist!',filename));
    end
    load(filename);
    
    %load sub behavioural file
    filename_behav = sprintf('%s/%01d_behav.mat',datapath_behav,sublist(isub));
    if ~exist(filename_behav,'file')
        error(sprintf('%s does not exist!',filename_behav));
    end
    load(filename_behav);
    
    %load list of rejected EEG trials
    filename_rej = sprintf('%s/TaskRep_%02d_rejectedTrials.mat',datapath,sublist(isub));
    if ~exist(filename_rej,'file')
        error(sprintf('%s does not exist!',filename_rej));
    end
    load(filename_rej);
    
    %load trial triggers
    filename_trig = sprintf('%s/ConditionLabels_S%02d',datapath_triggers,isub);
    if ~exist(filename_rej,'file')
        error(sprintf('%s does not exist!',filename_trig));
    end
    load(filename_trig);
    
    %take non-rejected EEG trials and correct trials only
    GoodTrials = 1-rejectedTrialz;
    bindx = data_matrix(:,7) == 1 & GoodTrials';
    eindx = bindx(logical(GoodTrials')) == 1; % eeg index
    eegdat = zscore(eeg.data(:,:,eindx),[],3);
    ntrials = size(eegdat,3);
    timelist = eeg.timepoints./1000;
    
    trigidx = [];
    trigidx = conditionmatrix(eindx, :);
    
    behav = [];
    behav(:,:) = data_matrix(bindx,:);
    
    dat = [];
    dat = eegdat;
    dat = permute(dat,[3 1 2]);
    
    %% baseline data 200-50ms before reward cue presentation
    fsample = 250;
    % the next step allows re-epoching by changing the value in elemtimes:
    elemtimes = [0.0];
    nelems = size(elemtimes,2);
    if size(elemtimes,1) == 1, elemtimes = repmat(elemtimes,[ntrials 1]); end
    
    % baselining options
    base_early = 1;
    
    %time window for analysis (relative to elemtimes!!)
    dtime = fix([-1.000 +5.000]*fsample);
    %time window for baselining (relative to elemtimes!!)
    % dbase = fix([-0.200 -0.050]*fsample);
    dbase = fix([-0.600 -0.300]*fsample);
    dstep = 1;
    %create new time vector based on selected window
    time = (dtime(1):dstep:dtime(2))/fsample;
    ntimes = length(time);
    
    %channels to analyse
    chanlist = 1:61;
    nchans   = length(chanlist);
    
    %preallocate new data array
    dataerp_elem = nan(ntrials,nchans,ntimes,nelems);
    %re-epoch data and baseline
    for itrial = 1:ntrials
        for ielem = 1:nelems
            ionset = find(timelist <= elemtimes(itrial,ielem),1,'last');
            itime = ionset+(dtime(1):dstep:dtime(2));
            itime = min(itime,length(timelist));
            dataerp_elem(itrial,:,:,ielem) = dat(itrial,chanlist,itime);
            if base_early
                izero = find(timelist <= 0,1,'last');
                ibase = izero + (dbase(1):dstep:dbase(2));
                dataerp_elem(itrial,:,:,ielem) = squeeze(dataerp_elem(itrial,:,:,ielem))-repmat(permute(mean(dat(itrial,chanlist,ibase),3),[2 1]),[1,ntimes]);
            end
        end
    end
    
    if base_early
        dat = dataerp_elem;
    end
    
    %% group trials into relevant conditions
    for i=1:size(dat,1)
        switch trigidx(i,1)
            %low reward
            case 5
                
                %identify colour trials versus shape trials
                %stim order: yellow square, blue square, yelllow circle, blue
                %circle
                %let 1=yellow, 2=blue, 3 = square, 4 = circle
                switch trigidx(i,2)
                    % colour task cue
                    case {7 8}
                        switch trigidx(i,3)
                            %yellow square
                            case 11
                                %order name as rel vs irrel feat.
                                LoyelSqr(i,1)   = 1;
                                %blue square
                            case 12
                                LobluSqr(i,1) = 1;
                                %yellow circle
                            case 13
                                LoyelCir(i,1) = 1;
                                %blue circle
                            case 14
                                LobluCir(i,1) = 1;
                        end
                        
                        %shape task cue
                    case {9 10}
                        switch trigidx(i,3)
                            % yellow square
                            case 11
                                LosqrYel(i,1) = 1;
                                % blue square
                            case 12
                                LosqrBlu(i,1) = 1;
                                % yellow circle
                            case 13
                                LocirYel(i,1) = 1;
                                % blue circle
                            case 14
                                LocirBlu(i,1) = 1;
                        end
                end
                
                % high reward
            case 6
                switch trigidx(i,2)
                    % colour task cue
                    case {7 8}
                        switch trigidx(i,3)
                            %yellow square
                            case 11
                                %order name as rel vs irrel feat.
                                HiyelSqr(i,1)   = 1;
                                %blue square
                            case 12
                                HibluSqr(i,1) = 1;
                                %yellow circle
                            case 13
                                HiyelCir(i,1) = 1;
                                %blue circle
                            case 14
                                HibluCir(i,1) = 1;
                        end
                        
                        %shape task cue
                    case {9 10}
                        switch trigidx(i,3)
                            % yellow square
                            case 11
                                HisqrYel(i,1) = 1;
                                % blue square
                            case 12
                                HisqrBlu(i,1) = 1;
                                % yellow circle
                            case 13
                                HicirYel(i,1) = 1;
                                % blue circle
                            case 14
                                HicirBlu(i,1) = 1;
                        end
                end
        end
    end
    
    LoyelSqrtrls = find(LoyelSqr==1);
    LobluSqrtrls = find(LobluSqr==1);
    LoyelCirtrls = find(LoyelCir==1);
    LobluCirtrls = find(LobluCir==1);
    LosqrYeltrls = find(LosqrYel==1);
    LosqrBlutrls = find(LosqrBlu==1);
    LocirYeltrls = find(LocirYel==1);
    LocirBlutrls = find(LocirBlu==1);
    
    HiyelSqrtrls = find(HiyelSqr==1);
    HibluSqrtrls = find(HibluSqr==1);
    HiyelCirtrls = find(HiyelCir==1);
    HibluCirtrls = find(HibluCir==1);
    HisqrYeltrls = find(HisqrYel==1);
    HisqrBlutrls = find(HisqrBlu==1);
    HicirYeltrls = find(HicirYel==1);
    HicirBlutrls = find(HicirBlu==1);
    
    LoyelSqrk(:, :, :) = dat(LoyelSqrtrls,:,:);
    LobluSqrk(:, :, :) = dat(LobluSqrtrls,:,:);
    LoyelCirk(:, :, :) = dat(LoyelCirtrls,:,:);
    LobluCirk(:, :, :) = dat(LobluCirtrls,:,:);
    LosqrYelk(:, :, :) = dat(LosqrYeltrls,:,:);
    LosqrBluk(:, :, :) = dat(LosqrBlutrls,:,:);
    LocirYelk(:, :, :) = dat(LocirYeltrls,:,:);
    LocirBluk(:, :, :) = dat(LocirBlutrls,:,:);
    
    HiyelSqrk(:, :, :) = dat(HiyelSqrtrls,:,:);
    HibluSqrk(:, :, :) = dat(HibluSqrtrls,:,:);
    HiyelCirk(:, :, :) = dat(HiyelCirtrls,:,:);
    HibluCirk(:, :, :) = dat(HibluCirtrls,:,:);
    HisqrYelk(:, :, :) = dat(HisqrYeltrls,:,:);
    HisqrBluk(:, :, :) = dat(HisqrBlutrls,:,:);
    HicirYelk(:, :, :) = dat(HicirYeltrls,:,:);
    HicirBluk(:, :, :) = dat(HicirBlutrls,:,:);
    
    
    fprintf('\nCalculating RDM for subject %d.\n',sublist(isub));
    
    %average over trials for each condition and store in a single matrix
    condsxchansxtime(1, :, :)  = squeeze(mean(LoyelSqrk(:,:,:), 1));
    condsxchansxtime(2, :, :)  = squeeze(mean(LobluSqrk(:,:,:), 1));
    condsxchansxtime(3, :, :)  = squeeze(mean(LoyelCirk(:,:,:), 1));
    condsxchansxtime(4, :, :)  = squeeze(mean(LobluCirk(:,:,:), 1));
    condsxchansxtime(5, :, :)  = squeeze(mean(LosqrYelk(:,:,:), 1));
    condsxchansxtime(6, :, :)  = squeeze(mean(LosqrBluk(:,:,:), 1));
    condsxchansxtime(7, :, :)  = squeeze(mean(LocirYelk(:,:,:), 1));
    condsxchansxtime(8, :, :)  = squeeze(mean(LocirBluk(:,:,:), 1));
    condsxchansxtime(9, :, :)  = squeeze(mean(HiyelSqrk(:,:,:), 1));
    condsxchansxtime(10, :, :) = squeeze(mean(HibluSqrk(:,:,:), 1));
    condsxchansxtime(11, :, :) = squeeze(mean(HiyelCirk(:,:,:), 1));
    condsxchansxtime(12, :, :) = squeeze(mean(HibluCirk(:,:,:), 1));
    condsxchansxtime(13, :, :) = squeeze(mean(HisqrYelk(:,:,:), 1));
    condsxchansxtime(14, :, :) = squeeze(mean(HisqrBluk(:,:,:), 1));
    condsxchansxtime(15, :, :) = squeeze(mean(HicirYelk(:,:,:), 1));
    condsxchansxtime(16, :, :) = squeeze(mean(HicirBluk(:,:,:), 1));
    timexchansxconds = permute(condsxchansxtime, [3 2 1]);
    
    %% create model representational dissimilarity matrices (RDMs)
    %reward coding model
    modeldismatrixrwd = squareform([ones(8)*0 ones(8)*1; ones(8)*1 ones(8)*0]);
    
    %task rule coding model
    modeldismatrixtsk= ([ones(4)*0 ones(4)*1 ones(4)*0 ones(4)*1; ones(4)*1 ones(4)*0 ones(4)*1 ones(4)*0; ...
        ones(4)*0 ones(4)*1 ones(4)*0 ones(4)*1; ones(4)*1 ones(4)*0 ones(4)*1 ones(4)*0]);
    modeldismatrixtsk=squareform(modeldismatrixtsk);
    
    %relevant feat coding model
    f1 = [0 1 0 1 1 1 1 1]; f2 = [1 0 1 0 1 1 1 1]; f3 = [1 1 1 1 0 0 1 1]; f4 = [1 1 1 1 1 1 0 0];
    modeldismatrixrelfeat = [f1 f1; f2 f2; f1 f1; f2 f2; f3 f3; f3 f3; f4 f4; f4 f4; ...
        f1 f1; f2 f2; f1 f1; f2 f2; f3 f3; f3 f3; f4 f4; f4 f4];
    modeldismatrixrelfeat = squareform(modeldismatrixrelfeat);
    
    %irrelevant feat coding model
    ir1=[0 0 1 1 1 1 1 1]; ir2 = [1 1 0 0 1 1 1 1]; ir3 = [1 1 1 1 0 1 0 1]; ir4 = [1 1 1 1 1 0 1 0];
    modeldismatrixirrelfeat = [ir1;ir1;ir2;ir2;ir3;ir4;ir3;ir4];
    modeldismatrixirrelfeat=repmat(modeldismatrixirrelfeat,2);
    modeldismatrixirrelfeat=squareform(modeldismatrixirrelfeat);
    
    %motor coding model
    m1 = [0 1 0 1 0 0 1 1 0 1 0 1 0 0 1 1]; m2 = [1 0 1 0 1 1 0 0 1 0 1 0 1 1 0 0];
    modeldismatrixhand = [m1;m2;m1;m2;m1;m1;m2;m2;m1;m2;m1;m2;m1;m1;m2;m2];
    modeldismatrixhand = squareform(modeldismatrixhand);
    
    X = [modeldismatrixrwd' modeldismatrixtsk' modeldismatrixrelfeat' modeldismatrixirrelfeat' modeldismatrixhand'];
    X = zscore(X);
    
    
    %% perform model fitting with model and data RDMs
    
    for m=1:length(timexchansxconds) %m = time
        make2D = squeeze(timexchansxconds(m,:,:)); %chans by conds
        stimDistances = pdist(make2D', 'mahalanobis', lib.fromToolbox.decoding.covdiag(dat(:,:,m)));
        
        %compute unzscored version for illustrative RDMs
        stimDistances_nozscore = stimDistances;
        stimdismatrix_nozscore = squareform(stimDistances_nozscore);
        RDM_output_nozscore(isub,:,:,m) = stimdismatrix_nozscore(:,:,:);
        
        %zscore condition distances to make coding models comparable
        stimDistances = zscore(stimDistances);
        stimdismatrix = squareform(stimDistances);
        RDM_output(isub,:,:,m) = stimdismatrix(:,:,:);
        
        b = glmfit(X,stimDistances','normal','constant','on');
        RDM_fit_glm(isub,m,:) = b;
        
    end
    
    clear('chansxtimexconds','condsxchansxtime', 'timexchansxconds');
    clear -regexp ^Lo ^Hi;
end

%% Make a quick plot

% data settings
plotsubs = [1:30]; %subjects to plot
plotconditions = 2; %conditions to plot

%select data for plotting
%datamatrix should be participants X time points X conditions (or
%regressors, etc.)
plottimewindow = [-0.200 +3.500];
itime = timelist>=plottimewindow(1)&timelist<=plottimewindow(2);
plottime = timelist(itime);
plotdata = RDM_fit_glm(:,itime,plotconditions);

% plotting settings
xlims = plottime([1 end]); %set horizontal axis limits
ylims = [-0.2 0.5]; %set vertical axis limits
eventonsets = [0]; %onsets of events (relative to plottime!)

% condition labels
condnames = {'Reward Coding'};
% permutation test settings
corrwin = [+0.000 +3.500]; %time window in which to test (against 0)
nperm = 10000; %in final analysis for a paper this should be 10000 or more
plottestnames = {'Ttest'};
pthresh = 0.05; %threshold at which p-values show up on the graph

% output settings
analysisname = 'Fig_3C_reward';
save_stats = true;
do_print = true;

%smooth data?
fsample = round(median(1./diff(timelist)));
ftype = 'gaussian';
fsize = 0.012*fsample;
plotdata = lib.fromToolbox.common.filtfast(plotdata,2,[],ftype,fsize);

% open figure
clc
figure
% color settings
plotcol = [0 0 0];
lib.fromToolbox.plt.SetupPlotting();

%add vertical and horizontal reference lines
hold on
plot(xlims,[0 0],'k-','linewidth',1.0,'color',[1 1 1]*0.25); %plot horizontal line
for ievent = 1:length(eventonsets)
    %plot vertical lines at times of interest
    plot(eventonsets(ievent)*ones(1,2),ylims,'k-','linewidth',1.0,'color',[1 1 1]*0.25)
end
set(gca,'box','off')

% plot sig bars and perform cluster-based permutation test
it = plottime>=corrwin(1) & plottime<= corrwin(2);
sigtime = plottime(it);
stats = struct;
stats.pthresh = pthresh;
stats.corrwin = corrwin;
stats.nperm = nperm;
stats.sigtime = sigtime;
stats.pvals = [];

sigcol = plotcol;
testdata = plotdata;

nsig = size(testdata,3);
for isig = 1:nsig
    [~,praw] = ttest(testdata(:,it,isig));
    [datobs, datrnd] = lib.fromToolbox.stat.cluster_test_helper(testdata(:,it,isig)', nperm,'t');
    [h, p, clusterinfo] = lib.fromToolbox.stat.cluster_test(datobs, datrnd, 0, 0.05,0.05, 'sum');
    pcorr = p';
    barcol = [0 0 0];
    yheight = min(ylims) + range(ylims)*0.05*(isig);
    pheight = min(ylims) + range(ylims)*0.05*(isig+0.5)+0.02;
    b = bwconncomp(pcorr(:,:)<pthresh);
    stepsize = median(diff(sigtime))/4;
    stats.pvals(isig).praw = praw;
    stats.pvals(isig).pcorr = pcorr;
    stats.pvals(isig).name = plottestnames{1};
    for ib = 1:b.NumObjects
        currid = b.PixelIdxList{ib};
        if nnz(currid)>0
            xvec = [sigtime(currid(1))-stepsize sigtime(currid(end))+stepsize];
            plot(xvec,ones(length(xvec),1)*yheight,'k-','linewidth',4,'color',barcol);
            addSignificanceText(pcorr(1,currid(1)),[mean(xvec),pheight]);
        end
    end
end
clear testdata;

% plot time courses (with shaded error bars)
cfg = [];
hplot = lib.fromToolbox.plt.plotpatch(plotdata,plottime,plotcol,cfg);
legend(hplot,condnames,'fontsize',16), legend boxoff
xlim(xlims); ylim(ylims);

% set axis ticks, axis labels, etc.
xstpsize = 1;
set(gca,'xtick',[0:xstpsize:4]);
set(gca,'ytick',unique([0 ylims]));

figpath =  'figures/erp_psd_correlation'; if ~exist(figpath,'dir'), mkdir(figpath); end
fname = sprintf('%s/%s',figpath,'reward_erp');
print(gcf,'-dpng',fname);
print(gcf,'-dpdf','-painters',fname);

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

%% Correlations
bothFreqs = [8 12; 20 30];
avgWindow = [0 1];
freqNames = {'alpha','beta'};
roiNames = {'all','frontal','central','posterior'};

for j = 1:length(channelsForFigures)
    decodingTime = timelist>=avgWindow(1) & timelist<=(avgWindow(2)+0.1);
    rewardDecoding = RDM_fit_glm(:,decodingTime,2);
    meanDecoding = squeeze(mean(rewardDecoding,2));

    varsForCorr = [medianDiffRT' meanDecoding];

    for i = 1:size(bothFreqs,1)
        psdTime = (lowReward.t>=avgWindow(1)) & (lowReward.t<=avgWindow(2));
        psdFreq = (lowReward.f>=bothFreqs(i,1)) & (lowReward.f<=bothFreqs(i,2));

        [~,selectedChannels,~] = intersect(loadedLayout.layout, ...
            channelsForFigures{j});
        dataForCorr = subjectAvgDiff(:,selectedChannels,psdFreq,psdTime);
        varsForCorr = [varsForCorr mean(mean(mean(dataForCorr,4),3),2)];
    end
    
    disp(roiNames{j});
    [r,p] = partialcorr([varsForCorr participantAge participantSex],'Type','Spearman')
    
    for m = 1:size(r,1)
        for n = 1:size(r,2)
            if (m<n)
                r(m,n) = 0;
            end
        end
    end
    
    figure; hold on;    
    imagesc(r(1:4,1:4),[-1 1]);
    colormap(myColormap);
    colorbar;

    for m = 1:size(r,1)
        for n = 1:size(r,2)
            if (m<n)
                pvalCharacter = getSignificantCharacter(p(m,n));
                if (isempty(pvalCharacter))
                    text(m-0.35,n,sprintf('%1.3f%s',p(m,n),pvalCharacter), ...
                        "FontSize",16);
                else
                    text(m-0.35,n,sprintf('%1.3f%s',p(m,n),pvalCharacter), ...
                        "FontWeight","bold","FontSize",16);
                end
            end
        end
    end
    
    xticks(1:4);
    xticklabels({'RT_{diff}','ERP','Alpha','Beta'});
    yticks(1:4);
    yticklabels({'RT_{diff}','ERP','Alpha','Beta'});
    set(gca,'TickLength',[0 0]);
    ax = gca;
    ax2 = axes('Position',ax.Position,'XColor',[1 1 1],...
        'YColor',[1 1 1],'Color','none','XTick',[],'YTick',[]);
    set(gca, 'YDir','reverse');
    
    fname = sprintf('%s/%s_%s',figpath,'correlation_all_mat',roiNames{j});
    print(gcf,'-dpng',fname);
    print(gcf,'-dpdf','-painters',fname);

    close all;
end

%% Additional functions

function significanceText = getSignificantCharacter(pvalue)
    significanceText = '';
    
    if (pvalue<0.05)
        significanceText = '*';
    end
end

function addSignificanceText(pvalue,textPosition)
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
        return;
    end

    text(textPosition(1),textPosition(2),significanceText, ...
        "FontSize",18,"HorizontalAlignment","center");
end