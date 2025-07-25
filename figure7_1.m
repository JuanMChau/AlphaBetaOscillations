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

% create arrays for average ERPs
averageERPs = [];

% load channel layout for later
loadedLayout = load('config/layout.mat');

% channels used for plotting
channelsForFigures = {loadedLayout.layout, ...
    {'AF3','AF4','AFz','F1','F2','F3','F4','Fz'}, ...
    {'C1','C2','C3','C4','Cz'}, ...
    {'O1','O2','Oz','PO3','PO4','POz'}};

% ROI names
roiNames = {'Frontal','Central','Posterior'};

% savepath spectrogram
path_spectrogram_psd_avg = 'spectrograms/psd/avg/laplacian/split/';

%% Load psd data (just to get the split indices)
filename = [path_spectrogram_psd_avg 'allSubjects_low.mat'];
if ~exist(filename,'file')
    error(sprintf('%s does not exist!',filename));
end
load(filename);

avgPSDLow = avgPSD;

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
    eegdat = eeg.data(:,:,eindx);
    ntrials = size(eegdat,3);
    timelist = eeg.timepoints./1000;
    
    trigidx = [];
    trigidx     = conditionmatrix(eindx, :);
    
    behav = [];
    behav(:,:) = data_matrix(bindx,:);
    
    dat = [];
    dat = eegdat;
    dat   = permute(dat,[3 1 2]);
    
    %% baseline data 200-50ms before reward cue presentation
    fsample = 250;
    % the next step allows re-epoching by changing the value in elemtimes:
    elemtimes = [0.0];
    nelems = size(elemtimes,2);
    if size(elemtimes,1) == 1, elemtimes = repmat(elemtimes,[ntrials 1]); end
    
    % baselining options
    base_early  = 1;
    
    %time window for analysis (relative to elemtimes!!)
    dtime = fix([-1.000 +5.000]*fsample);
    %time window for baselining (relative to elemtimes!!)
    % dbase = fix([-0.200 -0.050]*fsample);
    dbase = fix([-0.600 -0.300]*fsample);
    dstep   = 1;
    %create new time vector based on selected window
    time    = (dtime(1):dstep:dtime(2))/fsample;
    ntimes  = length(time);
    
    %channels to analyse
    chanlist = 1:61;
    nchans   = length(chanlist);
    
    %preallocate new data array
    dataerp_elem = nan(ntrials,nchans,ntimes,nelems);
    %re-epoch data and baseline
    for itrial = 1:ntrials
        for ielem = 1:nelems
            ionset = find(timelist <= elemtimes(itrial,ielem),1,'last');
            itime  = ionset+(dtime(1):dstep:dtime(2));
            itime  = min(itime,length(timelist));
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

    for roi = 2:4
        [~,selectedChannels,~] = intersect(loadedLayout.layout, ...
            channelsForFigures{roi});
        lowRewardERPs = squeeze(mean(mean(dat(trigidx(:,1)==5, ...
            selectedChannels,:),1),2));
        highRewardERPs = squeeze(mean(mean(dat(trigidx(:,1)==6, ...
            selectedChannels,:),1),2));
        averageERPs(isub,:,1,roi-1) = lowRewardERPs;
        averageERPs(isub,:,2,roi-1) = highRewardERPs;
    end
end

%% Plot average ERPs as well

% data settings
plotsubs = [1:30]; %subjects to plot
plotconditions = 2; %conditions to plot

%select data for plotting
%datamatrix should be participants X time points X conditions (or
%regressors, etc.)
plottimewindow = [-0.200 +3.500];
itime = timelist>=plottimewindow(1)&timelist<=plottimewindow(2);
plottime = timelist(itime);

% plotting settings
xlims = plottime([1 end]); %set horizontal axis limits
ylims = [-10 12]; %set vertical axis limits
eventonsets = [0]; %onsets of events (relative to plottime!)

% condition labels
condnames = {'Low reward','High reward'};
% permutation test settings
corrwin = [+0.000 +3.500]; %time window in which to test (against 0)
nperm = 10000; %in final analysis for a paper this should be 10000 or more
plottestnames = {'Ttest'};
pthresh = 0.05; %threshold at which p-values show up on the graph

for reward = unique(avgPSDLow.globalMedianRT)
    for roi = 1:size(averageERPs,4)
        % Select data to plot
        plotdata = averageERPs(avgPSDLow.globalMedianRT==reward,itime,:,roi);
    
        %smooth data?
        fsample = round(median(1./diff(timelist)));
        ftype = 'gaussian';
        fsize = 0.012*fsample;
        plotdata = lib.fromToolbox.common.filtfast(plotdata,2,[],ftype,fsize);
    
        % open figure
        clc
        figure
        % color settings
        plotcol = colormap(lib.fromToolbox.plt.linspecer(size(plotdata,3),'qualitative'));
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
        stats   = struct;
        stats.pthresh = pthresh;
        stats.corrwin = corrwin;
        stats.nperm   = nperm;
        stats.sigtime = sigtime;
        stats.pvals   = [];
        
        sigcol   = [0 0 0];
        testdata = diff(plotdata,[],3);
    
        nsig     = size(testdata,3);
        for isig = 1:nsig
            [~,praw] = ttest(testdata(:,it,isig));
            [datobs, datrnd] = lib.fromToolbox.stat.cluster_test_helper(testdata(:,it,isig)', nperm,'t');
            [h, p, clusterinfo] = lib.fromToolbox.stat.cluster_test(datobs, datrnd, 0, 0.05,0.05, 'sum');
            pcorr  = p';
            barcol  = sigcol(isig,:);
            yheight = min(ylims) + range(ylims)*0.05*(isig);
            pheight = min(ylims) + range(ylims)*0.05*(isig+0.5)+0.5;
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
        % clear testdata;
        
        % plot time courses (with shaded error bars)
        cfg = [];
        hplot = lib.fromToolbox.plt.plotpatch(plotdata,plottime,plotcol,cfg);
        legend(hplot,condnames,'fontsize',16), legend boxoff
        xlim(xlims); ylim(ylims);
        
        % set axis ticks, axis labels, etc.
        xstpsize = 1;
        set(gca,'xtick',[0:xstpsize:4]);
        set(gca,'ytick',unique([0 ylims]));
        
        figpath =  'figures/erp'; if ~exist(figpath,'dir'), mkdir(figpath); end
        fname = sprintf('%s/%s%s_%s',figpath,'reward_erp_',roiNames{roi}, ...
            condnames{2-reward});
        print(gcf,'-dpng',fname);
        print(gcf,'-dpdf','-painters',fname);
    end
end

%% Additional functions

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