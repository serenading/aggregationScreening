clear
close all

addpath('auxiliary/')

%% Script reads in featureTable and optionally plots 
% boxplot for all strains (40 and 5 worms) and/or 
% expanded feature values across worm categories and movie phases (40 worm blob feats only).

% author: @serenading. June 2020

%% Set analysis parameters
% set which feature extraction timestamp to use
extractStamp = '20200519_153722';

% which worm density to plot feature for?
wormNum = 40; % 40 or 5

% choose which feature to plot. Cell containing strings that match feature variable names
features2plot = {'foodEdge_notPhase1_fraction'};
% 'cluster_size_20_25min_50th','cluster_size_20_25min_90th','cluster_size_40_45min_50th','cluster_size_40_45min_90th'};
% {'cluster_perdurance_90prc','cluster_perdurance_50prc','cluster_perdurance_10prc'};
%{'blobSpeed_50th','blobSpeed_90th','blobHu0_50th','blobHu1_50th','blobHu2_50th','blobHu3_50th','blobHu4_50th','blobHu5_50th','blobHu6_50th'};%,'blobSpeed_50th','d_blobSpeed_abs_50th','d_blobSpeed_abs_90th'};
% {'blobCompactness_cluster_50th','blobCompactness_cluster_90th','blobSpeed_cluster_50th','blobSpeed_cluster_90th','blobArea_cluster_50th','blobArea_cluster_90th','blobQuirkiness_cluster_50th','blobQuirkiness_cluster_90th','blobSolidity_cluster_50th','blobSolidity_cluster_90th','blobHu0_cluster_50th','blobHu1_cluster_50th','blobHu2_cluster_50th','blobHu3_cluster_50th'};

% what to plot?
plotboxplot = true;
plotPhase4Divergent = false; % only works for 40 worms feature prior to expansion

%% Load the latest featureTable
featureTable = loadLatestFeatureTable(extractStamp,wormNum);

%% Extract value for strain names and the desired features
for featCtr = 1:numel(features2plot)
    featureName = features2plot{featCtr};
    featureVal = featureTable.(featureName);
    
    %% Boxplot
    if plotboxplot
        % extract group stats
        [means,strains,n] = grpstats(featureVal,featureTable.strain_name,{@median,'gname','numel'});
        % concatenate extracted stats
        stats = table(means,strains,n);
        % sort stats by median
        stats_sorted = sortrows(stats);
        
        % go through each strain generate new grouping variables
        strainnames_sorted = NaN(size(featureTable.strain_name));
        colour_sorted = zeros(numel(strainnames_sorted),3); % [0 0 0] gives black default colour
        
        % go through each group
        for grpCtr = 1:numel(stats_sorted.strains)
            strain = stats_sorted.strains(grpCtr);
            strainLogInd = strcmp(featureTable.strain_name,strain);
            % assert(nnz(strainLogInd) == stats_sorted.n(grpCtr)); % assertion will not work if some files contain NaN feature value
            strainnames_sorted(strainLogInd) = grpCtr;
            % generate colour labels of interest
            if strcmp(strain,'N2')
                colour_sorted(grpCtr,:) = [1 0 0]; % red
            elseif strcmp(strain,'CB4856')
                colour_sorted(grpCtr,:) = [0 1 1]; % cyan
            elseif strcmp(strain,'DA609')
                colour_sorted(grpCtr,:) = [0 0 1]; % blue
            elseif strcmp(strain,'NIC256')
                colour_sorted(grpCtr,:) = [1 1 0]; % yellow
            end
        end
        assert(nnz(isnan(strainnames_sorted))==0,'Some strainnames_sorted values are NaN')
        
        % plot
        figure;
        H = boxplot(featureVal,strainnames_sorted,...
            'PlotStyle','compact','BoxStyle','filled',...
            'Labels',stats_sorted.strains,'LabelOrientation','inline','Colors',colour_sorted);
        title(featureName,'Interpreter','none')
        %savefig(['/Users/sding/Desktop/AggregationScreening/fiveWorm/' char(feat2use) '.fig'])
    end
    
    %% Plot expanded feature values for different worm categories across phases
    if plotPhase4Divergent && wormNum==40
        % load divergent strain names
        load('strainsList/divergent.mat','strains');
        % initialise new figure. Each row of the subplot is a worm category, each column is a time phase
        phaseFig = figure;
        % generate colour map for plotting strains
        plotcolors = linspecer(numel(strains));
        % TODO: add something to check for regular expression pattern before splitting string
        featureNameSplit = strsplit(featureName,'_');
        % generate cell with expanded feature names
        expandedFeatureNames = {[featureNameSplit{1},'_phase1_' featureNameSplit{2}],...
            [featureNameSplit{1},'_phase2_' featureNameSplit{2}],...
            [featureNameSplit{1},'_phase3_' featureNameSplit{2}],...
            [featureNameSplit{1},'_notphase1_' featureNameSplit{2}],...
            [featureNameSplit{1},'_' featureNameSplit{2}],...
            [featureNameSplit{1},'_sw_phase1_' featureNameSplit{2}],...
            [featureNameSplit{1},'_sw_phase2_' featureNameSplit{2}],...
            [featureNameSplit{1},'_sw_phase3_' featureNameSplit{2}],...
            [featureNameSplit{1},'_sw_notphase1_' featureNameSplit{2}],...
            [featureNameSplit{1},'_sw_' featureNameSplit{2}],...
            [featureNameSplit{1},'_mw_phase1_' featureNameSplit{2}],...
            [featureNameSplit{1},'_mw_phase2_' featureNameSplit{2}],...
            [featureNameSplit{1},'_mw_phase3_' featureNameSplit{2}],...
            [featureNameSplit{1},'_mw_notphase1_' featureNameSplit{2}],...
            [featureNameSplit{1},'_mw_' featureNameSplit{2}],...
            [featureNameSplit{1},'_cluster_phase1_' featureNameSplit{2}],...
            [featureNameSplit{1},'_cluster_phase2_' featureNameSplit{2}],...
            [featureNameSplit{1},'_cluster_phase3_' featureNameSplit{2}],...
            [featureNameSplit{1},'_cluster_notphase1_' featureNameSplit{2}],...
            [featureNameSplit{1},'_cluster_' featureNameSplit{2}],...
            [featureNameSplit{1},'_pausedmw_phase1_' featureNameSplit{2}],...
            [featureNameSplit{1},'_pausedmw_phase2_' featureNameSplit{2}],...
            [featureNameSplit{1},'_pausedmw_phase3_' featureNameSplit{2}],...
            [featureNameSplit{1},'_pausedmw_notphase1_' featureNameSplit{2}],...
            [featureNameSplit{1},'_pausedmw_' featureNameSplit{2}]};
        % populate cell with expanded feature values
        expandedFeatureVals = NaN(numel(featureVal),numel(expandedFeatureNames));
        for expandedFeatureCtr = 1:numel(expandedFeatureNames)
            expandedFeatureVals(:,expandedFeatureCtr) = featureTable.(expandedFeatureNames{expandedFeatureCtr});
        end
        % subdivide feature by strain
        for strainCtr = 1:numel(strains)
            strain = strains{strainCtr};
            strainLogInd = strcmp(featureTable.strain_name,strain);
            % all worm category
            subplot(5,1,1); hold on
            vals2plot = expandedFeatureVals(strainLogInd,1:5);
            plot([1:5],nanmean(vals2plot,1),'-x','Color',plotcolors(strainCtr,:))
            title('all tracked objects')
            % sw category
            subplot(5,1,2); hold on
            vals2plot = expandedFeatureVals(strainLogInd,6:10);
            plot([1:5],nanmean(vals2plot,1),'-x','Color',plotcolors(strainCtr,:))
            title('single worms')
            % mw category
            subplot(5,1,3); hold on
            vals2plot = expandedFeatureVals(strainLogInd,11:15);
            plot([1:5],nanmean(vals2plot,1),'-x','Color',plotcolors(strainCtr,:))
            title('multi worms')
            % cluster category
            subplot(5,1,4); hold on
            vals2plot = expandedFeatureVals(strainLogInd,16:20);
            plot([1:5],nanmean(vals2plot,1),'-x','Color',plotcolors(strainCtr,:))
            title('clusters')
            % pausedmw category
            subplot(5,1,5); hold on
            vals2plot = expandedFeatureVals(strainLogInd,21:25);
            plot([1:5],nanmean(vals2plot,1),'-x','Color',plotcolors(strainCtr,:))
            title('paused multi worms')
        end
        % format plots
        allaxes = findall(phaseFig, 'type', 'axes');
        linkaxes(allaxes);
        for subplotCtr = 1:numel(allaxes) % findall gets subplot handles in the reverse order
            legend(allaxes(subplotCtr),strains)
            xticks(allaxes(subplotCtr),[1:5])
            xticklabels(allaxes(subplotCtr),{'0-15','15-30','30-45','15-45','0-45'})
            xlabel(allaxes(subplotCtr),{'Recording duration (min)'})
            ylabel(allaxes(subplotCtr),{featureName},'Interpreter','None')
        end
    end
end

%% some 40 worm results:
% ECA252 (high hu2,hu3, low blob speed) , XZ1514 (low hu0,hu1,low quirkiness)