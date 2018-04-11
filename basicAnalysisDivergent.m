clear
close all

%% script extracts perimeter and area multi-worm blob features for each of the strains of interest, normalise them against corresponding single worm blob feature within each replicate, and plots a histogram for each strain

%% set parameters
% load('strainsList/controls.mat')
load('strainsList/divergent.mat')
saveResults = true;

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',30,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',25,...
    'LineWidth',3);

% get list of file names for each strain
[strainFileList] = getFileList(strains);
% create empty figures
perimeterFig = figure; hold on
areaFig = figure; hold on
% generate colormap for plotting each strain
colorMap = distinguishable_colors(length(strains));

%% go through each strain
for strainCtr = 1:length(strains)
    filenames = strainFileList.([strains{strainCtr} 'List_40']);
    %% go through each recording
    for fileCtr = 1:length(filenames)
        
        %% load data
        filename = filenames{fileCtr}
        % ignore the two problematic files for the divergent set
        if ~strcmp(filename,'/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/agg_2.2_180118/2.2_5_cb4856_oo_Set0_Pos0_Ch5_18012018_133212_skeletons.hdf5')...
                & ~strcmp(filename,'/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/agg_2.2_180118/2.2_3_jt11398_77_Set0_Pos0_Ch5_18012018_112925_skeletons.hdf5')
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton'); % in pixels
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            % features = h5read(strrep(filename,'skeletons','featuresN'),'/features_timeseries/');
            
            %% look at proportion of worms with good skeleton
            % find worm clusters vs. single worms based on skeletonisation
            multiWormLogInd = logical(~trajData.is_good_skel);
            singleWormLogInd = logical(trajData.is_good_skel);
            % read features
            perimeter.(strains{strainCtr}){fileCtr} = blobFeats.perimeter(multiWormLogInd);
            area.(strains{strainCtr}){fileCtr} = blobFeats.area(multiWormLogInd);
            swPerimeter.(strains{strainCtr}){fileCtr} = blobFeats.perimeter(singleWormLogInd);
            swArea.(strains{strainCtr}){fileCtr} = blobFeats.area(singleWormLogInd);
            % normalise area and perimeter from this movie with sw features from this movie; store value for threshold box plot later
            perimeterNorm.(strains{strainCtr}){fileCtr} = perimeter.(strains{strainCtr}){fileCtr}/median(swPerimeter.(strains{strainCtr}){fileCtr});
            areaNorm.(strains{strainCtr}){fileCtr} = area.(strains{strainCtr}){fileCtr}/median(swArea.(strains{strainCtr}){fileCtr});
            % remove low values below 1 as clusters should by definition be larger than single worms thus min normalised value should be 1
            perimeterNorm.(strains{strainCtr}){fileCtr} =  perimeterNorm.(strains{strainCtr}){fileCtr}( perimeterNorm.(strains{strainCtr}){fileCtr}>1);
            areaNorm.(strains{strainCtr}){fileCtr} = areaNorm.(strains{strainCtr}){fileCtr}(areaNorm.(strains{strainCtr}){fileCtr}>1);
            
%             %% plot
%             set(0,'CurrentFigure',perimeterFig)
%             histogram(perimeterNorm.(strains{strainCtr}){fileCtr},'Normalization','pdf','DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
%             set(0,'CurrentFigure',areaFig)
%             histogram(areaNorm.(strains{strainCtr}){fileCtr},'Normalization','pdf','DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
        end
    end
    % combine the replicates
    perimeterNormPooled.(strains{strainCtr}) = perimeterNorm.(strains{strainCtr}){1}; % initialise new pooling variable with the first replicate
    areaNormPooled.(strains{strainCtr}) = areaNorm.(strains{strainCtr}){1};
    for fileCtr = 2:length(filenames)
        perimeterNormPooled.(strains{strainCtr}) = vertcat(perimeterNormPooled.(strains{strainCtr}),perimeterNorm.(strains{strainCtr}){fileCtr});
        areaNormPooled.(strains{strainCtr}) = vertcat(areaNormPooled.(strains{strainCtr}),areaNorm.(strains{strainCtr}){fileCtr});
    end
    set(0,'CurrentFigure',perimeterFig)
    histogram(perimeterNormPooled.(strains{strainCtr}),'Normalization','pdf','DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
    set(0,'CurrentFigure',areaFig)
    histogram(areaNormPooled.(strains{strainCtr}),'Normalization','pdf','DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
    
end

% format and export histogram
set(0,'CurrentFigure',perimeterFig)
set(gca, 'YScale', 'log')
legend(strains,'Location','eastoutside')
title('cluster perimeter')
xlabel('relative perimeter')
ylabel('probability')
figurename = 'figures/perimeter_log';
if saveResults
    exportfig(perimeterFig,[figurename '.eps'],exportOptions)
    %eps2pdf([figurename '.eps'])
    %system(['epstopdf ' figurename '.eps']);
    %system(['rm ' figurename '.eps']);
end

set(0,'CurrentFigure',areaFig)
set(gca, 'YScale', 'log')
legend(strains,'Location','eastoutside')
title('cluster area')
xlabel('relative area')
ylabel('probability')
figurename = 'figures/area_log';
if saveResults
    exportfig(areaFig,[figurename '.eps'],exportOptions)
    %eps2pdf([figurename '.eps'])
    %system(['epstopdf ' figurename '.eps']);
    %system(['rm ' figurename '.eps']);
end

% save data
if saveResults
    save('results/areaPerimeter.mat','perimeterNorm','areaNorm')
    savefig(perimeterFig,'figures/perimeter.fig')
    savefig(areaFig,'figures/area.fig')
end