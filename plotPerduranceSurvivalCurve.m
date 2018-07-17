

cumSurvivalFig = figure;
legendMatrix=cell(length(strains)*length(wormnums),1);

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        %% load data
        if dataset ==1
            filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_list.txt']);
        elseif dataset ==2
            filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
        end
        numFiles = length(filenames);
        frameDist = [];
        for fileCtr=1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            numCloseNeighbr = h5read(filename,'/num_close_neighbrs');
            neighbrDist = h5read(filename,'/neighbr_distances');
            %% filter data
            % filter green by blob size and intensity
            trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                intensityThresholds(wormnum),maxBlobSize);
            % filter green by small cluster status
            trajData.filtered = trajData.filtered&...
                ((numCloseNeighbr== 2 & neighbrDist(:,3)>=minNeighbrDist)...
                |(numCloseNeighbr== 3 & neighbrDist(:,4)>=(minNeighbrDist))...
                |(numCloseNeighbr== 4 & neighbrDist(:,5)>=(minNeighbrDist)));
            %% generate distribution of small cluster frames
            smallClusterFrames = trajData.frame_number(trajData.filtered)';
            if isempty(smallClusterFrames) == false
                magicInd = diff([0 diff([smallClusterFrames]) 0]==1);
                continuousFrames = find(magicInd == -1) - find(magicInd == 1) + 1; % list lengths of consecutive frames
                singleFrames = length(smallClusterFrames) - sum(continuousFrames(:)); % find number of single frames
                frameDist = [frameDist ones(1,singleFrames) continuousFrames]; % compile distribution of consecutive frames by adding each movie
            end
            % calculate how many clusters disappear due to broken trajectory
                % generate logical index for single frames and last frames
                % of a continuous frame run
            smallClusterFramesLogInd = [diff(smallClusterFrames)~=1, true];
                %list the next frame number after the end of a continuous run of (or a single) frames
            nextFrameList = smallClusterFrames(smallClusterFramesLogInd)+1;
                %list the worm index for those corresponding frames
            smallClusterWorms = trajData.worm_index_joined(trajData.filtered)';
            nextFrameWormList = smallClusterWorms(smallClusterFramesLogInd);
                %check that the worm index still exists in the next frame
                smallClusterContinuesCount=0;
            for nextFrameCtr = 1:length(nextFrameList)
                nextFrame = nextFrameList(nextFrameCtr);
                nextFrameWorm = nextFrameWormList(nextFrameCtr);
                if ismember(nextFrameWorm,trajData.worm_index_joined(find(trajData.frame_number==nextFrame)));
                    smallClusterContinuesCount = smallClusterContinuesCount+1;
                end
            end
            proportionToContinue = smallClusterContinuesCount/nnz(frameDist)*100
            strcat(strain, '\_', wormnum, filename)
        end
        % plot cumulative survival
        [ecdfy,ecdfx] = ecdf(frameDist);
        plot(ecdfx,1-ecdfy) % gives a smoother curve than the survival function
        %ecdf(frameDist,'function','survivor','alpha',0.01,'bounds','on')
        hold on
        legendMatrix{(numCtr-1)*length(strains)+(strainCtr)}= strcat(strain, '\_', wormnum);
    end
end
%% format graphs and export
xlabel('frames elapsed (at 9fps)')
xlim([0 40])
ylabel('remaining proportion')
legend(legendMatrix)
if dataset ==1
    epsFileName = ['figures/smallClusterPerdurance/green1/pdf/smallClusterPerduranceSurvivalPooled.eps'];
    figFileName = ['figures/smallClusterPerdurance/green1/fig/smallClusterPerduranceSurvivalPooled.fig'];
elseif dataset ==2
    epsFileName = ['figures/smallClusterPerdurance/green2/pdf/smallClusterPerduranceSurvivalPooled.eps'];
    figFileName = ['figures/smallClusterPerdurance/green2/fig/smallClusterPerduranceSurvivalPooled.fig'];
end
%savefig(figFileName)
%exportfig(cumSurvivalFig,epsFileName,exportOptions)
%system(['epstopdf ' epsFileName]);
%system(['rm ' epsFileName]);