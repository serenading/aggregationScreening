function [sw_perdurance,mw_perdurance,cluster_perdurance] = calculatePerdurance(sw_feature,mw_feature,cluster_feature,strain,fileCtr)

%% function calculates cluster perdurance

swUniqueInd = unique(sw_feature.(strain){fileCtr});
sw_perdurance = NaN(length(swUniqueInd),1);
for swCtr = 1:length(swUniqueInd)
    swInd = swUniqueInd(swCtr);
    sw_perdurance(swCtr) = numel(find(sw_feature.(strain){fileCtr} == swInd));
end

mwUniqueInd = unique(mw_feature.(strain){fileCtr});
mw_perdurance = NaN(length(mwUniqueInd),1);
for mwCtr = 1:length(mwUniqueInd)
    mwInd = mwUniqueInd(mwCtr);
    mw_perdurance(mwCtr) = numel(find(mw_feature.(strain){fileCtr} == mwInd));
end

clusterUniqueInd = unique(cluster_feature.(strain){fileCtr});
cluster_perdurance = NaN(length(clusterUniqueInd),1);
for clusterCtr = 1:length(clusterUniqueInd)
    clusterInd = clusterUniqueInd(clusterCtr);
    cluster_perdurance(clusterCtr) = numel(find(cluster_feature.(strain){fileCtr} == clusterInd));
end