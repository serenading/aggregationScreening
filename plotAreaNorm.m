f = '/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/fortyWormFeaturesTable_20200519_153722_new_20200620_blobFeats.csv';
featureTable = readtable(f,'Delimiter',',','preserveVariableNames',true);
normMethod = 'probability';

%% calculate normalised areas
blobArea_pausedmw_50th_norm = featureTable.blobArea_pausedmw_50th./featureTable.blobArea_sw_50th;
blobArea_pausedmw_90th_norm = featureTable.blobArea_pausedmw_90th./featureTable.blobArea_sw_50th;
blobArea_cluster_50th_norm = featureTable.blobArea_cluster_50th./featureTable.blobArea_sw_50th;
blobArea_cluster_90th_norm = featureTable.blobArea_cluster_90th./featureTable.blobArea_sw_50th;

% blobArea_pausedmw_foodEdge_50th_norm = featureTable.blobArea_pausedmw_foodEdge_50th./featureTable.blobArea_sw_50th;
% blobArea_pausedmw_foodEdge_90th_norm = featureTable.blobArea_pausedmw_foodEdge_90th./featureTable.blobArea_sw_50th;
% blobArea_cluster_foodEdge_50th_norm = featureTable.blobArea_cluster_foodEdge_50th./featureTable.blobArea_sw_50th;
% blobArea_cluster_foodEdge_90th_norm = featureTable.blobArea_cluster_foodEdge_90th./featureTable.blobArea_sw_50th;

%% plots

% 50th prctile plot
figure; hold on
histogram(blobArea_pausedmw_50th_norm,'BinWidth',0.1,'Normalization',normMethod,'DisplayStyle','stairs')
histogram(blobArea_cluster_50th_norm,'BinWidth',0.1,'Normalization',normMethod,'DisplayStyle','stairs')
legend({'blobArea_pausedmw_50th_norm','blobArea_cluster_50th_norm'},'Interpreter','none')
xlim([0 25]); xticks([0:1:25])
ylabel(normMethod)
% 90th prctile plot
figure; hold on
histogram(blobArea_pausedmw_90th_norm,'BinWidth',0.1,'Normalization',normMethod,'DisplayStyle','stairs')
histogram(blobArea_cluster_90th_norm,'BinWidth',0.1,'Normalization',normMethod,'DisplayStyle','stairs')
legend({'blobArea_pausedmw_90th_norm','blobArea_cluster_90th_norm'},'Interpreter','none')
xlim([0 25]); xticks([0:1:25])
ylabel(normMethod)

% % 50th prctile plot
% figure; hold on
% histogram(blobArea_pausedmw_foodEdge_50th_norm,'BinWidth',0.1,'Normalization','count','DisplayStyle','stairs')
% histogram(blobArea_cluster_foodEdge_50th_norm,'BinWidth',0.1,'Normalization','count','DisplayStyle','stairs')
% legend({'blobArea_pausedmw_foodEdge_50th_norm','blobArea_cluster_foodEdge_50th_norm'},'Interpreter','none')
% xlim([0 25]); xticks([0:1:25])
% ylabel('Probability')
% % 90th prctile plot
% figure; hold on
% histogram(blobArea_pausedmw_foodEdge_90th_norm,'BinWidth',0.1,'Normalization','count','DisplayStyle','stairs')
% histogram(blobArea_cluster_foodEdge_90th_norm,'BinWidth',0.1,'Normalization','count','DisplayStyle','stairs')
% legend({'blobArea_pausedmw_foodEdge_90th_norm','blobArea_cluster_foodEdge_90th_norm'},'Interpreter','none')
% xlim([0 25]); xticks([0:1:25])
% ylabel('Probability')