Data location:

Data are copied over two separate hard drives. All 5 worm dataset is on HD1. 40 worm data for the divergent strains are on HD1, and on HD2 for the rest of the strains.

Features summaries files to use:

20200630_171156 augmented feature summaries have xxxx features (without dorsal-ventral features) with short, matlab-friendly names. Augmented with 5 fold at 0.8 trajectory ratio, up to 30 min. Windows: window0: 0-15 min window1: 15-30 min, window2: 30-45 min, window3: 15-45 min, window 4: 0-45 min.
20200519_153722 feature summaries have 3016 features (without dorsal-ventral features) with short, matlab-friendly names.
20200511_162714 feature summaries have 3016 features (without dorsal-ventral features) with short, matlab-friendly names; for three windows of subdivided 15-minute movie. Windows: window0: 0-15 min, window1: 15-30 min, window2: 30-45 min).
20191024_122847 feature summaries have 4548 full features (including dorsal-ventral features). 


Useful scripts:

generateFeatSummary: First script to run. Combines Tierpsy tables with metadata table to generate a combined FeatureTable used for downstream analysis. Subdivides table into forty-worm and five-worm datasets. Appends missing "food_region" features for 20200519_153722 and 20200511_162714 feature summaries.
plotControlPCA: PCA and manova test on three control strains to assess the effects of experimental procedures. Found that strain and time in diapause have an effect; factors unlikely to have an effect include run number, (bleach prep?), camera number. Length of worm as a function of time showed little difference. 
plotAllPCA: PCA on all data using all features, coloured by worm density and strains. 
testForFeatNormality: test features extracted from three control strains at 5 worm density to see how normally distributed they are. Found that most are normal; non-Gaussian features include path curvatures, normalized body lengths, etc.
compareSwFeatDensityEffect: analyses density effects in isolated single (skeletonisable) worms using Tierpsy features from 5 vs. 40 worm experiments. Identifies key features that show density dependency in more than 10% of the strains(, and identifies key strains (none) where more than 10% of the Tierpsy features show density dependency). 
classifyVariable: uses supervised machine learning algorithms to train classifiers for a specified variable based on extracted Tierpsy features. It also has the option to apply sequantial feature selection to identify top features to use for classification.
compareGenoPhenoDM: compares distance matrices  based on genotype (SNPs) and phenotype (Tierpsy features) across strains. Useful option: makeMappingFile: outputs a .tsv file in the correct format for cegwas2-nf mapping

Feature calculation scripts and functions:

calculateBlobSpeed: Function calculates blobSpeed (smoothed over 1 second unless otherwise specified) and d_blobSpeed (over dT window of 1/3 second unless otherwise specified, calculated from smoothed speeds).
calculateClusterPerduranceAndSizeDist: Script extracts cluster perdurance and cluster size distribution statistics and optionally generates diagnostic plots for the divergent panel.
calculateBlobFeatures: Script extracts Tierpsy blob features from 40 worm tracking data to generate a feature table for joining onto the master feature table.
calculateDivStat: Script extracts div statistics from tracking data to generate a feature table for joining onto the master feature table. Script is based on calculateBlobFeatures.m.

expandBlobFeature: Function expands base blob feature based on stats, worm type, and movie phase.
expandBlobFeature2: Function expands base blob feature based on stats, worm type, and food region.
expandDivFeature: Function expands base div feature based on stats, worm type, food region, and movie phase.

Feature plotting scripts:
plotFeatsFromTable: reads in featureTable and optionally plots boxplots for all strains (40 and 5 worms) and/or expanded feature values across worm categories and movie phases (40 worm only).
plotDivStats: reads in featureTable and plots stacked bar graphs for blob category and food region statistics to show comprising fractions of tracking indices (40 worm only).


Useful functions:

loadLatestFeatureTable: loads the latest featureTable, because these get updated all the time as new features are added.
extractHDLocation: finds which 40 worm files are on which HD and appends the information to featureTable under "onHD1" heading. Because the full dataset is spread over two HD's. 
appendFeatsToFeatureTable: appends new features from newFeatureTable to the existing featureTable.
preprocessFeatMat: pro-processes features matrix with several steps such as dropping features with zero standard deviation or too many NaN's, imputing nan values to global mean, and z-normalising feature matrix.
filterFeatureTable: filters features table by specified strain and feature requirements before classification tasks. Uses dropFeats and dropStrains.
dropFeats
dropStrains