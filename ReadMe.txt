Features summaries files to use:

20200519_153722 feature summaries have 3016 features (without dorsal-ventral features) with short, matlab-friendly names;
20200511_162714 feature summaries have 3016 features (without dorsal-ventral features) with short, matlab-friendly names; for three windows of subdivided 15-minute movies;
20191024_122847 feature summaries have 4548 full features (including dorsal-ventral features). 

Useful scripts:

generateFeatSummary: First script to run. Combines Tierpsy tables with metadata table to generate a combined FeatureTable used for downstream analysis. Subdivides table into forty-worm and five-worm datasets. Appends missing "food_region" features for 20200519_153722 and 20200511_162714 feature summaries.
plotControlPCA: PCA and manova test on three control strains to assess the effects of experimental procedures. Found that strain and time in diapause have an effect; factors unlikely to have an effect include run number, (bleach prep?), camera number. Length of worm as a function of time showed little difference. 
plotAllPCA: PCA on all data using all features, coloured by worm density and strains. 
testForFeatNormality: test features extracted from three control strains at 5 worm density to see how normally distributed they are. Found that most are normal; non-Gaussian features include path curvatures, normalized body lengths, etc.
compareSwFeatDensityEffect: analyses density effects in isolated single (skeletonisable) worms using Tierpsy features from 5 vs. 40 worm experiments. Identifies key features that show density dependency in more than 10% of the strains(, and identifies key strains (none) where more than 10% of the Tierpsy features show density dependency). 
classifyVariable: uses supervised machine learning algorithms to train classifiers for a specified variable based on extracted Tierpsy features. It also has the option to apply sequantial feature selection to identify top features to use for classification.

Useful functions:

preprocessFeatMat: pro-processes features matrix with several steps such as dropping features with zero standard deviation or too many NaN's, imputing nan values to global mean, and z-normalising feature matrix.
filterFeatureTable: filters features table by specified strain and feature requirements before classification tasks. Uses dropFeats and dropStrains.
dropFeats
dropStrains