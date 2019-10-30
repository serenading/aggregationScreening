%% generate tif images for hand labeling
% open file
fid = fopen('strainsList/noFoodContourFiles_skel_40.txt');
filename = 'initialise';
% read text file one line at a time
while ischar(filename)
    % move onto the next line
    filename = fgetl(fid);
    % get path to the MaskedVideo file
    maskedVideoFileName = strrep(filename,'Results','MaskedVideos');
    maskedVideoFileName = strrep(maskedVideoFileName,'_featuresN.hdf5','.hdf5');
    % read full image from the MaskedVideo
    fullData = h5read(maskedVideoFileName,'/full_data');
    firstFullImage = fullData(:,:,1);
    % save first image
    splitMaskedVideoFileName = strsplit(maskedVideoFileName,'/');
    imageFileName1 = splitMaskedVideoFileName{end-1};
    imageFileName2 = splitMaskedVideoFileName{end};
    imageFileName2 = strrep(imageFileName2,'.hdf5','.jpg');
    imageFileName = ['/Volumes/behavgenom_archive$/Serena/AggregationScreening/Auxiliary/manualFoodContourImages/aggScreening' imageFileName1 '__' imageFileName2];
    imwrite(firstFullImage,imageFileName);
end
fclose(fid);

%% hand label food contour using VGG annotator (http://www.robots.ox.ac.uk/~vgg/software/via/via.html) and save annotations

%% interpolate xy coordinates to generate the correct coordinate format
addpath('auxiliary/')
% read xy coordinates
[~,~,annotations] = xlsread('/Volumes/behavgenom_archive$/Serena/AggregationScreening/Auxiliary/food_contour.xlsx');
for imageCtr = 1:length(annotations)
    imageCtr
    annotation = annotations{imageCtr};
    annotationSplit = strsplit(annotation,'.jpg');
    coordsRaw = strsplit(annotationSplit{2},'[');
    xcoordsRaw = strsplit(coordsRaw{3},']');
    ycoordsRaw = strsplit(coordsRaw{2},']');
    xcoords = xcoordsRaw{1};
    ycoords = ycoordsRaw{1};
    xcoords = strsplit(xcoords,',');
    ycoords = strsplit(ycoords,',');
    xcoords = cellfun(@str2double,xcoords)';
    ycoords = cellfun(@str2double,ycoords)';
    % interpolate xy coordinates using Euclidian distances
    format short g
    foodCntCoords = interparc(0:1/834:1,xcoords,ycoords,'linear')';
    
    %% write food contour coordinates to skeletons.hdf5
    filename = annotationSplit{1};
    filename = ['/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/' strrep([filename,'_skeletons.hdf5'],'__','/')];
    try
        foodContourCoords = h5read(filename,'/food_cnt_coord');
    catch ME
        if strcmp(ME.identifier,'MATLAB:imagesci:h5read:libraryError')
            h5create(filename,'/food_cnt_coord',size(foodCntCoords),'Datatype','double');
            h5write(filename,'/food_cnt_coord',double(foodCntCoords));
        end
    end
    
    %% delete bad feature files so they can be re-calculated
    filename = strrep(filename,'_skeletons','_featuresN');
    delete(filename)
end