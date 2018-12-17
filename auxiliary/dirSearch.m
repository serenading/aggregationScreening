function [fileList wormNames] = dirSearch(directory, fileSuffix)

% DIRSEARCH  Search the specified directory return all the files in its
% subdirectories that contain files ending in fileSuffix. Also extract worm
% names from the directory structure. NB - worm name extraction will of
% course only work for the specific structure we use:
% directory/mutant name/allele/on_food/.../fileName
% 
% Uses rdir function by Gus Brown downloaded from MATLAB Central
%
% Input
%   directory  - The directory to search for files
%   fileSuffix - The file ending unique to the set of files to be extracted
%
% Output
%   fileList   - The full paths to all files in directory ending in
%                fileSuffix
%   wormNames  - The worm names corresponding to each file found of the
%                form mutant-name_allele_on_food (or off_food)
% 
% 
% Copyright Medical Research Council 2013
% André Brown, andre.brown@csc.mrc.ac.uk, aexbrown@gmail.com
% 
% 
% The MIT License
% 
% Copyright (c)  Medical Research Council 2013
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.



% get all the projectedAmpsNoNaN files in the directory
% slashes must be treated differently on mac and pc so check if you're on a
% pc.
if ispc
    timeSeriesNames = rdir([directory '**\*' fileSuffix]);
else
    timeSeriesNames = rdir([directory '**/*' fileSuffix]);
end

% initialise file list and worm names
fileList = cell(length(timeSeriesNames), 1);
wormNames = cell(length(timeSeriesNames), 1);

for i = 1:length(timeSeriesNames)
    % get the name of the current time series
    timeSeriesName = timeSeriesNames(i).name;
    fileList{i} = timeSeriesName;
    
    % extract the mutant, allele, and whether it's on or off food (i.e.
    % everything in the name from the root directory up to the
    % specification of the worm side)
    if ispc
        wormNameDirectory = regexpi(timeSeriesName, ...
            ['(?<=' directory ').+(?=\\on\_food\\)'], 'match');
        wormNames{i} = char(strrep(wormNameDirectory, '\', '_'));
    else
        wormNameDirectory = regexpi(timeSeriesName, ...
            ['(?<=' directory ').+(?=/on\_food/)'], 'match');
        wormNames{i} = char(strrep(wormNameDirectory, '/', '_'));
    end
end