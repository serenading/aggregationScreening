function ac = calculateImageAutocorrelation(maskedImageStack2D,frameStd)

%% function calculates pixel image autocorrelation across half of the frames as specified by numFrames
% (so the last frame can extend to maximum lag time).
% Calculation based on formula in Zero-normalized cross-correlation (ZNCC) section of https://en.wikipedia.org/wiki/Cross-correlation
%% INPUT
% maskedImageStack2D: logical matrix, in the dimension of frame number x linearised boolean image pixels for each frame. 
%% OUTPUT
% ac: 1x1 scalar, autocorrelation coefficient.


%% FUNCTION

% get total number of frames contained in 2D image stack
numFrames = size(maskedImageStack2D,2);
% get the image frame total pixel number
pixelNum = size(maskedImageStack2D,1);
% get the total number of starting frames to calculate image autocorrelation for
numStartingFrames = floor(numFrames/2);
% get the number of lag times (tau)
numLags = numStartingFrames;
% generate matrix to hold standard deviation product values between corresponding frames
frameStdMat = frameStd.*frameStd';
% pre-allocate matrix to hold autocorrelation values
ac = NaN(numStartingFrames,numLags);
% go through each starting imaging frame
for frameCtr = 1:numStartingFrames
    % read the starting frame
    thisStartingFrame = maskedImageStack2D(:,frameCtr)'; % row vector 1 by npx
    % calculate image cross covariance
    autocovar = thisStartingFrame * maskedImageStack2D; % longer than necessary row vector. linear algebra magic - Faster than summimg the dot product of untransposed thisStartingFrame
    % divide the desired slice of autocovar by standard deviation to get autocorrelation
    autocorr = autocovar(frameCtr:frameCtr+numLags-1)./frameStdMat(frameCtr,frameCtr:frameCtr+numLags-1);
    % assign value to variable
    ac(frameCtr,:) = autocorr; 
end
% normalise by the number of total pixels from the first frame
ac = ac/pixelNum;