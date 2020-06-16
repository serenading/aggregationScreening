% Author: @serenading. June 2020

function [phase1LogInd, phase2LogInd, phase3LogInd, notPhase1LogInd] = restrictPhase(trajData,frameRate)

phase1LogInd = trajData.frame_number < frameRate*60*15; % first 15 min
phase2LogInd = (trajData.frame_number < frameRate*60*30) & ~ phase1LogInd; % middle 15 min
phase3LogInd = ~phase1LogInd & ~phase2LogInd; % final 15 min
notPhase1LogInd = ~phase1LogInd; % all but the first 15 min

end