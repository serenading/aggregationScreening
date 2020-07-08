f = '/Users/sding/Desktop/sampleData/1.1_1_da609_Set0_Pos0_Ch5_14012018_092702_featuresN.hdf5';
trajData = h5read(f,'/trajectories_data');
frameBins = 25;
n_frames = max(trajData.frame_number);
n_blob = NaN(1,ceil(n_frames/frameBins));
for frameCtr = 1:frameBins:n_frames
    frames = frameCtr
    frameLogInd = trajData.frame_number == frameCtr;
    n_blob(frameCtr) = numel(unique(trajData.worm_index_joined(frameLogInd)));
end
plot(1:n_frames,n_blob,'-.')
xlabel('frames')
ylabel('n_blob')