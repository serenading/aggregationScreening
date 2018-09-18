function trajData_worm_label = setBadFlag(filename)

%% function returns a logical index for bad objects (usually edge of a plate or worm tracks)
%% that are tracked by accident and not possible to manually flag nad for all frames, due to worm index changing too often

% three hand picked worm_index_joined indices for bad object, from different parts of the movie
if contains(filename,'3.2_6_ju1440_a3_Set0_Pos0_Ch3')
    badWormIds = [1,4183,11006];
elseif contains(filename,'9.2_9_qx1792_d7_ps2025_8c_Set0_Pos0_Ch6')
    badWormIds = [2389,4177,7243];
elseif contains(filename,'14.3_7_xz1514_1b_n2_6b_Set0_Pos0_Ch2')
    badWormIds = [20822,39355,61841];
elseif contains(filename,'2.1_5_nic261_f6_Set0_Pos0_Ch2')
    badWormIds = [1558,4346,7380];
elseif contains(filename,'9.2_5_cx11262_fe_Set0_Pos0_Ch4')
    badWormIds = [267,2514,4227];
elseif contains(filename,'12.1_3_my2741_55_Set0_Pos0_Ch2')
    badWormIds = [8,10807,22062];
elseif contains(filename,'13.2_1_nic1107_08_Set0_Pos0_Ch6')
    badWormIds = [678,2773,3827];
end

% find the maximum range of bounding box xy coordinates for the bad entry based on the three worm indices
plateWorms = h5read(filename,'/plate_worms');
for badWormIdCtr = 1:length(badWormIds)
    badWormId = badWormIds(badWormIdCtr);
    badWormIdx = find(plateWorms.worm_index_joined == badWormId);
    badWormIdx = badWormIdx(1);
    bad_xcoords(badWormIdCtr,:) = [plateWorms.bounding_box_xmin(badWormIdx) plateWorms.bounding_box_xmax(badWormIdx)];
    bad_ycoords(badWormIdCtr,:) = [plateWorms.bounding_box_ymin(badWormIdx) plateWorms.bounding_box_ymax(badWormIdx)];
end
bad_x = [min(bad_xcoords(:,1)) max(bad_xcoords(:,2))];
bad_y = [min(bad_ycoords(:,1)) max(bad_ycoords(:,2))];

% generate logical index for bad objects, if their centroid xy coordinates fall within the range of the bounding box
trajData = h5read(filename,'/trajectories_data');
trajData.worm_bad_set = trajData.coord_x >= bad_x(1) & trajData.coord_x <= bad_x(2) & trajData.coord_y >= bad_y(1) & trajData.coord_y <= bad_y(2);
trajData_worm_label = uint8(trajData.worm_bad_set)*3; % make the flag 3 so it matches "trajData.worm_label = 3" for manually set bad flag