function [] = ModifyDataFilesForLaneDiscrepantDetectors(filename,detectorsToAddZeroLanes,howManyZero,detectorsToRemoveFirstLane,detectorsToAddAvgLanes,howMany)

load(filename)
save([filename '_backupB4Modifications'])

k = 1;
for ID = detectorsToRemoveFirstLane
    
    ind = find(pems.vds == ID);
    pems.data(ind).flw = pems.data(ind).flw(:,2:end);
    pems.data(ind).occ = pems.data(ind).occ(:,2:end);
    pems.data(ind).spd = pems.data(ind).spd(:,2:end);
    k = k + 1;
    
end

k = 1;
for ID = detectorsToAddAvgLanes
    
    ind = find(pems.vds == ID);
    avgFlw = mean(pems.data(ind).flw,2);
    avgOcc = mean(pems.data(ind).occ,2);
    avgSpd = mean(pems.data(ind).spd,2);
    pems.data(ind).flw = [pems.data(ind).flw repmat(avgFlw,1,howMany(k))];
    pems.data(ind).occ = [pems.data(ind).occ repmat(avgOcc,1,howMany(k))];
    pems.data(ind).spd = [pems.data(ind).spd repmat(avgSpd,1,howMany(k))];
    k = k + 1;
    
end

k = 1;
for ID = detectorsToAddZeroLanes
    
    ind = find(pems.vds == ID);
    pems.data(ind).flw = [pems.data(ind).flw zeros(size(pems.data(ind).flw,1),howManyZero(k))];
    pems.data(ind).occ = [pems.data(ind).occ zeros(size(pems.data(ind).occ,1),howManyZero(k))];
    pems.data(ind).spd = [pems.data(ind).spd zeros(size(pems.data(ind).spd,1),howManyZero(k))];
    k = k + 1;
    
end

save(filename)