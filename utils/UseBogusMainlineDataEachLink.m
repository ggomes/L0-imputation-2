function CellData = UseBogusMainlineDataEachLink(CellData,BogusData,LinkStructure)
% Breaks down the celldata structure back to Link denomination and uses the
% artificial data as measurements for each link. Cells with major onramps
% have a different fundamental diagram assignment rule which ideally should
% be re-evaluated for this function. This is not done inside this function
% yet, so future users beware!! Ideally, fundamental diagrams should be 
% re-calibrated using the new artificial data

% Also, this function is written to work with the specific data file 
% structure. Another case not handled by this function is the HOV.

OutCellData = cell(1,length(BogusData.link_id));

newCellIndex = 1;

for i = 1:length(CellData)-1
    
   % Create a new cell for each link in a cell
   for j = 1:length(CellData{i}.LinkIDs)
       
       % first replicate all fields to the new cell
       OutCellData{newCellIndex} = CellData{i};
       
       % then adjust the fields that are supposed to change
       OutCellData{newCellIndex}.Links = CellData{i}.Links(j);
       OutCellData{newCellIndex}.LinkIDs = CellData{i}.LinkIDs(j);
       OutCellData{newCellIndex}.ORlinkIDs = CellData{i}.ORlinkIDs(j);
       OutCellData{newCellIndex}.FRlinkIDs = CellData{i}.FRlinkIDs(j);
       OutCellData{newCellIndex}.MLLanes = CellData{i}.MLLanes(j);
       OutCellData{newCellIndex}.HVLanes = CellData{i}.HVLanes(j);
       OutCellData{newCellIndex}.linklength = CellData{i}.linklength(j);
       OutCellData{newCellIndex}.orflwnonimputed = CellData{i}.orflwnonimputed(j);
       OutCellData{newCellIndex}.imputeor = CellData{i}.imputeor(j);
       OutCellData{newCellIndex}.frflwnonimputed = CellData{i}.frflwnonimputed(j);
       OutCellData{newCellIndex}.imputefr = CellData{i}.imputefr(j);
       OutCellData{newCellIndex}.orperlink = CellData{i}.orperlink(j);
       OutCellData{newCellIndex}.frperlink = CellData{i}.frperlink(j);
       
       % Insert Bogus Data
       ind = OutCellData{newCellIndex}.LinkIDs(1) == [BogusData.link_id];
       OutCellData{newCellIndex}.MLspeed.Data = BogusData.speeds(ind,:)';
       OutCellData{newCellIndex}.MLflow.Data = BogusData.inflows(ind,:)';
       OutCellData{newCellIndex}.MLdensity.Data = BogusData.densities(ind,:)';
       
       newCellIndex = newCellIndex + 1;

   end
    
end

CellData = OutCellData;

% Increase capacity of last cell
CellData{end}.MLFD.q_max = 10000;
CellData{end}.MLFD.vf = 70;
CellData{end}.MLFD.w = 10;
CellData{end}.MLFD.rj = 10000/70+10000/10;
CellData{end}.MLFD.rho_crit = 10000/70;

