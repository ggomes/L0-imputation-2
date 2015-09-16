function CellData = InsertBogusMainlineData(CellData,BogusData)

for i = 1:length(CellData)
    
   % find the first link in cell in the BogusData structure
   ind = CellData{i}.LinkIDs(1) == [BogusData.link_id];
   % insert values
   if ~isempty(find(ind,1))
       CellData{i}.MLspeed.Data = BogusData.speeds(ind,:)';
       CellData{i}.MLflow.Data = BogusData.inflows(ind,:)';
       CellData{i}.MLdensity.Data = BogusData.densities(ind,:)';
   end
    
end

% Increase capacity of last cell
CellData{end}.MLFD.q_max = 10000;
CellData{end}.MLFD.vf = 70;
CellData{end}.MLFD.w = 10;
CellData{end}.MLFD.rj = 10000/70+10000/10;
CellData{end}.MLFD.rho_crit = 10000/70;

