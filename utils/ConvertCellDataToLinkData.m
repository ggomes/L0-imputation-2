function [LinkData] = ConvertCellDataToLinkData(SplitCellData,AfterSplits,LinkStructure)

% Get data from the AfterSplits structure
names = fieldnames(AfterSplits);
for i = 1:length(names)
    eval([names{i} '= AfterSplits.' names{i} ';']);
end
tratio = round(DemandT/SimT);

% Find Total Number of Links in the CellData
totalNumOfLinks = 0;
for i = 1:length(SplitCellData)-1
    totalNumOfLinks = totalNumOfLinks + length(SplitCellData{i}.Links);
end

% Required Parts in the Output
IDs = [];
Lengths = [];
NumOfLanes = [];
Speeds = cell(totalNumOfLinks,1);
Densities = cell(totalNumOfLinks,1);
Flows = cell(totalNumOfLinks,1);
ORIDs = cell(totalNumOfLinks,1);
FRIDs = cell(totalNumOfLinks,1);
ORDemands = cell(totalNumOfLinks,1);
FRFlows = cell(totalNumOfLinks,1);
Capacities = [];
FreeFlowSpeeds = [];
CongestionSpeeds = [];
JamDensities = [];

linkCounter = 1;
for i = 1:length(SplitCellData)-1
    
   linkCounterUp = linkCounter;
   for j = 1:length(SplitCellData{i}.LinkIDs)
       IDs = [IDs; SplitCellData{i}.LinkIDs(j)];
       Lengths = [Lengths; SplitCellData{i}.linklength(j)];
       NumOfLanes = [NumOfLanes; SplitCellData{i}.MLLanes(j)];
       Speeds{linkCounter} = mean(reshape(Velocity(:,i),tratio,288));
       Densities{linkCounter} = mean(reshape(Nh(:,i),tratio,288))./LinkLength(i);
       Flows{linkCounter} = mean(reshape(FlowCompare(:,i),tratio,288))/SimT;
       Capacities(linkCounter) = Qmax(i)/SimT;
       FreeFlowSpeeds(linkCounter) = vf(i)/SimT*sum(SplitCellData{i}.linklength);
       CongestionSpeeds(linkCounter) = w(i)/SimT*sum(SplitCellData{i}.linklength);
       JamDensities(linkCounter) = rhojam(i)/sum(SplitCellData{i}.linklength);
       linkCounter = linkCounter + 1;
   end
   
   if ~isempty(onrampIDs{i}) && i~=1
       ORIDs{linkCounterUp} = onrampIDs{i};
       ORDemands{linkCounterUp} = mean(reshape(ORINP(:,i),tratio,288))/SimT;
   end
   
   if ~isempty(offrampIDs{i})
       FRIDs{linkCounter-1} = offrampIDs{i};
       FRFlows{linkCounter-1} = mean(reshape(FrFlow(:,i),tratio,288))/SimT;
   end
    
end

%% Last Cell
linkCounterUp = linkCounter;
for j = 1:length(SplitCellData{end}.LinkIDs)
    IDs = [IDs; SplitCellData{end}.LinkIDs(j)];
    Lengths = [Lengths; SplitCellData{end}.linklength(j)];
    NumOfLanes = [NumOfLanes; SplitCellData{end}.MLLanes(j)];
    Speeds{linkCounter} = mean(reshape(Velocity(:,end),tratio,288));
    Densities{linkCounter} = mean(reshape(Nh(:,end),tratio,288))./LinkLength(end);
    Flows{linkCounter} = mean(reshape(FlowCompare(:,end),tratio,288))/SimT;
    Capacities(linkCounter) = Qmax(end)/SimT;
    FreeFlowSpeeds(linkCounter) = vf(end)/SimT*sum(SplitCellData{end}.linklength);
    CongestionSpeeds(linkCounter) = w(end)/SimT*sum(SplitCellData{end}.linklength);
    JamDensities(linkCounter) = rhojam(end)/sum(SplitCellData{end}.linklength);
    linkCounter = linkCounter + 1;
end

if ~isempty(onrampIDs{end})
    ORIDs{linkCounterUp} = onrampIDs{end};
    ORDemands{linkCounterUp} = mean(reshape(ORINP(:,end),tratio,288))/SimT;
end

if ~isempty(offrampIDs{end})
    FRIDs{linkCounterUp} = offrampIDs{end};
    FRFlows{linkCounterUp} = mean(reshape(FrFlow(:,end),tratio,288))/SimT;
end

%% Zero'th cell
% find omitted links
for i = 1:length(LinkStructure)
    if LinkStructure(i).linkid == CellData{1}.LinkIDs(1)
        break;
    end
end
numberOfOmittedLinks = i-1;
% add zero slots to beginning of arrays
IDs = [zeros(numberOfOmittedLinks,1); IDs];
Lengths = [zeros(numberOfOmittedLinks,1); Lengths];
NumOfLanes = [zeros(numberOfOmittedLinks,1); NumOfLanes];
Speeds = [mat2cell(zeros(numberOfOmittedLinks,1),ones(numberOfOmittedLinks,1),1); Speeds];
Densities = [mat2cell(zeros(numberOfOmittedLinks,1),ones(numberOfOmittedLinks,1),1); Densities];
Flows = [mat2cell(zeros(numberOfOmittedLinks,1),ones(numberOfOmittedLinks,1),1); Flows];
ORIDs = [mat2cell(nan(numberOfOmittedLinks,1),ones(numberOfOmittedLinks,1),1); ORIDs];
FRIDs = [mat2cell(nan(numberOfOmittedLinks,1),ones(numberOfOmittedLinks,1),1); FRIDs];
ORDemands = [mat2cell(nan(numberOfOmittedLinks,1),ones(numberOfOmittedLinks,1),1); ORDemands];
FRFlows = [mat2cell(nan(numberOfOmittedLinks,1),ones(numberOfOmittedLinks,1),1); FRFlows];
Capacities = [zeros(1,numberOfOmittedLinks) Capacities];
FreeFlowSpeeds = [zeros(1,numberOfOmittedLinks) FreeFlowSpeeds];
CongestionSpeeds = [zeros(1,numberOfOmittedLinks) CongestionSpeeds];
JamDensities = [zeros(1,numberOfOmittedLinks) JamDensities];

for j = 1:numberOfOmittedLinks
    IDs(j) = LinkStructure(j).linkid;
    Lengths(j) = LinkStructure(j).Length;
    NumOfLanes(j) = LinkStructure(j).MLLanes;
    Speeds{j} = mean(reshape(Velocity(:,1),tratio,288));
    Densities{j} = mean(reshape(Nh(:,1),tratio,288))./Lengths(j);
    Flows{j} = mean(reshape(FlowCompare(:,1),tratio,288))/SimT;
    Capacities(j) = Qmax(1)/SimT;
    FreeFlowSpeeds(j) = vf(1)/SimT*sum(SplitCellData{1}.linklength);
    CongestionSpeeds(j) = w(1)/SimT*sum(SplitCellData{1}.linklength);
    JamDensities(j) = rhojam(1)/sum(SplitCellData{1}.linklength);
end

% if ~isempty(onrampIDs{end})
%     ORIDs{linkCounterUp} = onrampIDs{end};
%     ORDemands{linkCounterUp} = mean(reshape(ORINP(:,end-1),tratio,288))/SimT;
% end
% 
% if ~isempty(offrampIDs{end})
%     FRIDs{linkCounterUp} = offrampIDs{end};
%     FRFlows{linkCounterUp} = mean(reshape(FrFlow(:,end),tratio,288))/SimT;
% end

%% Data Structure Assembly
LinkData.LinkIDs = IDs;
LinkData.LinkLengths = Lengths;
LinkData.NumOfLanes = NumOfLanes;
LinkData.Speeds = Speeds;
LinkData.Densities = Densities;
LinkData.Flows = Flows;
LinkData.OnrampIDs = ORIDs;
LinkData.OnrampDemands = ORDemands;
LinkData.OfframpIDs = FRIDs;
LinkData.OfframpFlows = FRFlows;
LinkData.Capacities = Capacities;
LinkData.FreeFlowSpeeds = FreeFlowSpeeds;
LinkData.CongestionSpeeds = CongestionSpeeds;
LinkData.JamDensities = JamDensities;