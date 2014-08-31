clear 
close all

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src','def'))

load(fullfile(root,'analysis','fwy'))
day = 735633;
goodthresh = 0.7;
cell_array = build_cell_array(fwy,day,goodthresh,'use_first_segment');
load(fullfile(root,'analysis','AfterFilledCellData_test1_03-Feb-2014'),'imputation_celldata')
time = 0:300:86100;
% figure
% plot(cell_array.cells(1).segments(1).ml_link.sensor.get_5min_data(day,time,'flw'),'k','Linewidth',3)
% hold on
% plot(imputation_celldata{1}.MLflow.Data,'r--','LineWidth',2)

load(fullfile(root,'analysis','imputation_input_test1_03-Feb-2014'))
get_time = 3600*( Time(1):SimT:(Time(2)-SimT) );

%% Flow % [720x19] 1.2588 0.9783 1.2917 1.0691 1.3340
X = Flow;
x = nan(720,19);
for i=1:length(cell_array.cells)
    sensor = cell_array.cells(i).segments(1).ml_link.sensor;
    x(:,i) = sensor.get_5min_data_time_type(day,get_time,'flw')*SimT;
end
figure
plot(get_time,X(:,1),'k')
hold on 
plot(get_time,x(:,1),'r')
clear X x

%% OK Density % [720x19] 6.1420   13.8526    5.6753    4.7285   10.5255 
X = Density;
x = nan(720,19);
for i=1:length(cell_array.cells)
    cell = cell_array.cells(i);
    sensor = cell_array.cells(i).segments(1).ml_link.sensor;
    LinkLen = cell.get_total_length_miles;
    x(:,i) = sensor.get_5min_data_time_type(day,get_time,'dty')*LinkLen;
end
figure
plot(get_time,X(:,1),'k')
hold on 
plot(get_time,x(:,1),'r')
clear X x

%% OK Speed % [720x19] 60.1183   59.1879   72.5572   63.3952   63.2871
X = Speed;
x = nan(720,19);
for i=1:length(cell_array.cells)
    sensor = cell_array.cells(i).segments(1).ml_link.sensor;
    x(:,i) = sensor.get_5min_data_time_type(day,get_time,'spd');
end
figure
plot(get_time,X(:,1),'k')
hold on 
plot(get_time,x(:,1),'r')
clear X x

% OrPresent % [1x19] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0
X = OrPresent;
x = nan(1,19);
for i=1:length(cell_array.cells)
    x(i) = cell_array.cells(i).has_or;
end
if(any(X~=x))
    disp('failed: OrPresent')
end

% FrPresent % [1x19] 0 1 0 1 1 0 1 1 1 1 1 0 1 1 0 1 1 1 0
X = FrPresent;
x = nan(1,19);
for i=1:length(cell_array.cells)
    x(i) = cell_array.cells(i).has_fr;
end
if(any(X~=x))
    disp('failed: FPresent')
end

%% ImputeOR % [1x19] 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 1 1 0 0
%% ImputeFR % [1x19] 0 0 0 0 0 0 1 1 1 1 1 0 1 1 0 1 0 1 0

%% OK CellLengths % [1x19] 0.4056 1.1590 0.4425 0.3898 0.6946
X = CellLengths;
x = zeros(1,19);
for i=1:length(cell_array.cells)
    cell = cell_array.cells(i);
    x(i) = cell.get_total_length_miles;
end
if(max(abs(X-x))>0.015)
    disp('failure: CellLengths')
end

%% Demand_Giv % [720x19] 0 0.0367 0.1026 0.0543 0.0169
 % BETA_Giv % [720x19] 0 0.0874  0 0.0804 0.0216
% [Demand,BETA] = deriveDemandsAndSplitsFromRampFlows(Impute,Density,OrFlow,FrFlow,STime,Qmax,vf,w,rhojam);

%% ORBoundsLower % [720x18]  0 0 0 0 0
 % ORBoundsUpper % [720x18] 0.0367 0.1026 0.0543 0.0169 0.0304
 % FRBoundsLower % [720x18] 0 0 0 0 0
 % FRBoundsUpper % [720x18]  0 0.0977  0 0.0943 0.0323
 % DemandLower % [720x19]  0 0 0 0 0
 % DemandUpper % [720x19]  0 0.0367 0.1026 0.0543 0.0169
 % BETALower % [720x19] 0 0 0 0 0
 % BETAUpper % [720x19] 0 0.0874  0 0.0804 0.0216
% [ORBoundsLower,ORBoundsUpper,FRBoundsLower,FRBoundsUpper,DemandLower,DemandUpper,BETALower,BETAUpper] = ProcessBounds(Bounds,CellData,numtime,numcell,Impute,SimT,STime,OrFlow,FrFlow,Qmax,w,rhojam,vf,Density,ptr);

% OrFlow_Giv % [720x19] 0.0190 0.0367 0.1026 0.0543 0.0169
X = OrFlow_Giv;
x = nan(720,19);
for i=1:length(cell_array.cells)
    cell = cell_array.cells(i);
    xx = zeros(720,1);
    for j=1:length(cell.segments)
        or_link = cell.segments(j).or_links;
        if(~isempty(or_link) & ~isempty(or_link.sensor))
            xx = xx + or_link.sensor.get_5min_data_time_type(day,get_time,'flw')*SimT;
        end
    end
    x(:,i) = xx;
end
figure
plot(get_time,X(:,1),'k')
hold on 
plot(get_time,x(:,1),'r')
clear X x


%% FrFlow_Giv % [720x19]  0 0.0977  0 0.0943 0.0323

%% Capacities % [1x19] 10.3880 9.7347 12.6420 12.5277 11.4129
X = Capacities;
x = nan(1,19);
for i=1:length(cell_array.cells)
    cell_array.cells(i).has_major_or
    if(cell_array.cells(i).has_major_or)
        ii = min([i+1 length(cell_array.cells)]);
    else
        ii = i;
    end
    ml_link = cell_array.cells(ii).segments(1).ml_link;
    ml_fd = ml_link.fd;
    x(i)=ml_fd.f_max.*ml_link.lanes*SimT*3600;
end
find(abs(X-x)>0.01)

%% FreeFlowSpeeds % [1x19] 0.2259 0.0807 0.2262 0.2481 0.1419
X = FreeFlowSpeeds;
x = nan(1,19);
for i=1:length(cell_array.cells)
    cell_array.cells(i).has_major_or
%     if(cell_array.cells(i).has_major_or)
%         ii = min([i+1 length(cell_array.cells)]);
%     else
        ii = i;
%     end
    ml_links = [cell_array.cells(ii).segments.ml_link];
    ml_fd = ml_links(1).fd;
    x(i)=(SimT*3600)*ml_fd.vf/sum([ml_links.length_meters]);
end
find(abs(X-x)>0.01)

%% CongestionSpeeds % [1x19]  0.0434 0.0209 0.0785 0.0891 0.0251
X = CongestionSpeeds;
x = nan(1,19);
for i=1:length(cell_array.cells)
    cell_array.cells(i).has_major_or
%     if(cell_array.cells(i).has_major_or)
%         ii = min([i+1 length(cell_array.cells)]);
%     else
        ii = i;
%     end
    ml_links = [cell_array.cells(ii).segments.ml_link];
    ml_fd = ml_links(1).fd;
    x(i)=(SimT*3600)*ml_fd.w/sum([ml_links.length_meters]);
end
find(abs(X-x)>0.01)

%% JamDensities % [1x19]  285.2591  585.5744  217.0184  191.1262  534.7995
X = JamDensities;
x = nan(1,19);
for i=1:length(cell_array.cells)
    cell_array.cells(i).has_major_or
    if(cell_array.cells(i).has_major_or)
        ii = min([i+1 length(cell_array.cells)]);
    else
        ii = i;
    end
    ml_links = [cell_array.cells(ii).segments.ml_link];
    ml_fd = ml_links(1).fd;
    lanes = ml_links(1).lanes;
    x(i)=ml_fd.dens_jam*sum([ml_links.length_meters])*lanes;
end
find(abs(X-x)>0.01)

%% SimT % [1x1] 0.0014
%% DownstreamCongested % [1x1] 0
%% boundedImputation % [1x1] 1
 
 
 