function [Y]=make_imputation_input()

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src','def'))

Y = struct( 'Flow',[],...    % [720x19]	[1.2588 0.9783 1.2917 1.0691 1.3340 1.4452 1.3160 1.6315 0.9296 1.0141 1.0296 1.0188 1.0584 0.8337 0.7375 1.1132 1.3530 0.6638 0.4490]
'Density',[],...	         % [720x19]	[6.1420 13.8526 5.6753 4.7285 10.5255 6.0089 12.4033 37.9275 25.1520 11.8565 9.2270 1.6126 14.6944 8.4561 3.0921 49.8209 12.8484 5.0735 0.7417]		
'Speed',[],...		         % [720x19] [60.1183 59.1879 72.5572 63.3952 63.2871 52.4068 61.2939 50.3617 65.4231 65.6263 67.3529 65.5381 55.1602 61.0061 64.3730 64.2935 68.1924 67.6857 60.3645]	
'OrPresent',[],...           % [1x19] [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0]   
'FrPresent',[],...		     % [1x19] [0 1 0 1 1 0 1 1 1 1 1 0 1 1 0 1 1 1 0]
'ImputeOR',[],...		     % [1x19] [0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 1 1 0 0]
'ImputeFR',[],...		     % [1x19] [0 0 0 0 0 0 1 1 1 1 1 0 1 1 0 1 0 1 0]
'CellLengths',[],...	     % [1x19] [0.4056 1.1590 0.4425 0.3898 0.6946 0.3019 0.8009 1.6196 2.2637 1.0671 0.8355 0.1437 1.0629 0.8577 0.3743 3.9895 1.8192 0.7184 0.1387
'ORBoundsLower',[],...       % [720x18] zeros
'ORBoundsUpper',[],...	     % [720x18] [0.0367 0.1026 0.0543 0.0169 0.0304 2.7778 0.0539 0.2077 0.0196 0.0125 0.0455 0.0113 0.0492 0.0369 0.1346 0.0034  0  0
'FRBoundsLower',[],...	     % [720x18] zeros
'FRBoundsUpper',[],...	     % [720x18] [0 0.0977  0 0.0943 0.0323  0 2.7778 5.5556 2.7778 2.7778 2.7778  0 2.7778 2.7778  0 0.1084 0.6536 2.7778
'Demand_Giv',[],...	         % [720x19] [0 0.0367 0.1026 0.0543 0.0169 0.0304  0 0.0539 0.2077 0.0196 0.0125 0.0455 0.0113 0.0492 0.0369 0.1346 0.0034  0  0
'DemandLower',[],...	     % [720x19] zeros
'DemandUpper',[],...	     % [720x19] 0
'BETA_Giv',[],...		     % [720x19] [0 0.0367 0.1026 0.0543 0.0169 0.0304 2.7778 0.0539 0.2077 0.0196 0.0125 0.0455 0.0113 0.0492 0.0369 0.1346 0.0034  0  0
'BETALower',[],...		     % [720x19] zeros
'BETAUpper',[],...		     % [720x19] [0 0.0874  0 0.0804 0.0216  0 1.0000 1.0000 1.0000 1.0000 1.0000  0 1.0000 1.0000  0 0.0899 0.9542 1.0000  0
'OrFlow_Giv',[],...	         % [720x19] [0.0190 0.0367 0.1026 0.0543 0.0169 0.0304  0 0.0539 0.2077 0.0196 0.0125 0.0455 0.0113 0.0492 0.0369 0.1346 0.0034  0  0
'FrFlow_Giv',[],...	         % [720x19] [0 0.0977  0 0.0943 0.0323  0  0  0  0  0  0  0  0  0  0 0.1084 0.6536  0  0
'Capacities',[],...	         % [1x19] [10.3880 9.7347   12.6420   12.5277   11.4129   11.3680   12.2990   12.0867   10.5840   10.5840 9.9307 9.9960   10.5677 9.4897 9.9960   12.7890   11.7600   11.7600 9.9307
'FreeFlowSpeeds',[],...      % [1x19] [0.2259 0.0807 0.2262 0.2481 0.1419 0.2883 0.1212 0.0561 0.0440 0.0933 0.1208 0.7239 0.0819 0.1079 0.2694 0.0242 0.0533 0.1350 0.7045
'CongestionSpeeds',[],...    % [1x19] [0.0434 0.0209 0.0785 0.0891 0.0251 0.1150 0.0434 0.0214 0.0127 0.0270 0.0261 0.1507 0.0165 0.0374 0.0928 0.0087 0.0191 0.0483 0.1002
'JamDensities',[],...	     % [1x19] 1000*[0.2853 0.5856 0.2170 0.1911 0.5348 0.1383 0.3851 0.7791 1.0738 0.5062 0.4629 0.0801 0.7689 0.3415 0.1449 1.9977 0.8368 0.3304 0.1132
'SimT',[],...                % [1x1] 0.0014
'DownstreamCongested',[],... % [1x1] 0
'boundedImputation',[]);     % [1x1] 1

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

%% OK Flow % [720x19] 1.2588 0.9783 1.2917 1.0691 1.3340
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
Y.Flow = x;
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
Y.Density = x;
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
Y.Speed = x;
clear X x

%% OK OrPresent % [1x19] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0
X = OrPresent;
x = nan(1,19);
for i=1:length(cell_array.cells)    
    xx = false;
    for j=1:length(cell_array.cells(i).segments)
        xx = xx || ~isempty(cell_array.cells(i).segments(j).get_or_dwnstr);  % NOTICE DOWNSTREAM!!
    end
    x(i) = xx;
end
if(any(X~=x))
    disp('failed: OrPresent')
end
Y.OrPresent = x;

%% OK FrPresent % [1x19] 0 1 0 1 1 0 1 1 1 1 1 0 1 1 0 1 1 1 0
X = FrPresent;
x = nan(1,19);
for i=1:length(cell_array.cells)
    x(i) = cell_array.cells(i).has_fr;
end
if(any(X~=x))
    disp('failed: FPresent')
end
Y.FrPresent = x;

%% OK ImputeOR % [1x19] 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 1 1 0 0
X = ImputeOR;
x = nan(1,19);
for i=1:length(cell_array.cells) 
    xx = false;
    for j=1:length(cell_array.cells(i).segments)
        or_dwnstr = cell_array.cells(i).segments(j).get_or_dwnstr;
        if(~isempty(or_dwnstr))
            xx = xx || isempty(or_dwnstr.sensor) || ~or_dwnstr.sensor.is_good_on_day(day,goodthresh);
        end
    end
    x(i) = xx;
end
if(any(X~=x))
    disp('failed ImputeOR')
end

%% OK ImputeFR % [1x19] 0 0 0 0 0 0 1 1 1 1 1 0 1 1 0 1 0 1 0
X = ImputeFR;
x = nan(1,19);
for i=1:length(cell_array.cells) 
    xx = false;
    for j=1:length(cell_array.cells(i).segments)
        fr = cell_array.cells(i).segments(j).fr_links;
        if(~isempty(fr))
            xx = xx || isempty(fr.sensor) || ~fr.sensor.is_good_on_day(day,goodthresh);
        end
    end
    x(i) = xx;
end
if(any(X~=x))
    disp('failed ImputeOR')
end

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
Y.CellLengths = x;

%% Demand_Giv % [720x19] 0 0.0367 0.1026 0.0543 0.0169
 % BETA_Giv % [720x19] 0 0.0874  0 0.0804 0.0216
Impute=ImputeFR|ImputeOR;
[xDemand,xBETA] = deriveDemandsAndSplitsFromRampFlows( Impute,...
                                                     Density, ...
                                                     OrFlow_Giv, ...
                                                     FrFlow_Giv, ...
                                                     get_time, ...
                                                     Capacities, ...
                                                     FreeFlowSpeeds, ...
                                                     CongestionSpeeds, ...
                                                     JamDensities);
                                                 
X = Demand_Giv
x = xDemand;


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
Y.OrFlow_Giv = x;
clear X x

%% FrFlow_Giv % [720x19]  0 0.0977  0 0.0943 0.0323

%% OK Capacities % [1x19] 10.3880 9.7347 12.6420 12.5277 11.4129
X = Capacities;
x = nan(1,19);
for i=1:length(cell_array.cells)
    cell_array.cells(i).has_major_or;
    if(cell_array.cells(i).has_major_or)
        ii = min([i+1 length(cell_array.cells)]);
    else
        ii = i;
    end
    ml_link = cell_array.cells(ii).segments(1).ml_link;
    ml_fd = ml_link.fd;
    x(i)=ml_fd.f_max.*ml_link.lanes*SimT*3600;
end
if(any(abs(X-x)>0.01))
    disp('failed Capacities')
end
Y.Capacities = x;

%% OK FreeFlowSpeeds % [1x19] 0.2259 0.0807 0.2262 0.2481 0.1419
X = FreeFlowSpeeds;
x = nan(1,19);
for i=1:length(cell_array.cells)
%     if(cell_array.cells(i).has_major_or)
%         ii = min([i+1 length(cell_array.cells)]);
%     else
        ii = i;
%     end
    ml_links = [cell_array.cells(ii).segments.ml_link];
    ml_fd = ml_links(1).fd;
    x(i)=(SimT*3600)*ml_fd.vf/sum([ml_links.length_meters]);
end
if(any(abs(X-x)>0.01))
    disp('failed FreeFlowSpeeds')
end
Y.FreeFlowSpeeds = x;

%% OK CongestionSpeeds % [1x19]  0.0434 0.0209 0.0785 0.0891 0.0251
X = CongestionSpeeds;
x = nan(1,19);
for i=1:length(cell_array.cells)
%     if(cell_array.cells(i).has_major_or)
%         ii = min([i+1 length(cell_array.cells)]);
%     else
        ii = i;
%     end
    ml_links = [cell_array.cells(ii).segments.ml_link];
    ml_fd = ml_links(1).fd;
    x(i)=(SimT*3600)*ml_fd.w/sum([ml_links.length_meters]);
end
if(any(abs(X-x)>0.01))
    disp('CongestionSpeeds')
end
Y.CongestionSpeeds = x;

%% OK JamDensities % [1x19]  285.2591  585.5744  217.0184  191.1262  534.7995
X = JamDensities;
x = nan(1,19);
for i=1:length(cell_array.cells)
    if(cell_array.cells(i).has_major_or)
        fd_ii = min([i+1 length(cell_array.cells)]);
    else
        fd_ii = i;
    end
    ml_fd    = cell_array.cells(fd_ii).segments(1).ml_link.fd;      % downstream
    ml_lanes = cell_array.cells(fd_ii).segments(1).ml_link.lanes;   % downstream
    ml_links = [cell_array.cells(i).segments.ml_link];              % not downstream
    x(i)=ml_fd.dens_jam*sum([ml_links.length_meters])*ml_lanes;
end
if(any(abs(X-x)>0.01))
    disp('failed JamDensities')
end
Y.JamDensities = x;

%% OK SimT % [1x1] 0.0014
Y.SimT = 5/3600;

%% OK DownstreamCongested % [1x1] 0
Y.DownstreamCongested = 0;

%% OK boundedImputation % [1x1] 1
Y.boundedImputation = 1;
