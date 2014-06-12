root = fileparts(fileparts(mfilename('fullpath')));

cfg.config.folder = fullfile(root,'cfg');
cfg.config.initial = fullfile(cfg.config.folder,'210W_initial.xml');
cfg.config.prefix = fullfile(cfg.config.folder,'210W_v');
cfg.data.fileshare = 'Z:\traffic_data\pems';
cfg.data.generated = fullfile(root,'generated');
cfg.fwy.district = '7';

cfg.days.daily = (datenum(2014,2,1):datenum(2014,2,7))';
cfg.days.hourly = (datenum(2014,2,1):datenum(2014,2,7))';
cfg.days.fivemin = (datenum(2014,2,1):datenum(2014,2,7))';
cfg.days.good = (datenum(2014,2,1):datenum(2014,2,7))';

cfg.fdcalibration.algorithm = 'simplefit';

cfg.boundaryflows.doHistorical=true;
cfg.boundaryflows.time = datenum(2014,2,3);
cfg.boundaryflows.predictionHorizon =  datenum(0,0,0,23,55,0);
cfg.boundaryflows.data.days =  datenum(2014,2,3);
cfg.boundaryflows.data.dtAgg = datenum(0, 0, 0, 0, 5, 0);
cfg.boundaryflows.data.doForcePemsToMat = false;
datafolder = fullfile(cfg.data.generated,'data');
cfg.boundaryflows.options.dataClean = fullfile(datafolder,'dataBF.mat');
cfg.boundaryflows.data.saveFilenameLanes = fullfile(datafolder,'BF_lanes.mat');
cfg.boundaryflows.data.processed_data_folder = datafolder;
cfg.boundaryflows.options.dirInput = datafolder;
cfg.boundaryflows.options.doPlot = false;

cfg.splitratios.data.district=7;
cfg.splitratios.data.historicalDays=datenum(2014,2,3);
cfg.splitratios.data.targetDays=datenum(2014,2,3);
cfg.splitratios.options.prediction=false;
cfg.splitratios.options.default_fr = 0.1;

cfg.imputation.day = cfg.days.fivemin(cfg.days.fivemin==datenum(2014,2,3));
cfg.imputation.goodthresh = 0.7;
cfg.imputation.boundarycondition = 1;  % 2 for congested downstream conditions (2 needs to be tested)
cfg.imputation.DemandTimeStep = 5/60; % 5 min default demand profile granularity
cfg.imputation.SimulationTimeStep = 5/60/60; % 5 sec default simulation time step
cfg.imputation.MeasurementTimeStep = 5/60; % 5 min measurements
cfg.imputation.FirstAndLastDetectorPostmiles = [42 23]; % make sure these correspond to healthy detectors, otherwise the result plots will be off by several miles, doesn't affect the results only affects plots.

% HOVinfo given as user input (this is only needed when the HOV lane detection is not a separate VDS)
cfg.imputation.HOVinfo.HOVLinks = [];
cfg.imputation.HOVinfo.HOVLanes = zeros(1,127);
cfg.imputation.HOVinfo.HOVLanes(cfg.imputation.HOVinfo.HOVLinks) = 1;
cfg.imputation.HOVinfo.Hours = [];
cfg.imputation.HOVinfo.PPV = 1.2;
cfg.imputation.HOVinfo.HOVexists = zeros(1,127);
cfg.imputation.HOVinfo.HOVexists(cfg.imputation.HOVinfo.HOVLinks) = 1;
cfg.imputation.HOVinfo.HOVseparate = 0;

% Major Ramps
cfg.imputation.MajorOR = [31 96];
cfg.imputation.MajorFR = [25 92];

% Ramp Flow Bounds (temporary measure)
cfg.imputation.Bounds = [];

% Data file suffix
cfg.imputation.suffix = 'test1';

% Artificial Mainline Data (temporary measure)
cfg.imputation.artificialDataFileName = [];
cfg.imputation.insertArtificialData = 0;

% overrides given as user input (eventually result of fault detection)
switch cfg.imputation.day
    case datenum(2014,2,3)
        cfg.imputation.overrides = [772918 717682 717673 761374 717657 717649 717642 717624];%772902];
    otherwise
        cfg.imputation.overrides = [];
end

% Optional imputation inputs like fudged mainline data, ramp flow
% bounds. Any other optional inputs that come up can be added root
cfg.imputation.optionalInputs = {cfg.imputation.Bounds, cfg.imputation.insertArtificialData, cfg.imputation.artificialDataFileName};


%% folder definitions
cfg.folders.pems_daily = fullfile(cfg.data.fileshare,'daily');
cfg.folders.pems_hourly = fullfile(cfg.data.fileshare,'hourly');
cfg.folders.pems_5min = fullfile(cfg.data.fileshare,'5min');
cfg.folders.processed = fullfile(cfg.data.generated,'data');
cfg.folders.reports = fullfile(cfg.data.generated,'reports');
cfg.folders.imputation = fullfile(cfg.data.generated,'imputation');

%% search path
wrkspace = fileparts(root);
addpath(genpath_nogit(fullfile(wrkspace,'L0-utilities')))
addpath(genpath_nogit(fullfile(wrkspace,'L0-imputation')))
addpath(genpath_nogit(fullfile(wrkspace,'L0-fd-calibration')))
addpath(genpath_nogit(fullfile(wrkspace,'L0-split-ratios')))
% addpath(fullfile(wrkspace,'L0-boundary-flows')); setupBoundaryFlows;
clear wrkspace

se = ScenarioPtr;
se.load([cfg.config.prefix '8.xml']);

% run_imputation_2(se,cfg.folders.processed,imputation_day)
se = run_imputation(se, ...
    cfg.folders.processed, ...
    cfg.imputation.day, ...
    cfg.imputation.overrides, ...
    cfg.imputation.goodthresh, ...
    cfg.imputation.boundarycondition, ...
    cfg.imputation.HOVinfo, ...
    cfg.imputation.MajorOR, ...
    cfg.imputation.MajorFR,...
    cfg.folders.imputation,...
    cfg.imputation.suffix,...
    cfg.imputation.MeasurementTimeStep,...
    cfg.imputation.SimulationTimeStep,...
    cfg.imputation.DemandTimeStep,...
    cfg.imputation.FirstAndLastDetectorPostmiles,...
    cfg.imputation.optionalInputs);

se.save([cfg.config.prefix '10.xml']);

%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post processing for vehicle types and omitted ramps
% addNominalDemandsAndSplitRatiosToOmittedRamps(se); % also assign first and last fundamental diagrams to snipped links upstream and downstream
% AdjustVehicleTypeDemandsAndSplits(se, cfg.imputation.HOVinfo);
%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

se.save([cfg.config.prefix '11.xml']);
% NOTE: Need to finish this.