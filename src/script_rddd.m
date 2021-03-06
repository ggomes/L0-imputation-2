clear
close all

root = fileparts(fileparts(mfilename('fullpath')));

load(fullfile(root,'analysis','imputation_input_test1_03-Feb-2014'))

X.Flow = Flow;    % [720x19]	[1.2588 0.9783 1.2917 1.0691 1.3340 1.4452 1.3160 1.6315 0.9296 1.0141 1.0296 1.0188 1.0584 0.8337 0.7375 1.1132 1.3530 0.6638 0.4490]
X.Density = Density;	         % [720x19]	[6.1420 13.8526 5.6753 4.7285 10.5255 6.0089 12.4033 37.9275 25.1520 11.8565 9.2270 1.6126 14.6944 8.4561 3.0921 49.8209 12.8484 5.0735 0.7417]		
X.Speed = Speed;		         % [720x19] [60.1183 59.1879 72.5572 63.3952 63.2871 52.4068 61.2939 50.3617 65.4231 65.6263 67.3529 65.5381 55.1602 61.0061 64.3730 64.2935 68.1924 67.6857 60.3645]	
X.OrPresent = OrPresent;           % [1x19] [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0]   
X.FrPresent = FrPresent;		     % [1x19] [0 1 0 1 1 0 1 1 1 1 1 0 1 1 0 1 1 1 0]
X.ImputeOR = ImputeOR;		     % [1x19] [0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 1 1 0 0]
X.ImputeFR = ImputeFR;		     % [1x19] [0 0 0 0 0 0 1 1 1 1 1 0 1 1 0 1 0 1 0]
X.CellLengths = CellLengths;	     % [1x19] [0.4056 1.1590 0.4425 0.3898 0.6946 0.3019 0.8009 1.6196 2.2637 1.0671 0.8355 0.1437 1.0629 0.8577 0.3743 3.9895 1.8192 0.7184 0.1387
X.ORBoundsLower = ORBoundsLower;       % [720x18] zeros
X.ORBoundsUpper = ORBoundsUpper;	     % [720x18] [0.0367 0.1026 0.0543 0.0169 0.0304 2.7778 0.0539 0.2077 0.0196 0.0125 0.0455 0.0113 0.0492 0.0369 0.1346 0.0034  0  0
X.FRBoundsLower = FRBoundsLower;	     % [720x18] zeros
X.FRBoundsUpper = FRBoundsUpper;	     % [720x18] [0 0.0977  0 0.0943 0.0323  0 2.7778 5.5556 2.7778 2.7778 2.7778  0 2.7778 2.7778  0 0.1084 0.6536 2.7778
X.Demand_Giv = Demand_Giv;	         % [720x19] [0 0.0367 0.1026 0.0543 0.0169 0.0304  0 0.0539 0.2077 0.0196 0.0125 0.0455 0.0113 0.0492 0.0369 0.1346 0.0034  0  0
X.DemandLower = DemandLower;	     % [720x19] zeros
X.DemandUpper = DemandUpper;	     % [720x19] 0
X.BETA_Giv = BETA_Giv;		     % [720x19] [0 0.0367 0.1026 0.0543 0.0169 0.0304 2.7778 0.0539 0.2077 0.0196 0.0125 0.0455 0.0113 0.0492 0.0369 0.1346 0.0034  0  0
X.BETALower = BETALower;		     % [720x19] zeros
X.BETAUpper = BETAUpper;		     % [720x19] [0 0.0874  0 0.0804 0.0216  0 1.0000 1.0000 1.0000 1.0000 1.0000  0 1.0000 1.0000  0 0.0899 0.9542 1.0000  0
X.OrFlow_Giv = OrFlow_Giv;	         % [720x19] [0.0190 0.0367 0.1026 0.0543 0.0169 0.0304  0 0.0539 0.2077 0.0196 0.0125 0.0455 0.0113 0.0492 0.0369 0.1346 0.0034  0  0
X.FrFlow_Giv = FrFlow_Giv;	         % [720x19] [0 0.0977  0 0.0943 0.0323  0  0  0  0  0  0  0  0  0  0 0.1084 0.6536  0  0
X.Capacities = Capacities;	         % [1x19] [10.3880 9.7347   12.6420   12.5277   11.4129   11.3680   12.2990   12.0867   10.5840   10.5840 9.9307 9.9960   10.5677 9.4897 9.9960   12.7890   11.7600   11.7600 9.9307
X.FreeFlowSpeeds = FreeFlowSpeeds;      % [1x19] [0.2259 0.0807 0.2262 0.2481 0.1419 0.2883 0.1212 0.0561 0.0440 0.0933 0.1208 0.7239 0.0819 0.1079 0.2694 0.0242 0.0533 0.1350 0.7045
X.CongestionSpeeds = CongestionSpeeds;    % [1x19] [0.0434 0.0209 0.0785 0.0891 0.0251 0.1150 0.0434 0.0214 0.0127 0.0270 0.0261 0.1507 0.0165 0.0374 0.0928 0.0087 0.0191 0.0483 0.1002
X.JamDensities = JamDensities;	     % [1x19] 1000*[0.2853 0.5856 0.2170 0.1911 0.5348 0.1383 0.3851 0.7791 1.0738 0.5062 0.4629 0.0801 0.7689 0.3415 0.1449 1.9977 0.8368 0.3304 0.1132
X.SimT = SimT;                % [1x1] 0.0014
X.DownstreamCongested = DownstreamCongested; % [1x1] 0
X.boundedImputation = boundedImputation;     % [1x1] 1

clear Flow Density Speed *Present *Bounds* Demand* BETA* Capacities
clear *Speeds JamDensities SimT DownstreamCongested bounded*

old_Imputation_Core_Algorithm(X);
                                     