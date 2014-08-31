function [FlowCompare,...
          Nh,...
          Velocity,...
          OrFlow,...
          FrFlow,...
          ORINP,...
          BETAF,...
          DJ,...
          c,...
          ImputeOR,...
          ImputeFR,...
          OrPresent,...
          FrPresent,...
          OrFlow_Giv,...
          FrFlow_Giv,...
          DensityError,...
          FlowError] = old_Imputation_Core_Algorithm(X)

%IMPUTATION_CORE_ALGORITHM % Takes in the data and topology matrices and
%produces Demand, Split Ratio, Onramp Flow, Offramp flow profiles along
%with simulated mainline flows and densities, which are inputs to the
%mega-cell splitting algorithm later on
%
%
%   INPUTS:
%
%   (For the following, n is the number of cells in the network and t is
%   the number of timesteps for the given day)
%   1) Flow: txn matrix of flow measurement values. The MeasureT time-step
%   is upsampled to the SimT granularity by linear interpolation
% 
%   2) Density: txn matrix of density measurement values.
% 
%   3) ActVelocity: txn matrix of speed measurement values.
% 
%   4) OrPresent: 1xn vector of 1s and 0s, specifying for each cell whether
%   the cell has an onramp upstream or not.
% 
%   5) FrPresent: 1xn vector of 1s and 0s, specifying for each cell whether
%   the cell has an offramp downstream or not.
% 
%   6) ImputeOR: 1xn vector of 1s and 0s, specifying if for each cell
%   whether its onramp needs to be imputed.
% 
%   7) ImputeFR: 1xn vector of 1s and 0s, specifying if for each cell
%   whether its offramp needs to be imputed.
% 
%   8) CellLengths: 1xn vector of cell lengths.
% 
%   9) ORBoundsLower: tx(n-1) matrix of onramp flow lower bounds. The very 
%   first cell is assumed to have no onramps and is therefore clipped.
% 
%   10) ORBoundsUpper: tx(n-1) matrix of onramp flow upper bounds. The very
%   first cell is assumed to have no onramps and is therefore clipped.
% 
%   11) FRBoundsLower: tx(n-1) matrix of offramp flow lower bounds. The 
%   very last cell is assumed to have no onramps and is therefore clipped.
% 
%   12) FRBoundsUpper: tx(n-1) matrix of offramp flow upper bounds. The 
%   very last cell is assumed to have no onramps and is therefore clipped.
% 
%   13) Demand: txn matrix of known demands, initialized to zero if none
%   available.
% 
%   14) DemandLower: txn matrix of demand lower bounds. The first one is
%   not clipped here because if there is known upstream onramp flows, they
%   are added to the upstream mainline flow and used as the upstream
%   boundary flow.
% 
%   15) DemandUpper: txn matrix of demand upper bounds.
% 
%   16) BETA: txn matrix of known split ratios, initialized to zero if none
%   available.
% 
%   17) BETALower: txn matrix of split ratio lower bounds. The last one is
%   not clipped but not used either since the very last offramp is assumed
%   unbounded.
% 
%   18) BETAUpper: txn matrix of split ratio upper bounds.
% 
%   19) OrFlow_Giv: txn matrix of known onramp flows. The very first onramp
%   flow is added to the upstream mainline flow and the total is used as
%   the upstream boundary flow in imputation and simulation.
% 
%   20) FrFlow_Giv: txn matrix of known offramp flows. The very last
%   offramp is not clipped but not used either.
% 
%   21) Qmax: 1xn vector of total capacities in [veh]. The [veh/h] value is
%   scaled by the simulation time-step and total number of lanes to obtain
%   a value in [veh/timestep]. Any additional dedicated lanes are 
%   aggregated here. I.e. if there is an HOV lane and 4 mainline lanes, the
%   capacity is given by: (4*mainlineCapacity + 1*HOVLaneCapacity)*SimT
% 
%   22) vf: 1xn vector of free flow speeds in [veh/timestep]
% 
%   23) w: 1xn vector of congestion wave speeds in [veh/timestep]
% 
%   24) rhojam: 1xn vector of jam densities in [vehicles]
% 
%   25) SimT: Simulation time-step in hours. 
% 
%   26) DownstreamCongested: a logical specifying downstream congestion
%   with the value 1 and free flow with the value 0
%
%   27) boundedImputation: 1 if there exist ramp bounds, 0 otherwise
%
%   OUTPUTS:
%
%   1) FlowCompare: txn matrix of simulated mainline flows corresponding to
%   assumed detector locations (i.e. flows entering each cell on the
%   mainline, excluding the onramp
% 
%   2) Nh: txn matrix of simulated densities (Nh is for n^hat)
% 
%   3) Velocity: txn matrix of simulated speeds
% 
%   4) OrFlow: txn matrix of simulated onramp flows
% 
%   5) FrFlow: txn matrix of simulated offramp flows
% 
%   6) ORINP: txn matrix of estimated freeway demands (i.e. from onramps to
%   the freeway mainline. Zeros column for nodes that don't have onramps.
% 
%   7) BETAF: txn matrix of estimated split ratios. Zeros column for nodes 
%   that don't have offramps.
% 
%   8) DJ: txn matrix of estimated onramp demands (i.e. from side streets
%   to onramps)
% 
%   9) c: txn matrix of estimated cumulative demands into each cell.
% 
%   10) ImputeOR: 1x(n-1) vector of indicators for active imputation at each
%   node, starting with the 1st node at the end of the first cell. Hence
%   the size is one shorter than the input
% 
%   11) ImputeFR: 1x(n-1) vector of indicators for active offramp imputation at
%   each node, excluding the very last node. Hence the shorter vector
% 
%   12) OrPresent: 1xn vector of 1s and 0s, specifying for each cell whether
%   the cell has an onramp upstream or not.
% 
%   13) FrPresent: 1xn vector of 1s and 0s, specifying for each cell whether
%   the cell has an offramp downsstream or not.
% 
%   14) OrFlow_Giv: shifted index for input to megacell splitting 
%           
%   15) FrFlow_Giv: shifted index for input to megacell splitting 
%
%   16) DensityError: Average percentage error in density measurements.
%   Calculated by taking the average relative difference between
%   measurements and simulation at each sensor location throughout the 
%   whole dayand then taking the average of these values over all sensors.
%
%   17) FlowError: Average percentage error in flow measurements.
%   Calculation same as Density Error.
%
%   COMMENTS:
%
%       + In the LN-CTM imputation, the ramps are indexed by node rather
%       than by cell. Hence, node 1, which is the node at the end of
%       the first cell, has the FIRST cell's offramp and the SECOND cell's
%       onramp attached to it. Therefore the onramp indices need to be
%       shifted up by one. This is done in the pre-processing portion.
%       
%       + Since the very last and very first nodes do not need imputation, 
%       the ImputeOR, ImputeFR, OrPresent, FrPresent, OrFlow_Giv and 
%       FrFlow_Giv variables are snipped by one entry. Hence, the 
%       imputation loops over n-1 nodes for imputation. To keep the data
%       matrix sizes consistent, zero entries and columns are added for
%       omitted ramps at the very end.




Flow = X.Flow;
Density = X.Density;
ActVelocity = X.Speed;
OrPresent = X.OrPresent;
FrPresent = X.FrPresent;
ImputeOR = X.ImputeOR;
ImputeFR = X.ImputeFR;
CellLengths = X.CellLengths;
ORBoundsLower = X.ORBoundsLower;
ORBoundsUpper = X.ORBoundsUpper;
FRBoundsLower = X.FRBoundsLower;
FRBoundsUpper = X.FRBoundsUpper;
Demand = X.Demand_Giv;
DemandLower = X.DemandLower;
DemandUpper = X.DemandUpper;
BETA = X.BETA_Giv;
BETALower = X.BETALower;
BETAUpper = X.BETAUpper;
OrFlow_Giv = X.OrFlow_Giv;
FrFlow_Giv = X.FrFlow_Giv;
Qmax = X.Capacities;
vf = X.FreeFlowSpeeds;
w = X.CongestionSpeeds;
rhojam = X.JamDensities;
SimT = X.SimT;
DownstreamCongested = X.DownstreamCongested;
boundedImputation = X.boundedImputation;

%% Checks and Pre-processing
if FrPresent(end) && ImputeFR(end)
    display('+ The last cell has offramps which need to be imputed: Zero flows are assumed for this offramp')
end 

if DownstreamCongested

    BoundaryVel=ActVelocity(:,end);   
    Boundaryvf=vf(end)/SimT*CellLengths(end);    
    
else
    
    BoundaryVel=ActVelocity(:,end);  

end

% Check CFL conditions
if any(vf>1)    
    error('TOPL:User','Some of the sections are too short: Increase section lengths')
end

numtime = size(Flow,1);
numcell = size(Flow,2);
    
InputFlw = Flow(:,1) + OrFlow_Giv(:,1);  

Flow(:,1) = InputFlw;  % This is for comparison purposes, i.e. we add any onramp flow to the mainline if there is an onramp flow in the most upstream cell

% Note that the n'th cell offramp is paired with the n+1'th cell onramp for
% imputation. Also, we do not impute the first cell onramp and the last
% cell offramp. Hence we shift up the indices of onramps by using only 2nd
% onramp to last onramp and we snip the last offramp. Note that we end up
% with ramps indexed by nodes and of length one less than the number of
% cells.

ImputeOR=ImputeOR(2:end);  
ImputeFR=ImputeFR(1:end-1);
OrPresent=OrPresent(2:end);  
FrPresent=FrPresent(1:end-1);

OrFlow = OrFlow_Giv(:,2:end);
FrFlow = FrFlow_Giv(:,1:end-1);

OrFlow_Giv = OrFlow;
FrFlow_Giv = FrFlow;

Impute=ImputeFR|ImputeOR;
NoNodes=size(Density,2)+1;

%% Check that ramps are present whenever Imputation is enabled
% Specify whether onramps / offramps are only present etc.

LBounds=zeros(NoNodes-2,1);
UBounds=zeros(NoNodes-2,1);
for Ind=1:NoNodes-2
    LBounds(Ind)=~FrPresent(Ind);
    UBounds(Ind)=~OrPresent(Ind);
    if LBounds(Ind) && UBounds(Ind)
        if Impute(Ind)
            error('Both onramp and offramp not present, and Imputation is enabled')
        end
    end
end
if DownstreamCongested
    LBounds(NoNodes-1)=1;
    UBounds(NoNodes-1)=0;
end

%% LEARNING PORTION 

%Set Learning Parameters
% Parameter and variable names are consistent with the corresponding
% publications
GM=40;
G1=1*GM;
G2=0.001*GM;
PercTol=0.001;
PercTol2=0.2;
MAXLIM=100;
ITERMAX=25;
IterTrigger=[5 9 12 16 20]; %Trigger initiated after this Iteration.

% Change of c values compared to previous values is bounded by DerivativeBound number of vehicles per SimT seconds.. Choose this value with care..
DerivativeBound=2;
StartBound=1;

clear c cj;
c = zeros(numtime,numcell-1);
% THIS CODE WAS ELIMINATED BY R. HOROWITZ - PLEASE CHECK WHY IT WAS HERE
% for j=2:NoNodes-1
%     cj=1.5*w(j)*(rhojam(j)-Density(:,j));%+randn(size(STime))'*50+500;
%     c(:,j-1)=cj;
% end
% WE NEED TO CHECK IF THIS IS CORRECT
if DownstreamCongested
    c(:,NoNodes-1)=Qmax(NoNodes-1);
    Impute(NoNodes-1)=1;
end
% clear cj

disp('      o starting imputation loop')
cbest=c;
if boundedImputation
    IterBound_D_Beta = 2;
else
    IterBound_D_Beta = ITERMAX+1;  % CURRENTLY DISABLED
end
StartIterBound = 10 ; % make startIterBound greater than Iterbound_D_Beta

% Outer Loop -Perform Multiple iterations for convergence
% Middle Loop -Run the learning algorithm over 24 hours
% Innermost loop - Cycle through various links and learn the parameters corresponding to that link
%
% First link is in freeflow.. also the sink is assumed to be in freeflow.
% However, first link is indicated to be in congestion... This is just to
% decrease coding complexity.
% cj_LowerMatrix = zeros(size(c));
% cj_UpperMatrix = zeros(size(c));
for Iter=1:ITERMAX;  % maximum iterations to converge..
    
    Flags(1)=1;
    
    Mode=c(:,1:end-1)*0;
    
    InQ=0;
    Nh=zeros(numtime,numcell);
    Nh(1,:)=Density(1,:);
    cj=c(1,:);
    
    for i=1:length(Flow)-1
        
        Limit=min(Qmax,w.*(rhojam-Nh(i,:)));
        cjprev=cj;
        cj=c(i,:);
        cj=max([cjprev-DerivativeBound;min([cj;cjprev+DerivativeBound])]);

        for j = 1 : NoNodes-2
            
            NjVfj=min(Qmax(j),Nh(i,j)*vf(j));
            
            if ~Impute(j)
                cj(j)=NjVfj*(1-BETA(i,j))+Demand(i,j);
            else
                if Iter>StartIterBound                    
                    if LBounds(j)
                        cj(j)= max(cj(j),NjVfj);
                    end
                    if UBounds(j)
                        cj(j)= min(cj(j),NjVfj);
                    end
                end
            end
            
        end
        
        for j=2:NoNodes-1;
            if Impute(j-1)
                Flags(j)=cj(j-1)>=Limit(j)*(1-PercTol);
            else
                Flags(j)=cj(j-1)>=Limit(j);
            end
        end
        Flags(NoNodes)=0;
        
        InQ=InQ+InputFlw(i);
        Inflow=min(Limit(1),InQ);
        InQ=InQ-Inflow(1);
        
        j=0;
        
        NeedBoundaryImpute=(DownstreamCongested && BoundaryVel(i)<Boundaryvf * 0.9); % Arbitrary velocity for activating boundary imputation.
        if NeedBoundaryImpute
            Flags(NoNodes)=1;
        end
            
        
        while j<NoNodes-1
            j=j+1;
            NjVfj=min(Nh(i,j)*vf(j),Qmax(j));
            if j>1
            Njm1Vfjm1=min(Qmax(j-1),Nh(i,j-1)*vf(j-1));
            end
            Wj=min(Qmax(j),w(j)*(rhojam(j)-Nh(i,j)));
            if j==NoNodes-1
                Wjp1=Qmax(j);
            else
                Wjp1=min(Qmax(j+1),w(j+1)*(rhojam(j+1)-Nh(i,j+1)));
            end
           
            if ~Flags(j) && ~Flags(j+1)
                Mode(i,j)=1;
                
                BoundL=LBounds(j-1);
                BoundU=UBounds(j-1);              



                
                Nh0=Density(i+1,j)-(Nh(i,j)-NjVfj+cj(j-1)); %For second link
                NHt=Nh0/(1+G1);
                if Impute(j-1)
                    cj(j-1)=cj(j-1)+G1*NHt;
                    
                    if Iter > IterBound_D_Beta
                        UBound=min(Wj,Njm1Vfjm1*(1-BETALower(i,j-1))+DemandUpper(i,j));
                        LBound=max(0.01,min(Wj,Njm1Vfjm1*(1-BETAUpper(i,j-1))+DemandLower(i,j)));
                    else
                        UBound = Wj;
                        LBound = 0.01;
                    end
                    
                    if BoundU==1
                        UBound=min(UBound,Njm1Vfjm1);
                    elseif BoundL==1
                        LBound=max(LBound,Njm1Vfjm1);
                    end
                    if i>StartBound
                        LBound = max(LBound,cjprev(j-1)-DerivativeBound);
                        UBound = min(UBound,cjprev(j-1)+DerivativeBound);
                    end
                    
                    cj(j-1) = max(min(cj(j-1),UBound),LBound);
                end
                Nh(i+1,j)=Nh(i,j)-NjVfj+cj(j-1);
                
            elseif ~Flags(j) && Flags(j+1)
                Mode(i,j)=2;
                               
                % Check for selecting mode
                if Impute(j) && cj(j)>Wjp1*(1-PercTol) && cj(j)<Wjp1*(1+PercTol) && (Nh(i,j)+cj(j-1)-NjVfj)>Density(i+1,j)
                    cj(j)=min(cj(j),Wjp1*(1-PercTol/100));
                    Flags(j+1)=0;
                    j=j-1;
                    continue;
                end
                cj(j)=max(cj(j),Wjp1);
                
                BoundL1=LBounds(j-1);BoundU1=UBounds(j-1);
                BoundL2=LBounds(j);BoundU2=UBounds(j);
                
                
                Nh0=Density(i+1,j)-(Nh(i,j)+cj(j-1)-Wjp1*NjVfj/cj(j));
                NHt=Nh0/(1+G1*Impute(j-1)+G2*Wjp1*NjVfj*Impute(j));
                
                % Adapt second term
                if Impute(j)
                    cj(j)=1/(1/cj(j)-min(G2*NHt,1/cj(j)*0.99));
                    
                    
                    if BoundU2
                        UBound2=max(min(MAXLIM*Wjp1,NjVfj),0.01);
                    else
                        if Iter > IterBound_D_Beta
                            UBound2=min(MAXLIM*Wjp1,max(Wjp1,NjVfj*(1-BETALower(i,j))+DemandUpper(i,j)));
                        else
                            UBound2=MAXLIM*Wjp1;
                        end
                    end
                    
                    LBound2=max(NjVfj.*BoundL2,Wjp1);
                    if Iter > IterBound_D_Beta
                        LBound2=max(LBound2,NjVfj*(1-BETAUpper(i,j))+DemandLower(i,j));
                    else
                        LBound2=max(LBound2,0.01);
                    end                                  
                    
                    if i>StartBound
                        LBound2 = max(LBound2,cjprev(j)-DerivativeBound);
                        UBound2 = min(UBound2,cjprev(j)+DerivativeBound);
                    end
                                       
                    cj(j) = max(min(cj(j),UBound2),LBound2);
                end
                % Adapt first term
                if Impute(j-1)
                    cj(j-1)=cj(j-1)+G1*NHt; 
                    
                    if Iter > IterBound_D_Beta
                        UBound1=min(Wj,Njm1Vfjm1*(1-BETALower(i,j-1))+DemandUpper(i,j));
                        LBound1=max(0.01,min(Wj,Njm1Vfjm1*(1-BETAUpper(i,j-1))+DemandLower(i,j)));
                    else
                        UBound1 = Wj;
                        LBound1 = 0.01;
                    end

                    if BoundU1==1
                        UBound1=Njm1Vfjm1;
                    elseif BoundL1==1
                        LBound1=max(LBound1,Njm1Vfjm1);
                    end
                    
                    if i>StartBound
                        LBound1 = max(LBound1,cjprev(j-1)-DerivativeBound);
                        UBound1 = min(UBound1,cjprev(j-1)+DerivativeBound);
                    end
                    
                    cj(j-1) = max(min(cj(j-1),UBound1),LBound1);
                end
                Nh(i+1,j)=Nh(i,j)+cj(j-1)-Wjp1*NjVfj/cj(j);
                
            elseif Flags(j) && ~Flags(j+1)
                
                Mode(i,j)=3;
                if j==1
                    Nh(i+1,j)=Nh(i,j)+Inflow-min(Nh(i,j)*vf(j),Qmax(j));
                else
                    Nh(i+1,j)=Nh(i,j)+min(Qmax(j),w(j)*(rhojam(j)-Nh(i,j)))-min(Nh(i,j)*vf(j),Qmax(j));
                end
                
            elseif Flags(j) && Flags(j+1)
                
                Mode(i,j)=4;
                if Impute(j)
                    BoundL=LBounds(j);
                    BoundU=UBounds(j);
                    
                    if j==1
                        Nh0=Density(i+1,j)-(Nh(i,j)+Inflow-Wjp1*NjVfj/cj(j));
                    else
                        Nh0=Density(i+1,j)-(Nh(i,j)+Wj-Wjp1*NjVfj/cj(j));
                    end;
                    
                    % Check for selecting mode
                    if  cj(j)>Wjp1*(1-PercTol) && cj(j)<Wjp1*(1+PercTol) && Nh0<0
                        cj(j)=min(cj(j),Wjp1*(1-PercTol/100));
                        Flags(j+1)=0;
                        j=j-1;
                        continue;
                    end
                    cj(j)=max(cj(j),Wjp1);
                    % Otherwise continue
                    
                    NHt=Nh0/(1+G2*Wjp1*NjVfj);
                    
                    cj(j)=1/(1/cj(j)-min(G2*NHt,1/cj(j)*0.99));
                    
                    
                    if BoundU
                        UBound=max(min(MAXLIM*Wjp1,NjVfj),0.01);
                    else
                        if Iter > IterBound_D_Beta
                            UBound=min(MAXLIM*Wjp1,max(Wjp1,NjVfj*(1-BETALower(i,j))+DemandUpper(i,j)));
                        else
                            UBound=MAXLIM*Wjp1;
                        end
                    end
                    
                    LBound=max(NjVfj.*BoundL,Wjp1);
                    if Iter > IterBound_D_Beta
                        LBound=max(LBound,NjVfj*(1-BETAUpper(i,j))+DemandLower(i,j));
                    else
                        LBound=max(LBound,0.01);
                    end                                          
                    
                    if i>StartBound
                        LBound = max(LBound,cjprev(j)-DerivativeBound);
                        UBound = min(UBound,cjprev(j)+DerivativeBound);
                    end
                    
                    cj(j) = max(min(cj(j),UBound),LBound);
                    
                end
                if j==1
                    Nh(i+1,1)=Nh(i,j)+Inflow-Wjp1*NjVfj/cj(j);
                else
                    Nh(i+1,j)=Nh(i,j)+Wj-Wjp1*NjVfj/cj(j);
                end
            end
            %             clear UBound UBound1 UBound2 LBound LBound1 LBound2
            
        end
        cj=max(0.01,cj);
         
        c(i,:)=cj;
    end
    
    if any(any(Mode(1:end-1,:)==0))
        error('Check')
    end
    
    MError(Iter)=mean(mean(abs(Nh-Density)))./mean(mean(Density));
    
    disp(['        Error (' num2str(Iter) ':' num2str(ITERMAX) '): ' num2str(MError(Iter)*100,'%.2f') '%' ]);
    
    if Iter == StartIterBound || (Iter>StartIterBound && any(MError(Iter)<=MError(1:Iter-1)))
        cbest=c;
    end
    
    if any(Iter==IterTrigger)
        
        %  Find regions having bottlenecks  --- Mode is 3 in this case
        if DownstreamCongested
            ModeFilter=(Mode==3) | (Mode==1).*(ones(size(Mode,1),1)*[Impute]==0).*(Density-Nh>0) | (Mode==4).*(ones(size(Mode,1),1)*[Impute]==0).*(Density-Nh<0);
        else
            ModeFilter=(Mode==3) | (Mode==1).*(ones(size(Mode,1),1)*[0 Impute]==0).*(Density-Nh>0) | (Mode==4).*(ones(size(Mode,1),1)*[Impute 0]==0).*(Density-Nh<0);
        end
        BottleError=ModeFilter.*(Density-Nh)./max(Density,0.001).*(Density>0.001);
        
        %%You are trying to increase the ValCj, so.. only worry about the
        %%upperbound
        Filter=Mode*0+1;
        Filter(:,end)=0;
        [Row,Col]=find(BottleError.*Filter>0.05);
        
        for i=1:length(Row)
            if Impute(Col(i))
                Limit=min(Qmax(Col(i)+1),w(Col(i)+1)*(rhojam(Col(i)+1)-Nh(Row(i),Col(i)+1)))*(1+PercTol2); % Safety feature...
                c(Row(i),Col(i))=Limit;
            end
        end
        
        [Row,Col]=find(BottleError<-0.05);
        %%You are trying to reduce the ValCj, so.. only worry about the
        %%lowerbound
        for i=1:length(Row)
            if Col(i)>1 && Impute(Col(i)-1)
                Limit=min(Qmax(Col(i)),w(Col(i))*(rhojam(Col(i))-Nh(Row(i),Col(i))))*(1-PercTol2); % Safety feature...
                c(Row(i),Col(i)-1)=Limit;
            end
        end
        
    end
    
end
c=cbest;

%% LINEAR PROGRAM PORTION
% Determine the parameters dj and Beta from the estimated c's (cumulative demands)

Override=0;
Nha=Nh;
clear Nh
% Estimating parameters
BETA1=zeros(numtime,numcell);
dj1=zeros(numtime,numcell);
OrInp=zeros(numtime,numcell);
OnrampFlw=zeros(numtime,numcell);
FlowBet=Flow(:,2:end);
dprev=zeros(NoNodes-2,1);

%% ===========================LP WITH BOUNDS===============================
for i=2:length(Flow)
    dprev=dprev*0;
    for j=1:NoNodes-2;
        if Impute(j)
            Term1=min(Qmax(j+1),w(j+1)*(rhojam(j+1)-Nha(i,j+1)));
            Term2= min(Qmax(j),Nha(i,j)*vf(j));
            if ~OrPresent(j)
                 if all(FRBoundsUpper(:,j)==0)
                    BETA1(i,j)=max(0,1-c(i,j)/Term2); %(NO BOUNDS)
                 else
                     BETA1(i,j) = min(BETAUpper(i,j),max(BETALower(i,j),1-c(i,j)/Term2)); % BOUNDS APPLIED
                 end
                OnrampFlw(i,j)=0;
                continue;
            end
            
            if c(i,j)>Term1; %(CONGESTION)
                
                if ~FrPresent(j)
                     if all(ORBoundsUpper(:,j)==0)
                        dj1(i,j)=c(i,j)-Term2;  %(NO BOUNDS)
                     else                    
                         dj1(i,j) = min(DemandUpper(i,j+1),max(DemandLower(i,j+1),c(i,j)-Term2)); % BOUNDS APPLIED
                     end
                    OrInp(i,j)=dj1(i,j)-dprev(j);
                    OnrampFlw(i,j)=dj1(i,j)*Term1/c(i,j); 
                    dprev(j)=dj1(i,j)-OnrampFlw(i,j);
                    continue;
                end
                
                Term1=min(Qmax(j+1),w(j+1)*(rhojam(j+1)-Nha(i,j+1)));
                Term2= min(Qmax(j),Nha(i,j)*vf(j));
                fac = Term1/c(i,j);
                % Explicit solution
                if ImputeOR(j) && ImputeFR(j) || Override % both onramp and offramp need to be imputed
                    
                    if all(ORBoundsUpper(:,j)==0)
                        d_soln = min(c(i,j),max([0,dprev(j),Term1 - FlowBet(i,j),c(i,j)-Term2])); %(NO BOUNDS)
                    else
                        d_soln = min(DemandUpper(i,j+1),max([DemandLower(i,j+1),dprev(j),Term1 - FlowBet(i,j),c(i,j)-Term2])); % BOUNDS APPLIED
                    end
                    
                     if all(FRBoundsUpper(:,j)==0)
                        fr_soln = Term2 + d_soln - c(i,j); %(NO BOUNDS)
                     else
                         fr_soln = min(FRBoundsUpper(i,j),max(FRBoundsLower(i,j),Term2 + d_soln - c(i,j))); % BOUNDS APPLIED
                     end
                    
                elseif ~ImputeOR(j) % onramp is known, offramp needs to be imputed
                    
                    MultFac = abs(FlowBet(i,j)/max(100*SimT,OrFlow_Giv(i,j)))/2;  %Higher weighting for ramp flows.
                    if MultFac<1
                        d_vertex = (Term1 - FlowBet(i,j))/fac;
                    else
                        d_vertex = OrFlow_Giv(i,j)/fac;
                    end
                    
                    if all(ORBoundsUpper(:,j)==0)
                        d_soln = min(c(i,j),max([0,dprev(j),d_vertex,c(i,j)-Term2]));
                    else
                        d_soln = min(DemandUpper(i,j+1),max([DemandLower(i,j+1),dprev(j),d_vertex,c(i,j)-Term2])); % BOUNDS APPLIED
                    end
                    
                     if all(FRBoundsUpper(:,j)==0)
                        fr_soln = Term2 + d_soln - c(i,j);
                     else
                         fr_soln = min(FRBoundsUpper(i,j),max(FRBoundsLower(i,j),Term2 + d_soln - c(i,j))); % BOUNDS APPLIED
                     end
                    
                else % offramp is known, onramp needs to be imputed
                    
                    MultFac = abs(FlowBet(i,j)/max(100*SimT,FrFlow_Giv(i,j)))/2;  %Higher weighting for ramp flows.
                    if MultFac<1
                        d_vertex = (Term1 - FlowBet(i,j))/fac;
                    else
                        d_vertex = (Term1 + FrFlow_Giv(i,j) - fac*Term2)/fac;
                    end
                    
                    if all(ORBoundsUpper(:,j)==0)
                        d_soln = min(c(i,j),max([0,dprev(j),d_vertex,c(i,j)-Term2]));
                    else
                        d_soln = min(DemandUpper(i,j+1),max([DemandLower(i,j+1),dprev(j),d_vertex,c(i,j)-Term2])); % BOUNDS APPLIED
                    end
                    
                     if all(FRBoundsUpper(:,j)==0)
                        fr_soln = Term2 + d_soln - c(i,j);
                     else
                         fr_soln = min(FRBoundsUpper(i,j),max(FRBoundsLower(i,j),Term2 + d_soln - c(i,j))); % BOUNDS APPLIED
                     end
                    
                end
                dj1(i,j)=d_soln;
                BETA1(i,j)=fr_soln/Term2;
                OrInp(i,j)=d_soln-dprev(j);
                dprev(j)=d_soln-Term1/c(i,j)*dj1(i,j);
                
            else %(FREE FLOW)
                
                %%% Old alg
                if ~FrPresent(j)
                     if all(ORBoundsUpper(:,j)==0)
                        dj1(i,j)=c(i,j)-Term2;
                     else
                         dj1(i,j) = min(DemandUpper(i,j+1),max(DemandLower(i,j+1),c(i,j)-Term2)); % BOUNDS APPLIED
                     end
                    OrInp(i,j)=dj1(i,j)-dprev(j);
                    OnrampFlw(i,j)=dj1(i,j);
                    dprev(j)=0;
                    continue;
                end
                Term2= min(Qmax(j),Nha(i,j)*vf(j));
                % Explicit solution
                if ImputeOR(j) && ImputeFR(j) || Override
                    
                    if all(FRBoundsUpper(:,j)==0)
                        s_soln = min(Term2,min(Term2,max([0,Term2 - FlowBet(i,j),Term2-c(i,j)])));
                    else
                        s_soln = min(FRBoundsUpper(i,j),max([FRBoundsLower(i,j),Term2 - FlowBet(i,j),Term2-c(i,j)]));
                    end
                    
                     if all(ORBoundsUpper(:,j)==0)
                        r_soln = c(i,j)-Term2+ s_soln;
                     else
                         r_soln = min(ORBoundsUpper(i,j),max(ORBoundsLower(i,j),c(i,j)-Term2+ s_soln));
                     end
                    
                elseif ~ImputeOR(j)
                    
                    if all(FRBoundsUpper(:,j)==0)
                        s_bnd = max([0,Term2-c(i,j)]);
                    else
                        s_bnd = max([FRBoundsLower(i,j),Term2-c(i,j)]);
                    end
                    MultFac = abs(FlowBet(i,j)/max(100*SimT,OrFlow_Giv(i,j)))/2;  %Higher weighting for ramp flows.
                    if MultFac<1
                        s_vertex = Term2 - FlowBet(i,j);
                    else
                        s_vertex = Term2+OrFlow_Giv(i,j)-c(i,j);
                    end
                    
                    if all(FRBoundsUpper(:,j)==0)
                        s_soln = min(max(s_vertex,s_bnd),Term2);
                    else
                        s_soln = min(max(s_vertex,s_bnd),FRBoundsUpper(i,j));
                    end
                    
                     if all(ORBoundsUpper(:,j)==0)
                        r_soln = c(i,j)-Term2+ s_soln;
                     else
                         r_soln = min(ORBoundsUpper(i,j),max(ORBoundsLower(i,j),c(i,j)-Term2+ s_soln));
                     end
                    
                else
                    
                    if all(FRBoundsUpper(:,j)==0)
                        s_bnd = max([0,Term2-c(i,j)]);
                    else
                        s_bnd = max([FRBoundsLower(i,j),Term2-c(i,j)]);
                    end
                    MultFac = abs(FlowBet(i,j)/max(100*SimT,FrFlow_Giv(i,j)))/2;  %Higher weighting for ramp flows.
                    if MultFac<1
                        s_vertex = Term2 - FlowBet(i,j);
                    else
                        s_vertex = FrFlow_Giv(i,j);
                    end
                    
                    if all(FRBoundsUpper(:,j)==0)
                        s_soln = min(max(s_vertex,s_bnd),Term2);
                    else
                        s_soln = min(max(s_vertex,s_bnd),FRBoundsUpper(i,j));
                    end
                    
                     if all(ORBoundsUpper(:,j)==0)
                        r_soln = c(i,j)-Term2+ s_soln;
                     else
                         r_soln = min(ORBoundsUpper(i,j),max(ORBoundsLower(i,j),c(i,j)-Term2+ s_soln));
                     end
                    
                end
                dj1(i,j)=r_soln;
                BETA1(i,j)=s_soln/Term2;
                OrInp(i,j)=r_soln;
                dprev(j)=0;
            end
        end
    end
    
end
% ===========================LP WITH BOUNDS===============================

%% SIMULATION
% Simulate with the estimated demands and split ratios
% The previous step gave us demands from side streets to onramps. This step
% is to find the demands from onramps into freeway.

BETAF=max(0,min(BETA+BETA1,1));
clear BETA BETA1
DJ=dj1+Demand;
clear dj1 Demand

if DownstreamCongested
    DJBound=c(:,end)-min(Qmax(end),Nha(:,end)*vf(end));
    DJBound(c(:,end)<=1e-4+Qmax(end))=0;
end

% Using the DJ, run the model to get onramp input flows
Nh=zeros(numtime,numcell);

qj=zeros(numtime,numcell);
Inflow=zeros(1,numcell);
Outflow=zeros(1,numcell);

InQ=0; %Queue for the inflow
warning('off')
Nh(1,:)=Density(1,:);
dprev=zeros(1,numcell);
OnrampInput=zeros(numtime,numcell);
Bounddprev=0;
for i=2:length(Flow)
    Nh(i,:)=Nh(i-1,:)+Inflow-Outflow;
    
    Outflow(end)=min(Nh(i,end)*vf(end),Qmax(end));
    
    if DownstreamCongested
        if (DJBound(i)+Outflow(end))>Qmax(end)            
            BoundaryRampFlw(i)= Qmax(end)/(DJBound(i)+Outflow(end))*DJBound(i);
            BoundaryRampInp(i)=DJBound(i)-Bounddprev;

            Bounddprev=DJBound(i)-BoundaryRampFlw(i);
            Outflow(end)= Qmax(end)/(DJBound(i)+Outflow(end))*Outflow(end);
        else
            BoundaryRampFlw(i)= DJBound(i);
            BoundaryRampInp(i)=DJBound(i)-Bounddprev;
            Bounddprev=DJBound(i)-BoundaryRampFlw(i);
        end
    end
    
    
    InQ=InQ+InputFlw(i);
    
    Capacity=w(1)*(rhojam(1)-Nh(i,1));
    Inflow(1)=min(Capacity,InQ);
    InQ=InQ-Inflow(1);
    
    for j=2:NoNodes-1
        Beta=[1-BETAF(i,j-1) 1;BETAF(i,j-1) 0];
        Beta=round(Beta*1e4)/1e4;
        Demands=[min(Qmax(j-1),Nh(i,j-1)*vf(j-1)) DJ(i,j-1)]';
        
        NoOut=2;
        DemandOut=(Beta*Demands);
        Capacities=[min(Qmax(j),w(j)*(rhojam(j)-Nh(i,j))) 10000]';
        Dij=Beta.*(ones(NoOut,1)*Demands');
        if DemandOut(1)~=0
            Dij(1,:)=Dij(1,:)*min(DemandOut(1),Capacities(1))./DemandOut(1);
        end
        if DemandOut(2)~=0
            Dij(2,:)=Dij(2,:)*min(DemandOut(2),Capacities(2))./DemandOut(2);
        end
        Demands=sum(Dij,1);
        AdjFact=min(Dij./((ones(NoOut,1)*Demands).*Beta),[],1);
        AdjFact(isnan(AdjFact))=0;
        flow=Demands'.*AdjFact';
        Outflow(j-1)=flow(1);
        OrFlow(i,j-1)=flow(2);
        
        qj(i,j-1)=qj(i,j-1)-flow(2);
        INFlow=Beta*flow;
        Inflow(j)=INFlow(1);
        FrFlow(i,j-1)=INFlow(2);
        FlowBet(i,j-1)=INFlow(1)-flow(2);
        OnrampInput(i,j-1)=DJ(i,j-1)-dprev(j-1);
        dprev(j-1)=DJ(i,j-1)-flow(2);
    end
end
DensityError=num2str(mean(mean(abs(Density-Nh)))./mean(mean(Nh))*100);
display(['Final Density error = ' DensityError '%']);

% Check Onramps flows calculated in the previous section.
Nh=zeros(numtime,numcell);

qj=zeros(numtime,numcell);
Inflow=zeros(1,numcell);
Outflow=zeros(1,numcell);
ORINP=max(0,OnrampInput);
if DownstreamCongested
    BoundaryRampInp=max(0,BoundaryRampInp);
    DJBoundSav=DJBound*0;
end
BoundaryQ=0;
InQ=0; %Queue for the inflow
warning('off')
Nh(1,:)=Density(1,:);

clear OnrampInput 
for i=2:length(Flow)
    InFl(i,:) = Inflow;
    OutFl(i,:) = Outflow;
    Nh(i,:)=Nh(i-1,:)+Inflow-Outflow;
    
    Outflow(end)=min(Nh(i,end)*vf(end),Qmax(end));
    if DownstreamCongested
        BoundaryQ=BoundaryQ+BoundaryRampInp(i);
        DJBoundSav(i)=BoundaryQ;
        if (BoundaryQ+Outflow(end))>Qmax(end)
            BoundaryRampFlw(i)= Qmax(end)/(BoundaryQ+Outflow(end))*BoundaryQ;
            Outflow(end)= Qmax(end)/(BoundaryQ+Outflow(end))*Outflow(end);
            BoundaryQ=BoundaryQ-BoundaryRampFlw(i);        
        else
            BoundaryRampFlw(i)= BoundaryQ;
            Outflow(end)= Outflow(end);
            BoundaryQ=BoundaryQ-BoundaryRampFlw(i);                    
        end
    end
    
    qj(i,:)=qj(i-1,:)+ORINP(i,:);
    
    InQ=InQ+InputFlw(i);
    Capacity=w(1)*(rhojam(1)-Nh(i,1));
    Inflow(1)=min(Capacity,InQ);
        InQ=InQ-Inflow(1);
    for j=2:NoNodes-1                               % Cant I just do j=2?
        Beta=[1-BETAF(i,j-1) 1;BETAF(i,j-1) 0];
        Beta=round(Beta*1e4)/1e4;
        Demands=[min(Qmax(j-1),Nh(i,j-1)*vf(j-1)) qj(i,j-1)]';
        
        NoOut=2;
        DemandOut=(Beta*Demands);
        Capacities=[min(Qmax(j),w(j)*(rhojam(j)-Nh(i,j))) 10000]';
        Dij=Beta.*(ones(NoOut,1)*Demands');
        if DemandOut(1)~=0
            Dij(1,:)=Dij(1,:)*min(DemandOut(1),Capacities(1))./DemandOut(1);
        end
        if DemandOut(2)~=0
            Dij(2,:)=Dij(2,:)*min(DemandOut(2),Capacities(2))./DemandOut(2);
        end
        Demands=sum(Dij,1);
        AdjFact=min(Dij./((ones(NoOut,1)*Demands).*Beta),[],1);
        AdjFact(isnan(AdjFact))=0;
        flow=Demands'.*AdjFact';
       
        qj(i,j-1)=qj(i,j-1)-flow(2);
        INFlow=Beta*flow;
        Inflow(j)=INFlow(1);
        Outflow(j-1)=flow(1);
        OrFlow(i,j-1)=flow(2);
        FrFlow(i,j-1)=INFlow(2);
        FlowBet(i,j-1)=INFlow(1)-flow(2);
    end
    
end

FrFlow(:,NoNodes-1)=0;
FrPresent(NoNodes-1)=0;
OrPresent(NoNodes-1)=0;

for i=1:size(Density,2)
    Velocity(:,i)=min(70,InFl(:,i)./Nh(:,i)/SimT*CellLengths(i));
end;

DensityError=num2str(mean(mean(abs(Density-Nh)))./mean(mean(Nh))*100);
display(['Final Density error = ' DensityError '%']);

% Get Simulated Flow at expected detector locations
FlowCompare=[InFl(:,1) FlowBet];

FlowError=num2str(mean(mean(abs(FlowCompare-Flow)))./mean(mean(Flow))*100);
display(['Final Flow error = ' FlowError '%']);

%% Add onramp flow and demand profiles to end node
% This is to keep matrix sizes consistent
if DownstreamCongested
    OrFlow(:,numcell) = BoundaryRampFlw;
    ORINP(:,numcell) = BoundaryRampInp;
else
    OrFlow(:,numcell) = zeros(numtime,1);
    ORINP(:,numcell) = zeros(numtime,1);
end