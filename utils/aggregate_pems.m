function [agg_data] = aggregate_pems(data)

agg_data = struct('vds',data.vds,'lanes',size(data.flw,2),'flw',[],'dty',[],'spd',[],'time',data.time);
agg_data.flw = sum(data.flw,2);
agg_data.spd = sum(data.flw.*data.spd,2)./agg_data.flw;
agg_data.dty = agg_data.flw./agg_data.spd;
