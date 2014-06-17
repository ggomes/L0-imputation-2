classdef FundamentalDiagram < hgsetget
    
    properties
        units               % si|us
        vf                  % [meter/sec] freeflow speed
        w                   % [meter/sec] congestion speed
        f_max               % [veh/sec] capacity
        dens_jam            % [veh/sec] jam density
        dens_crit           % [veh/meter] ciritical density
    end
    
    methods
        
        function [obj]=FundamentalDiagram(A,units)
            switch(units)
                case 'si'
                    convt=1;
                    convd=1;
                case 'us'
                    convt=3600;         % seconds in one hour
                    convd=1609.34;      % meters in one mile
                otherwise
                    error('unknown units')
            end   
            obj.units = units;         
            obj.vf = A.free_flow_speed*convd/convt;
            obj.w = A.congestion_speed*convd/convt;
            obj.f_max = A.capacity/convt;
            obj.dens_crit =A.capacity/A.free_flow_speed/convd;
            obj.dens_jam = A.capacity*(1/A.free_flow_speed + 1/A.congestion_speed)/convd;
        end

    end
    
end

