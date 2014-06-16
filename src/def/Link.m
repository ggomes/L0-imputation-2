classdef Link < hgsetget
    
    properties ( Access = public )
        
        id
        mySegment           % reference to containing segment
        type                % ml,hov,or,fr
        lanes
        length_meters       % length in meters
                
        vf                  % [mph] freeflow speed
        w                   % [mph] congestion speed
        f_max               % [vph] capacity
        dens_jam            % [veh/mile] jam density
        dens_crit           % [veh/mile] ciritical density
                
        sensor              % reference to sensor
        
    end
    
    methods ( Access = public )
        
        function [obj] = set_att(obj,segment,xlink,type)
            obj.mySegment = segment;
            obj.type = type;
            obj.id = xlink.ATTRIBUTE.id;
            obj.lanes = xlink.ATTRIBUTE.lanes;
            obj.length_meters = xlink.ATTRIBUTE.length;

        end
        
        function [obj]=add_sensor(obj,S)
            obj.sensor = [obj.sensor S];
        end
        
        function [obj]=assign_fundamental_diagram_from_xml(obj,fd)
                A = fd.ATTRIBUTE;
                obj.vf = A.free_flow_speed; %*convt/convd;
                obj.w = A.congestion_speed; %*convt/convd;
                obj.f_max = A.capacity; %*convt;
                obj.dens_crit =A.capacity/A.free_flow_speed; %*convd;
                obj.dens_jam = A.capacity*(1/A.free_flow_speed + 1/A.congestion_speed); %*convd;
        end
            
    end
    
end

