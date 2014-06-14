classdef Sensor < hgsetget
    
    properties
        id
        myLink
        lanes
        vds      
    end
    
    methods
        
        function [obj] = set_att(obj,xsensor)
            obj.id = xsensor.ATTRIBUTE.id;
            obj.lanes = xsensor.ATTRIBUTE.lane_number;
            obj.vds = xsensor.ATTRIBUTE.sensor_id_original;
        end
        
    end
    
end

