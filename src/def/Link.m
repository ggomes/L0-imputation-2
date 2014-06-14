classdef Link < hgsetget
    
    properties
        
        id
        mySegment           % reference to containing segment
        type                % ml,hov,or,fr
        lanes
        length_meters       % length in meters
        
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
            
    end
    
end

