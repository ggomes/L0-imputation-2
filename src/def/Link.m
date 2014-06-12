classdef Link < hgsetget
    
    properties
        
        id
        mySegment           % reference to containing segment
        lanes
        length_meters       % length in meters
        
    end
    
    methods
        
        function [obj] = set_att(obj,segment,xlink)
            obj.mySegment = segment;
            obj.id = xlink.ATTRIBUTE.id;
            obj.lanes = xlink.ATTRIBUTE.lanes;
            obj.length_meters = xlink.ATTRIBUTE.length;
        end

            
    end
    
end

