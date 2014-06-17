classdef Cell < hgsetget
    
    properties
        
        segments
        
    end
    
    methods
        
        function [obj]=add_segment(obj,seg)
            obj.segments = [obj.segments seg];
        end
                    
    end
    
end

