classdef Cell < hgsetget
    
    properties
        
        segments
        
        % imputation specific
        ml_fd         % fd to be applied to the cell
        ml_lanes      % mainline lanes
        
    end
    
    methods
        
        function [obj]=add_segment(obj,seg)
            obj.segments = [obj.segments seg];
        end
                    
    end
    
end

