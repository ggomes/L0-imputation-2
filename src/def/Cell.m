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
        
        function [x]=get_total_length_miles(obj)
            x = obj.get_total_length_meters/1604;
        end
        
        function [x]=get_total_length_meters(obj)
            ml_links = [obj.segments.ml_link];
            x = sum([ml_links.length_meters]);
        end
        
        function [x]=has_or(obj)
            x = false;
            for i=1:length(obj.segments)
                x = x | ~isempty(obj.segments(i).or_links);
            end
            return
        end
        
        function [x]=has_fr(obj)
            x = false;
            for i=1:length(obj.segments)
                x = x | ~isempty(obj.segments(i).fr_links);
            end
            return
        end
        
        function [x]=has_major_or(obj)
            x = false;
            for i=1:length(obj.segments)
                x = x | obj.segments(i).has_major_or;
            end
        end
        
        function [x]=has_major_fr(obj)
            x = false;
            for i=1:length(obj.segments)
                x = x | obj.segments(i).has_major_fr;
            end
        end
        
          
                    
    end
    
end

