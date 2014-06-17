classdef CellArray < hgsetget
    
    properties
        cells
    end
    
    methods
        
        function [obj]=add_cell(obj,cell)
            obj.cells = [obj.cells cell];
        end
        
        function [obj]=add_segment_to_last_cell(obj,seg)
            c = length(obj.cells);
            if(c>0)
                obj.cells(c).add_segment(seg);
            end
        end
        
        function [x]=get_num_seg_per_cell(obj)
            x = zeros(1,length(obj.cells));
            for i=1:length(obj.cells)
               x(i) = length(obj.cells(i).segments);
            end
        end
        
    end
    
end

