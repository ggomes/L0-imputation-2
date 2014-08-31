classdef Segment < hgsetget
    
    properties
        
        myFreeway
        index         % order in the freeway, starting from 1
        
        up_node;      % reference to uptream node
        dn_node;      % reference to downstream node
        ml_link;      % reference to mainline link
        hv_link;      % reference to hov link
        or_links;     % reference to onramp links   
        fr_links;     % reference to offramp links   
        
    end
    
    methods
        
        function [obj]=Segment(fwy,ind)
            obj.index = ind;
            obj.myFreeway = fwy;
        end
        
        function [x]=has_good_sensor_on_day(obj,link_type,day,threshold)
            switch(link_type)
                case 'ml'
                    link = obj.ml_link;
                case 'hv'
                    link = obj.hv_link;
                case 'or'
                    link = obj.or_links;
                case 'fr'
                    link = obj.fr_links;
                otherwise
                    error('bad inptu parameter')
            end
            if(~isempty(link))
                x = link.has_good_sensor_on_day(day,threshold);
            else
                x = false;
            end
        end
        
        function [x]=has_major_or(obj)
           x = any([obj.or_links.is_major]);
        end
        
        function [x]=has_major_fr(obj)
           x = any([obj.fr_links.is_major]);
        end
        
        function [x]=get_or_dwnstr(obj)
            x = [];
            if(obj.index==length(obj.myFreeway.seg))
                return;
            end
            x = obj.myFreeway.seg(obj.index+1).or_links;
        end

    end
    
end

