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
                    link = obj.or_link;
                case 'fr'
                    link = obj.fr_link;
                otherwise
                    error('bad inptu parameter')
            end
            if(~isempty(link))
                x = link.has_good_sensor_on_day(day,threshold);
            else
                x = false;
            end
        end
        
        
    end
    
end

