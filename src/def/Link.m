classdef Link < hgsetget
    
    properties ( Access = public )
        
        id
        mySegment           % reference to containing segment
        type                % ml,hov,or,fr
        lanes
        length_meters       % length in meters
        
        fd                  % reference to FundamentalDiagram
        sensor              % reference to sensor
        
    end
    
    methods ( Access = public )
        
        function [obj] = set_att(obj,segment,xlink,type,units)
            obj.mySegment = segment;
            obj.type = type;
            obj.id = xlink.ATTRIBUTE.id;
            obj.lanes = xlink.ATTRIBUTE.lanes;
            switch(units)
                case 'si'
                    obj.length_meters = xlink.ATTRIBUTE.length;
                case 'us'
                    obj.length_meters = xlink.ATTRIBUTE.length*1609.34;
                otherwise
                    error('unknown units')
            end
        end
        
        function [obj]=add_sensor(obj,S)
            obj.sensor = [obj.sensor S];
        end
        
        function [obj]=assign_fundamental_diagram_from_xml(obj,fd,units)
            obj.fd = FundamentalDiagram(fd.ATTRIBUTE,units);
        end
        
        function [x]=has_good_sensor_on_day(obj,day,threshold)
            if(~isempty(obj.sensor))
                x = obj.sensor.is_good_on_day(day,threshold);
            else
                x = false;
            end
            if(isnan(x))
                x = false;
            end
        end
        
    end
    
end

