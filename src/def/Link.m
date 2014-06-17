classdef Link < hgsetget
    
    properties ( Access = public )
        
        id
        mySegment           % reference to containing segment
        type                % ml,hv,or,fr
        lanes
        length_meters       % length in meters
        
        fd                  % reference to FundamentalDiagram
        sensor              % reference to sensor
        
    end
    
    methods(Static)
      function x = type_to_int(xtype)
         switch(xtype)
             case 'ml'
                 x = 0;
             case 'hv'
                 x = 1;
             case 'or'
                 x = 2;
             case 'fr' 
                 x = 3;
             otherwise
                 x = nan;
         end
      end
      function x = int_to_type(xint)
          switch(xint)
              case 0
                  x = 'ml';
              case 1
                  x = 'hv';
              case 2
                  x = 'or';
              case 3 
                  x = 'fr';
              otherwise
                  x = '';
          end        
      end
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

