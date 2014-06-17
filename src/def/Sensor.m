classdef Sensor < hgsetget
    
    properties
        id
        myLink
        lanes
        vds  
        
        disabled
        
        five_min
        health
        
    end
    
    methods
        
        function [obj] = Sensor()
            obj.disabled = false;
            obj.health = struct( ...
                        'days',[], ...
                        'percent_observed',[] );
                    
            obj.five_min = repmat(struct( 'day',nan, ...
                                          'time',[],...
                                          'flw_vph',[],...
                                          'occ',[],...
                                          'spd_mph',[]) , ...
                                  1,0);       
        end
        
        function [obj] = set_att(obj,xsensor)
            obj.id = xsensor.ATTRIBUTE.id;
            obj.lanes = xsensor.ATTRIBUTE.lane_number;
            obj.vds = xsensor.ATTRIBUTE.sensor_id_original;
        end
        
        function [obj]=append_health_info(obj,days,percent_observed)
            add_days = ~ismember(days,obj.health.days);
            obj.health.days = [obj.health.days days(add_days)];
            obj.health.percent_observed = [obj.health.percent_observed percent_observed(add_days)];
        end
        
        function [obj]=append_5min_data_day(obj,day,data)
            if(isempty(obj.five_min) || ~ismember(day,[obj.five_min.day]) )
                obj.five_min = [obj.five_min ...
                                struct( 'day',day,'time',data.time,...
                                        'flw_vph',data.flw,...
                                         'occ',data.occ,...
                                         'spd_mph',data.spd) ];
            end
        end
        
        function [x]=is_good_on_day(obj,day,threshold)
            ind = day==obj.health.days;
            if(any(ind))
                x = ~obj.disabled & obj.health.percent_observed(ind)>=100*threshold;
            else
               x=nan; 
            end
        end
        
        function [obj]=disable(obj)
            obj.disabled = true;
        end
        function [obj]=enable(obj)
            obj.disabled = false;
        end
        
    end
    
end

