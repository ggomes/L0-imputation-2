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
    
    %  construction / setters .............................................
    methods ( Access = public )
        
        function [obj] = Sensor()
            obj.disabled = false;
            obj.health = struct( ...
                        'days',[], ...
                        'percent_observed',[] );
                    
            obj.five_min = repmat(struct( 'day',nan, ...
                                          'time',[],...
                                          'flw',[],...          % veh/hr
                                          'occ',[],...
                                          'spd',[]) , ...       % mile/hr
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
                                struct( 'day',day,...
                                        'time',86400*(data.time-day),...
                                        'flw',data.flw,...
                                        'occ',data.occ,...
                                        'spd',data.spd) ];
            end
        end
                
        function [obj]=disable(obj)
            obj.disabled = true;
        end
        
        function [obj]=enable(obj)
            obj.disabled = false;
        end
        
    end
    
    % getters .............................................................
    methods ( Access = public )
        
        function [x]=is_good_on_day(obj,day,threshold)
            ind = day==obj.health.days;
            if(any(ind))
                x = ~obj.disabled & obj.health.percent_observed(ind)>=100*threshold;
            else
               x=nan; 
            end
        end

        function [x]=get_5min_totals(obj,day,time,data_type)
            x = obj.get_5min_perlane(day,time,data_type);
            switch(data_type)
                case 'flw'
                    x = sum(x,2);
                case 'occ'
                    x = mean(x,2);
                case 'spd'
                    x = mean(x,2);
                otherwise
                    error('wrong data type')
            end
        end
        
        function [x]=get_5min_perlane(obj,day,time,data_type)
            x = [];
            dind = [obj.five_min.day]==day;
            if(~any(dind))
                return
            end
            D = obj.five_min(dind);
            ind = nan(1,length(time));
            for i=1:length(time)
                ind(i) = find(time(i)>D.time,1,'last');
            end
            x = D.(data_type)(ind,:);            
        end
         
    end
    
end

