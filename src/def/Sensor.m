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
                                          'dty',[],...          % veh/mile
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
                
                time = 86400*(data.time-day);
                flw = sum(data.flw,2);
                spd = sum(data.flw.*data.spd,2)./flw;
                dty = flw./spd;
                
                % smooth
                span = 11;
                degree = 2;
                flw = smooth(flw,span,'sgolay',degree);
                spd = smooth(spd,span,'sgolay',degree);
                dty = smooth(dty,span,'sgolay',degree);

                obj.five_min = [obj.five_min ...
                                struct( 'day',day,...
                                        'time',time,...
                                        'flw',flw,...
                                        'dty',dty,...
                                        'spd',spd) ];
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

        function [x]=get_5min_data(obj,day)
            x = [];
            dind = [obj.five_min.day]==day;
            if(~any(dind))
                return
            end
            D = obj.five_min(dind);
            x.flw = D.flw; 
            x.dty = D.dty; 
            x.spd = D.spd; 
        end
        
        function [x]=get_5min_data_time(obj,day,time)
            x = [];
            dind = [obj.five_min.day]==day;
            if(~any(dind))
                return
            end
            D = obj.five_min(dind);
            ind = nan(1,length(time));
            dt = min(diff(D.time));
            for i=1:length(time)
                ind(i) = find(time(i)>=D.time-dt/4,1,'last');
            end
            x.flw = D.flw(ind,:); 
            x.dty = D.dty(ind,:); 
            x.spd = D.spd(ind,:); 
        end
        
        function [x]=get_5min_data_time_type(obj,day,time,data_type)
            x = [];
            dind = [obj.five_min.day]==day;
            if(~any(dind))
                return
            end
            D = obj.five_min(dind);
            ind = nan(1,length(time));
            dt = min(diff(D.time));
            for i=1:length(time)
                ind(i) = find(time(i)>=D.time-dt/4,1,'last');
            end
            x = D.(data_type)(ind,:); 
        end
       
    end
    
end

