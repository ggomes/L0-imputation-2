classdef Freeway < hgsetget
    
    properties ( Access = public )
        
        units           % [us|si]
        seg
        
    end
    
    properties ( Access = private )
        all_links
        all_sensors
    end
    
    % construction ........................................................
    methods ( Access = public )
        
        function obj = Freeway(sc_ptr)
            
            str_fr_type = 'Off-Ramp';
            str_or_type = 'On-Ramp';
            
            sc_ptr.change_units_to('si');
            
            obj.units = 'si';
            
            link_ids = sc_ptr.get_link_ids;
            link_types = sc_ptr.get_link_types;
            node_ids = sc_ptr.get_node_ids;
            ordered_ind = extract_linear_fwy_indices(sc_ptr);
            ordered_fwy_ind = ordered_ind(strcmp(link_types(ordered_ind),'Freeway'));
            
            num_segments = length(ordered_fwy_ind);
            
            nodes(num_segments+1) = Node;
            
            % create segments
            for i=1:num_segments
                
                S = Segment(obj,i);
                
                % get mainline link
                xmllink = sc_ptr.scenario.NetworkSet.network.LinkList.link(ordered_fwy_ind(i));
                
                if(~strcmpi(xmllink.link_type.ATTRIBUTE.name,'Freeway'))
                    continue;
                end
                
                % get upstream node
                node_ind = xmllink.begin.ATTRIBUTE.node_id==node_ids;
                if(~any(node_ind))
                    error('bad node id')
                end
                xupnode = sc_ptr.scenario.NetworkSet.network.NodeList.node(node_ind);
                
                % get downstream node
                node_ind = xmllink.end.ATTRIBUTE.node_id==node_ids;
                if(~any(node_ind))
                    error('bad node id')
                end
                xdnnode = sc_ptr.scenario.NetworkSet.network.NodeList.node(node_ind);
                
                % get onramp links
                xorlinks = [];
                if(~isempty(xupnode.inputs))
                    n = length(xupnode.inputs.input);
                    in_ids = nan(1,n);
                    for j=1:n
                        in_ids(j) = xupnode.inputs.input(j).ATTRIBUTE.link_id;
                    end
                    [~,in_ind] = ismember(in_ids,link_ids);
                    or_ind = in_ind(strcmp(link_types(in_ind),str_or_type));
                    xorlinks = sc_ptr.scenario.NetworkSet.network.LinkList.link(or_ind);
                end
                
                % get offramp links
                xfrlinks = [];
                if(~isempty(xdnnode.outputs))
                    n = length(xdnnode.outputs.output);
                    out_ids = nan(1,n);
                    for j=1:n
                        out_ids(j) = xdnnode.outputs.output(j).ATTRIBUTE.link_id;
                    end
                    [~,out_ind] = ismember(out_ids,link_ids);
                    fr_ind = out_ind(strcmp(link_types(out_ind),str_fr_type));
                    xfrlinks = sc_ptr.scenario.NetworkSet.network.LinkList.link(fr_ind);
                end
                
                % create mainline link
                ml_link = Link;
                ml_link.set_att(S,xmllink,'ml','si');
                
                % create hov link
                hv_link = [];
                
                % create onramp links
                if(~isempty(xorlinks))
                    clear or_links
                    or_links(length(xorlinks)) = Link;
                    for j=1:length(xorlinks)
                        or_links(j).set_att(S,xorlinks(j),'or','si');
                    end
                else
                    or_links = [];
                end
                
                % create offramp links
                if(~isempty(xfrlinks))
                    clear fr_links
                    fr_links(length(xfrlinks)) = Link;
                    for j=1:length(xfrlinks)
                        fr_links(j).set_att(S,xfrlinks(j),'fr','si');
                    end
                else
                    fr_links = [];
                end
                
                
                % set node attributes
                set(nodes(i),'dnSegment',S);
                set(nodes(i),'id',xupnode.ATTRIBUTE.id);
                set(nodes(i+1),'upSegment',S);
                set(nodes(i+1),'id',xdnnode.ATTRIBUTE.id);
                
                % set segment attributes
                set(S,'up_node',nodes(i));
                set(S,'dn_node',nodes(i+1));
                set(S,'ml_link',ml_link);
                set(S,'hv_link',hv_link);
                set(S,'or_links',or_links);
                set(S,'fr_links',fr_links);
                
                obj.seg = [obj.seg S];
                
            end
            
            obj.all_links = [[obj.seg.ml_link] [obj.seg.hv_link] [obj.seg.or_links] [obj.seg.fr_links]];
            
            % assign sensors
            if(isfieldRecursive(sc_ptr.scenario,'SensorSet','sensor'))
                sensor_link = sc_ptr.get_sensor_link_map;
                for i=1:size(sensor_link,1)
                    sensor = Sensor;
                    sensor.set_att(sc_ptr.scenario.SensorSet.sensor(i));
                    myLink = obj.get_link_for_id(sensor_link(i,2));
                    if(~isempty(myLink))
                        set(sensor,'myLink',myLink);
                        myLink.add_sensor(sensor);
                    end
                end
            end
            obj.all_sensors = [obj.all_links.sensor];
            
            % assign fundamental diagrams
            if(isfieldRecursive(sc_ptr.scenario,'FundamentalDiagramSet','fundamentalDiagramProfile'))
                fdp = sc_ptr.scenario.FundamentalDiagramSet.fundamentalDiagramProfile;
                for i=1:length(fdp)
                    link = obj.get_link_for_id(fdp(i).ATTRIBUTE.link_id);
                    if(~isempty(link))
                        link.assign_fundamental_diagram_from_xml(fdp(i).fundamentalDiagram,'si');
                    end
                end
            end
            
            
        end
        
        % set is_major to true for links with id link_id if they are of 
        % type link_type
        function obj = set_major(obj,link_id,link_type)
           [have,ind]=ismember(link_id,[obj.all_links.id]);
           ind = ind( have & strcmp({obj.all_links(ind).type},link_type) );
           for i=ind
               obj.all_links(i).is_major = true;
           end
        end
    end
    
    % load data ...........................................................
    methods ( Access = public )
        
        function obj = load_loop_health_for_days(obj,folder,days)
            
            filename = fullfile(folder,'loop_health');
            if(~exist([filename '.mat'],'file'))
                error('loop health mat file does not exist')
            end
            
            load(filename)

            keep_days = intersect(days,loop_health.days);
            keep_vds = intersect([obj.all_sensors.vds],loop_health.vds);
            
            [~,days_ind]=ismember(keep_days,loop_health.days);
            for i=1:length(keep_vds)
                vds = keep_vds(i);
                sensor = obj.get_sensor_with_vds(vds);
                percent_observed = loop_health.percent_observed(loop_health.vds==vds,days_ind);
                sensor.append_health_info( keep_days' , percent_observed );
            end
            
        end
        
        function obj = load_5min_data_for_days(obj,folder,days)
            for i=1:length(days)
                filename = fullfile(folder,['pems5min_' datestr(days(i),'yyyy_mm_dd')]);
                if(~exist([filename '.mat'],'file'))
                    warning('loop 5-min mat file does not exist')
                    continue
                end
                load(filename)
                for j=1:length(obj.all_sensors)
                    ind = obj.all_sensors(j).vds==pems.vds;
                    if(~any(ind))
                        continue
                    end
                    obj.all_sensors(j).append_5min_data_day(days(i),pems.data(ind));
                end       
            end
        end
        
    end
    
    % getters .............................................................
    methods ( Access = public )
        
        % links ..................
        function L = get_link_for_id(obj,link_id)
            ind = [obj.all_links.id]==link_id;
            if(any(ind))
                L = obj.all_links(ind);
            else
                L = [];
            end
        end
        
        % sensors ................
        function S = get_sensor_with_vds(obj,vds)
        % sensor object corresponding to given vds
            ind=[obj.all_sensors.vds]==vds;
            if(any(ind))
                S=obj.all_sensors(ind);
            else
                S=[];
            end
        end
        
        function S = get_sensor_with_linkid(obj,link_id)
        % sensor object for a given link id
            ind = [obj.all_links.id]==link_id;
            if(any(ind))
                S=obj.all_links(ind).sensor;
            else
                S=[];
            end
        end
        
        function [x] = get_sensor_health_for_day(obj,link_type,day,threshold)  
        % [1xnumsegment] array of sensor health values for a given day and link type      
            x = false(1,length(obj.seg));
            for i=1:length(obj.seg)
                x(i)=obj.seg(i).has_good_sensor_on_day(link_type,day,threshold);
            end
        end
        
        function [x] = get_all_sensor_health_for_day(obj,day,threshold)
        % returns table of sensor health information. Each row is a sensor,
        % columns are [segment index,link id,link type,sensor id,sensor vds,sensor health]
            x = nan(length(obj.all_sensors),6);
            for i=1:length(obj.all_sensors)
               S = obj.all_sensors(i);
               if(~isempty(S.myLink))
                   x(i,1) = S.myLink.mySegment.index;
                   x(i,2) = S.myLink.id;
                   x(i,3) = Link.type_to_int(S.myLink.type);
               end
               x(i,4) = S.id;
               x(i,5) = S.vds;
               x(i,6) = S.is_good_on_day(day,threshold);
            end
        end
    end
        
    % setters .............................................................
    methods ( Access = public )
        
        % sensors ................
        function [obj] = disable_sensor_vds(obj,vdslist)
            for i=1:length(vdslist)
                S = obj.get_sensor_with_vds(vdslist(i));
                if(~isempty(S))
                    S.disable();
                end
            end
        end
        
        function [obj] = enable_sensor_vds(obj,vdslist)
            for i=1:length(vdslist)
                S = obj.get_sensor_with_vds(vdslist(i));
                if(~isempty(S))
                    S.enable();
                end
            end
        end


        
    end
    
end

