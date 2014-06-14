classdef Freeway
    
    properties ( Access = public )
        
        seg;
        
    end
    
    properties ( Access = private )
        all_links
    end
    
    methods ( Access = public )
        
        function obj = Freeway(sc_ptr)
            
            str_fr_type = 'Off-Ramp';
            str_or_type = 'On-Ramp';
            
            link_ids = sc_ptr.get_link_ids;
            link_types = sc_ptr.get_link_types;
            node_ids = sc_ptr.get_node_ids;
            ordered_ind = extract_linear_fwy_indices(sc_ptr);
            num_segments = length(ordered_ind);
            
            nodes = repmat(Node,1,num_segments+1);
            
            % create segments
            for i=1:num_segments
                
                S = Segment;
                
                % get mainline link
                xmllink = sc_ptr.scenario.NetworkSet.network.LinkList.link(ordered_ind(i));
                
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
                ml_link.set_att(S,xmllink,'ml');
                
                % create hov link
                hv_link = repmat(Link,1,0);
                
                % create onramp links
                or_links = repmat(Link,1,length(xorlinks));
                for j=1:length(xorlinks)
                    or_links(j).set_att(S,xorlinks(j),'or');
                end
                
                % create offramp links
                fr_links = repmat(Link,1,length(xfrlinks));
                for j=1:length(xfrlinks)
                    fr_links(j).set_att(S,xfrlinks(j),'fr');
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
                    ind=sensor.id==sensor_link(:,1);
                    if(any(ind))
                        myLink = obj.all_links(ind);
                        set(sensor,'myLink',myLink);
                        myLink.add_sensor(sensor);
                    end
                end
            end
            
        end
        
        function L = get_link_for_id(obj,link_id)
            ind = [obj.all_links.id]==link_id;
            if(any(ind))
                L = obj.all_links(ind);
            else
                L = repmat(Link,1,0);
            end
        end
       
        function S = get_sensor_with_vds(obj,vds)
            all_sensors = [obj.all_links.sensor];
            ind=[all_sensors.vds]==vds;
            if(any(ind))
                S=all_sensors(ind);
            else
                S=repmat(Sensor,1,0);
            end
        end
        
        function S = get_sensor_with_linkid(obj,link_id)
            ind = [obj.all_links.id]==link_id;
            if(any(ind))
                S=obj.all_links(ind).sensor;
            else
                S=repmat(Sensor,1,0);
            end
        end
        
    end
end

