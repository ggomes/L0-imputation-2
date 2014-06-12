classdef Freeway
    
    properties
        
        seg;
        
    end
    
    methods
        
        function obj = Freeway(sc_ptr)
            
            str_fr_type = 'Off-Ramp';
            str_or_type = 'On-Ramp';
            
            link_ids = sc_ptr.get_link_ids;
            link_types = sc_ptr.get_link_types;
            node_ids = sc_ptr.get_node_ids;
            ordered_ind = extract_linear_fwy_indices(sc_ptr);
            num_segments = length(ordered_ind);
            
            nodes = repmat(Node,1,num_segments+1);
            obj.seg = repmat(Segment,1,num_segments);
            
            for i=1:num_segments
                
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
                ml_link.set_att(obj.seg(i),xmllink);

                % create hov link
                hv_link = repmat(Link,1,0);
                
                % create onramp links
                or_links = repmat(Link,1,length(xorlinks));
                for j=1:length(xorlinks)
                    or_links(j).set_att(obj.seg(i),xorlinks(j));
                end
                   
                % create offramp links
                fr_links = repmat(Link,1,length(xfrlinks));
                for j=1:length(xfrlinks)
                    fr_links(j).set_att(obj.seg(i),xfrlinks(j));
                end
                
                % set node attributes
                set(nodes(i),'dnSegment',obj.seg(i));
                set(nodes(i),'id',xupnode.ATTRIBUTE.id);
                set(nodes(i+1),'upSegment',obj.seg(i));
                set(nodes(i+1),'id',xdnnode.ATTRIBUTE.id);
                
                % set segment attributes
                set(obj.seg(i),'up_node',nodes(i));
                set(obj.seg(i),'dn_node',nodes(i+1));
                set(obj.seg(i),'ml_link',ml_link);
                set(obj.seg(i),'hv_link',hv_link);
                set(obj.seg(i),'or_links',or_links);
                set(obj.seg(i),'fr_links',fr_links);

            end
            

            
            
        end
        
        
    end
    
end

