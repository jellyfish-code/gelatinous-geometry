%% Remesh
function [jelly, lim_reached] = remesh_SLM_butterfly(jelly, muscle_length)
lim_reached = 0;
jelly_temp = jelly;
mus_length_cur = sum(jelly.Edges.d_current(jelly.Edges.muscle == 1));
lim = numedges(jelly);
count = 0;
    %% Delete nodes when needed
    while min(jelly.Edges.d_current) < 0.2
        count = count + 1;
        if count >= lim
            lim_reached = 1;
            return
        end
        [~, idx] = sortrows(jelly.Edges, {'d_current'});
        i = idx(1);    
        a = jelly.Edges.EndNodes(i,:);
        node1 = a(1);
        node2 = a(2);
        nodes_for_deletion = [node1 node2];
        

        %First, define the new node as halfway between nodes1 and 2 (same
        %as above)
        new_Name = (jelly.Nodes.Node_Names(node1,:) + jelly.Nodes.Node_Names(node2,:))./2;
        new_x = (jelly.Nodes.x_coord(node1) + jelly.Nodes.x_coord(node2))./2;
        new_y = (jelly.Nodes.y_coord(node1) + jelly.Nodes.y_coord(node2))./2;
        %Is the new node on the edge? It is replacing one edge node or two
        %edge nodes
        if jelly.Nodes.edges(node1) > 0 && jelly.Nodes.edges(node2) > 0
            if jelly.Nodes.edges(node1) == 1 || jelly.Nodes.edges(node2) == 1
                new_edge = max(jelly.Nodes.edges) + 1;
            else
                new_edge = (jelly.Nodes.edges(node1) + jelly.Nodes.edges(node2))/2;
            end
        elseif jelly.Nodes.edges(node1) > 0
            new_edge = jelly.Nodes.edges(node1);
        elseif jelly.Nodes.edges(node2) > 0
            new_edge = jelly.Nodes.edges(node2);
        else
            new_edge = 0;
        end
        
        if (jelly.Nodes.outmus(node1) > 0 || jelly.Nodes.outmus(node2) > 0) && (jelly.Nodes.inmus(node1) == 0 && jelly.Nodes.inmus(node2) == 0)
            if jelly.Nodes.outmus(node1) > 0 && jelly.Nodes.outmus(node2) > 0
                new_outmus = (jelly.Nodes.outmus(node1) + jelly.Nodes.outmus(node2))/2;
            elseif jelly.Nodes.outmus(node1) > 0
                new_outmus = jelly.Nodes.outmus(node1);
            elseif jelly.Nodes.outmus(node2) > 0
                new_outmus = jelly.Nodes.outmus(node2);
            else
                new_outmus = 0;
            end
            new_inmus = 0;
            if mod(jelly.Nodes.mus_num(node1),1) == 0 && jelly.Nodes.mus_num(node1) ~= 0
                new_mus_num = jelly.Nodes.mus_num(node1);
            elseif mod(jelly.Nodes.mus_num(node2),1)== 0 && jelly.Nodes.mus_num(node2)~= 0
                new_mus_num = jelly.Nodes.mus_num(node2);
            end
        elseif (jelly.Nodes.inmus(node1) > 0 || jelly.Nodes.inmus(node2) > 0) && (jelly.Nodes.outmus(node1) == 0 && jelly.Nodes.outmus(node2) == 0)
            if jelly.Nodes.inmus(node1) > 0 && jelly.Nodes.inmus(node2) > 0
                new_inmus = (jelly.Nodes.inmus(node1) + jelly.Nodes.inmus(node2))/2;
            elseif jelly.Nodes.inmus(node1) > 0
                new_inmus = jelly.Nodes.inmus(node1);
            elseif jelly.Nodes.inmus(node2) > 0
                new_inmus = jelly.Nodes.inmus(node2);
            else
                new_inmus = 0;
            end
            new_outmus = 0;
            if  mod(jelly.Nodes.mus_num(node1),1) == 0 && jelly.Nodes.mus_num(node1) ~= 0
                new_mus_num = jelly.Nodes.mus_num(node1);
            elseif mod(jelly.Nodes.mus_num(node2),1)== 0 && jelly.Nodes.mus_num(node2)~= 0
                new_mus_num = jelly.Nodes.mus_num(node2);
            end
        else
            new_inmus = 0;
            new_outmus = 0;
            new_mus_num = 0;
        end
        
        props = table(new_Name, new_x, new_y, [0,0], [0,0], [0,0], [0,0], [0,0], new_outmus, new_inmus, new_mus_num, new_edge, ...
            'VariableNames', {'Node_Names', 'x_coord', 'y_coord', 'velocity', 'F_net', 'F_elastic', 'F_pressure', 'F_muscle', 'outmus', 'inmus', 'mus_num', 'edges'});
        jelly_temp = addnode(jelly_temp, props);
        newid = numnodes(jelly_temp); %Always added to the back, so the new node is the last one.

        %%Connect it up to connections of parents
        for node_search = 1:numnodes(jelly)
            found = findedge(jelly, node1, node_search);
            found2 = findedge(jelly, node2, node_search);
            if found ~= 0 || found2 ~= 0                
                dx = jelly.Nodes.x_coord(node_search) - jelly_temp.Nodes.x_coord(newid);
                dy = jelly.Nodes.y_coord(node_search) - jelly_temp.Nodes.y_coord(newid);
                d = (dx^2 + dy^2)^(1/2);
                if found ~= 0 && found2 ~= 0
                    s0 = (jelly.Edges.strain0(found) + jelly.Edges.strain0(found2))/2;
                    s1 = (jelly.Edges.strain1(found) + jelly.Edges.strain1(found2))/2;
%                     thick = (jelly.Edges.thickness(found) + jelly.Edges.thickness(found2))/2;
                elseif found ~= 0
                    s0 = jelly.Edges.strain0(found);
                    s1 = jelly.Edges.strain1(found);
%                     thick = jelly.Edges.thickness(found);
                else
                    s0 = jelly.Edges.strain0(found2);
                    s1 = jelly.Edges.strain1(found2);
%                     thick = jelly.Edges.thickness(found2);
                end
                    
                if found ~= 0 && jelly.Edges.muscle(found) ~= 0
                    m = jelly.Edges.muscle(found);
                elseif found2 ~= 0 && jelly.Edges.muscle(found2) ~= 0
                    m = jelly.Edges.muscle(found2);
                else
                    m = 0;
                end

                %%What should the strain on these new edges be?
                props = table([node_search, newid], 1, d, d/(s0+1), d/(s1+1), s0, s1, m, ...
                    'VariableNames', {'EndNodes', 'Weight', 'd_current', 'd_rel0', 'd_rel1', 'strain0', 'strain1', 'muscle'});
                %% new stuff
%                   props = table([node_search, newid], 1, d, d/(s0+1), s0, m, thick, ...
%                       'VariableNames', {'EndNodes', 'Weight', 'd_current', 'd_rel0', 'strain0', 'muscle', 'thickness'});
                %% end new stuff
                jelly_temp = addedge(jelly_temp, props); 
            end
        end
        %%Now delete the parent nodes
        jelly_temp = rmnode(jelly_temp, nodes_for_deletion);
        %This also deletes the edges that parents were connected to
        jelly = jelly_temp;
    
        if new_edge ~= 0
            %%Reorder the edges
            [trash, idx] = sortrows(jelly.Nodes, {'edges'});
            edge_idx = idx(trash.edges~=0); 
            a =  1;
            for i = 1:length(edge_idx)
                jelly.Nodes.edges(edge_idx(i)) = a;
                a = a+1;
            end
            jelly_temp = jelly;
        end
    end
    
    %% Add nodes where needed
    %If the distance is above a certain threshold, add a node in between
    %Node should connect to prev. nodes
    %If node is a muscle, it is still a muscle. If it's an edge, it is still an
    %edge.
    jelly_temp = jelly;
    count = 0;
    while max(jelly.Edges.d_current) > 1.5
        count = count + 1;
        if count >= lim
            lim_reached = 1;
            return
        end
        [~, idx] = sortrows(jelly.Edges, {'d_current'});
        i = idx(length(idx));    
        a = jelly.Edges.EndNodes(i,:);
        node1 = a(1);
        node2 = a(2);
        new_Name = (jelly.Nodes.Node_Names(node1,:) + jelly.Nodes.Node_Names(node2,:))./2;
        %new node is halfway between old nodes
        new_x = (jelly.Nodes.x_coord(node1) + jelly.Nodes.x_coord(node2))./2;
        new_y = (jelly.Nodes.y_coord(node1) + jelly.Nodes.y_coord(node2))./2;
        %Is the new node on the edge? Only if both parents are edges
        if jelly.Nodes.edges(node1) > 0 && jelly.Nodes.edges(node2) > 0  
           if min(jelly.Nodes.edges(node1), jelly.Nodes.edges(node2)) == 1 && max(jelly.Nodes.edges(node1), jelly.Nodes.edges(node2)) == max(jelly.Nodes.edges)
               new_edge = max(jelly.Nodes.edges(node1), jelly.Nodes.edges(node2)) + 1;
           elseif max(jelly.Nodes.edges(node1), jelly.Nodes.edges(node2)) > min(jelly.Nodes.edges(node1), jelly.Nodes.edges(node2)) + 1
               new_edge = 0;
           else
               new_edge = (jelly.Nodes.edges(node1) + jelly.Nodes.edges(node2))/2;
           end
        else
            new_edge = 0;
        end
        
        if (jelly.Nodes.outmus(node1) > 0 || jelly.Nodes.outmus(node2) > 0) && (jelly.Nodes.inmus(node1) == 0 && jelly.Nodes.inmus(node2) == 0)
            if jelly.Nodes.outmus(node1) > 0 && jelly.Nodes.outmus(node2) > 0
                new_outmus = (jelly.Nodes.outmus(node1) + jelly.Nodes.outmus(node2))/2;
            else
                new_outmus = 0;
            end
            new_inmus = 0;
            if mod(jelly.Nodes.mus_num(node1),1) == 0 && jelly.Nodes.mus_num(node1) ~= 0
                new_mus_num = jelly.Nodes.mus_num(node1);
            elseif mod(jelly.Nodes.mus_num(node2),1)== 0 && jelly.Nodes.mus_num(node2)~= 0
                new_mus_num = jelly.Nodes.mus_num(node2);
            end
        elseif (jelly.Nodes.inmus(node1) > 0 || jelly.Nodes.inmus(node2) > 0) && (jelly.Nodes.outmus(node1) == 0 && jelly.Nodes.outmus(node2) == 0)
            if jelly.Nodes.inmus(node1) > 0 && jelly.Nodes.inmus(node2) > 0
                new_inmus = (jelly.Nodes.inmus(node1) + jelly.Nodes.inmus(node2))/2;
            else
                new_inmus = 0;
            end
            new_outmus = 0;
            if  mod(jelly.Nodes.mus_num(node1),1) == 0 && jelly.Nodes.mus_num(node1) ~= 0
                new_mus_num = jelly.Nodes.mus_num(node1);
            elseif mod(jelly.Nodes.mus_num(node2),1)== 0 && jelly.Nodes.mus_num(node2)~= 0
                new_mus_num = jelly.Nodes.mus_num(node2);
            end
        else
            new_inmus = 0;
            new_outmus = 0;
            new_mus_num = 0;
        end
        
        props = table(new_Name, new_x, new_y, [0,0], [0,0], [0,0], [0,0], [0,0], new_outmus, new_inmus, new_mus_num, new_edge, ...
            'VariableNames', {'Node_Names', 'x_coord', 'y_coord', 'velocity', 'F_net', 'F_elastic', 'F_pressure', 'F_muscle', 'outmus', 'inmus', 'mus_num', 'edges'});
        jelly_temp = addnode(jelly_temp, props);
        newid = numnodes(jelly_temp); %Always added to the back, so the new node is the last one.

        %Connect it up to the rest of the graph
        %first its connected to parent nodes
        if jelly.Edges.muscle(i) ~= 0
            m = jelly.Edges.muscle(i); %if parent Edge is muscle, new edges are muscle
        else
            m=0;
        end
        props = table([newid node1; newid node2], [1, 1]', [jelly.Edges.d_current(i)/2 jelly.Edges.d_current(i)/2]', [jelly.Edges.d_rel0(i)/2 jelly.Edges.d_rel0(i)/2]', ...
            [jelly.Edges.d_rel1(i)/2 jelly.Edges.d_rel1(i)/2]', [jelly.Edges.strain0(i) jelly.Edges.strain0(i)]', [jelly.Edges.strain1(i) jelly.Edges.strain1(i)]', [m, m]', ...
            'VariableNames', {'EndNodes', 'Weight', 'd_current', 'd_rel0', 'd_rel1', 'strain0', 'strain1', 'muscle'});
        %% new stuff %%
%         props = table([newid node1; newid node2], [1, 1]', [jelly.Edges.d_current(i)/2 jelly.Edges.d_current(i)/2]', [jelly.Edges.d_rel0(i)/2 jelly.Edges.d_rel0(i)/2]', ...
%             [jelly.Edges.strain0(i) jelly.Edges.strain0(i)]', [m, m]', [jelly.Edges.thickness(i) jelly.Edges.thickness(i)]', ...
%             'VariableNames', {'EndNodes', 'Weight', 'd_current', 'd_rel0', 'strain0', 'muscle', 'thickness'});
        %% end new stuff %%
        jelly_temp = addedge(jelly_temp, props);

        %new node should also connect to any nodes that both parents are
        %connected to. 
        for node_search = 1:numnodes(jelly)
            found = findedge(jelly, node1, node_search);
            if found ~= 0
                found2 = findedge(jelly, node2, node_search);
                if found2 ~= 0
                    dx = jelly.Nodes.x_coord(node_search) - jelly_temp.Nodes.x_coord(newid);
                    dy = jelly.Nodes.y_coord(node_search) - jelly_temp.Nodes.y_coord(newid);
                    d = (dx^2 + dy^2)^(1/2);
                    s0 = (jelly.Edges.strain0(found) + jelly.Edges.strain0(found2))/2;
                    s1 = (jelly.Edges.strain1(found) + jelly.Edges.strain1(found2))/2;
%                     thick = (jelly.Edges.thickness(found) + jelly.Edges.thickness(found2))/2;

                    %%What should the strain on these new edges be?
                    props = table([node_search, newid], 1, d, d/(s0+1), d/(s1+1), s0, s1, 0, ...
                        'VariableNames', {'EndNodes', 'Weight', 'd_current', 'd_rel0', 'd_rel1', 'strain0', 'strain1', 'muscle'});

                    %% new stuff
%                     props = table([node_search, newid], 1, d, d/(s0+1), s0, 0, thick,...
%                         'VariableNames', {'EndNodes', 'Weight', 'd_current', 'd_rel0', 'strain0', 'muscle', 'thickness'});
                    %% end new stuff
                    
                    jelly_temp = addedge(jelly_temp, props); 

                end
            end
        end

        %%Parents are no longer connected to each other
        jelly_temp = rmedge(jelly_temp, node1, node2);
        %plot(jelly_temp, 'XData', jelly_temp.Nodes.x_coord, 'YData', jelly_temp.Nodes.y_coord);

        jelly = jelly_temp;
        if new_edge ~= 0
            %%Reorder the edges
            [trash, idx] = sortrows(jelly.Nodes, {'edges'});
            edge_idx = idx(trash.edges~=0);
            a =  1;
            for i = 1:length(edge_idx)
                jelly.Nodes.edges(edge_idx(i)) = a;
                a = a+1;
            end
            jelly_temp = jelly;
        end
    end
  
%% If two edge nodes get close enough to each other, connect them
    %%How do I make sure this is a concavity?? o.O
    [trash, idx] = sortrows(jelly.Nodes, {'edges'});
    edge_idx = idx(trash.edges~=0);
    new_con = [];
    for i = 1:length(edge_idx)
        node1 = edge_idx(i);
        node2 = edge_idx(mod(i+1, length(edge_idx))+1);

        dx = jelly.Nodes.x_coord(node1) - jelly.Nodes.x_coord(node2);
        dy = jelly.Nodes.y_coord(node1) - jelly.Nodes.y_coord(node2);
        d = (dx^2 + dy^2)^(1/2);
        
        new_con = cat(1, new_con, d);
    end
    
    track = min(new_con);
    
    while track < 1
        for i = 1:length(edge_idx)
            if new_con(i) == track
                break
            end
        end
        node1 = edge_idx(i);
        node_btw = edge_idx(mod(i, length(edge_idx))+1);
        node2 = edge_idx(mod(i+1, length(edge_idx))+1);

        if findedge(jelly, node1, node2) == 0
            
            prev_edge = findedge(jelly, node1, node_btw);
            prev_edge2 = findedge(jelly, node_btw, node2);

            if jelly.Edges.muscle(prev_edge) ~= 0 || jelly.Edges.muscle(prev_edge2) ~= 0
                m = max(jelly.Edges.muscle(prev_edge), jelly.Edges.muscle(prev_edge2)); %if parent Edge is muscle, new edges are muscle
            else
                m=0;
            end

            %Find the strain
            if findedge(jelly, node1, node2) == 0
                s0 = (jelly.Edges.strain0(prev_edge) + jelly.Edges.strain0(prev_edge2))/2;
                s1 = (jelly.Edges.strain1(prev_edge) + jelly.Edges.strain1(prev_edge2))/2;
                %% new stuff %%
    %             thick = (jelly.Edges.thickness(prev_edge) + jelly.Edges.thickness(prev_edge2))/2;
    %             props = table([node1, node2], 1, d, d/(s0+1), s0, m, thick,...
    %                     'VariableNames', {'EndNodes', 'Weight' 'd_current', 'd_rel0', 'strain0', 'muscle', 'thickness'});
                %% end new stuff %%
                props = table([node1, node2], 1, d, d/(s0+1), d/(s1+1), s0, s1, m, ...
                        'VariableNames', {'EndNodes', 'Weight' 'd_current', 'd_rel0', 'd_rel1', 'strain0', 'strain1', 'muscle'});
                jelly = addedge(jelly, props); 
            end

            %%node in between is no longer an edge
            jelly.Nodes.edges(node_btw) = 0;

            [trash, idx] = sortrows(jelly.Nodes, {'edges'});
            edge_idx = idx(trash.edges~=0);
            new_con = [];
            for i = 1:length(edge_idx)
                node1 = edge_idx(i);
                node2 = edge_idx(mod(i+1, length(edge_idx))+1);

                dx = jelly.Nodes.x_coord(node1) - jelly.Nodes.x_coord(node2);
                dy = jelly.Nodes.y_coord(node1) - jelly.Nodes.y_coord(node2);
                d = (dx^2 + dy^2)^(1/2);

                new_con = cat(1, new_con, d);
            end

            track = min(new_con);
        else
            new_con(i) = 1;
            track = min(new_con);
        end

    end

    %%Reorder the edges
    a =  1;
    for i = 1:length(edge_idx)
        jelly.Nodes.edges(edge_idx(i)) = a;
        a = a+1;
    end
    
    %% Maintain the length of the muscle
    new_mus = [];
    while mus_length_cur < 0.9*muscle_length
        [trash, idx] = sortrows(jelly.Nodes, {'edges'});
        edge_idx = idx(trash.edges~=0);
        edge_idx = cat(1, edge_idx, edge_idx(1:2));
        for i = 1:length(edge_idx)-2
            node1 = edge_idx(i);
            node2 = edge_idx(i+1);
            node3 = edge_idx(i+2);
            edge1 = findedge(jelly, node1, node2);
            edge2 = findedge(jelly, node2, node3);
            if jelly.Edges.muscle(edge1) == 1 && jelly.Edges.muscle(edge2) == 0
                new_mus = cat(1, new_mus, edge2);
            end
            if jelly.Edges.muscle(edge2) == 1 && jelly.Edges.muscle(edge1) == 0
                new_mus = cat(1, new_mus, edge1);
            end
        end
        for i = 1:length(new_mus)
            jelly.Edges.muscle(new_mus(i)) = 1;
        end
        mus_length_cur = sum(jelly.Edges.d_current(jelly.Edges.muscle == 1));
    end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        
