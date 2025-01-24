%{
======================================================================
    Function to remesh jellyfish 
    Remeshing is required as biological materials can undergo large
    deformations, which cannot be captured by static meshes. [See Chen 2000?]
======================================================================

INPUT:
        jelly (graph):          Graph containing information on Nodes and Edges of jellyfish.                     
        muscle_length (double): Sum of lengths of muscles in current geometry. Used to preserve total length of muscles during remeshing. 

OUTPUT:
        jelly (graph):          Updated graph containing information on Nodes and Edges of jellyfish.                     
        lim_reached (double):   Used to prevent infinite loop of remeshing. If returned 1, simulation exits (check?).
%}

%% Remesh
function [jelly, lim_reached] = remesh_SLM_butterfly(jelly, muscle_length)
lim_reached = 0;
jelly_temp = jelly;
mus_length_cur = sum(jelly.Edges.d_current(jelly.Edges.muscle == 1));
lim = numedges(jelly);
count = 0;
    %% Delete nodes when needed
    while min(jelly.Edges.d_current) < 0.35
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
            'VariableNames', {'Node_Names', 'x_coord', 'y_coord', 'velocity', 'F_net', 'stress_elastic', 'pressure', 'stress_muscle', 'outmus', 'inmus', 'mus_num', 'edges'});
        jelly_temp = addnode(jelly_temp, props);
        newid = numnodes(jelly_temp); %Always added to the back, so the new node is the last one.

        %%Connect it up to connections of parents
        for node_search = 1:numnodes(jelly)
            found = findedge(jelly, node1, node_search);
            found2 = findedge(jelly, node2, node_search);
            if length(found) >1 || length(found2)>1
                lim_reached = 1;
                return
            end
            if found ~= 0 || found2 ~= 0                
                dx = jelly.Nodes.x_coord(node_search) - jelly_temp.Nodes.x_coord(newid);
                dy = jelly.Nodes.y_coord(node_search) - jelly_temp.Nodes.y_coord(newid);
                d = (dx^2 + dy^2)^(1/2);
                if found ~= 0 && found2 ~= 0
                    s0 = (jelly.Edges.strain0(found) + jelly.Edges.strain0(found2))/2;
                    s1 = (jelly.Edges.strain1(found) + jelly.Edges.strain1(found2))/2;
                    %s_viscous = (jelly.Edges.strainviscous(found) + jelly.Edges.strainviscous(found2))/2;

                elseif found ~= 0
                    s0 = jelly.Edges.strain0(found);
                    s1 = jelly.Edges.strain1(found);
                    %s_viscous = jelly.Edges.strainviscous(found);

                else
                    s0 = jelly.Edges.strain0(found2);
                    s1 = jelly.Edges.strain1(found2);
                    %s_viscous = jelly.Edges.strainviscous(found2);

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
            'VariableNames', {'Node_Names', 'x_coord', 'y_coord', 'velocity', 'F_net', 'stress_elastic', 'pressure', 'stress_muscle', 'outmus', 'inmus', 'mus_num', 'edges'});
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
                    s0 = (jelly.Edges.strain0(found(1)) + jelly.Edges.strain0(found2(1)))/2;
                    s1 = (jelly.Edges.strain1(found(1)) + jelly.Edges.strain1(found2(1)))/2;
                    %s_viscous = (jelly.Edges.strainviscous(found) + jelly.Edges.strainviscous(found2))/2;

                    %%What should the strain on these new edges be?
                    props = table([node_search, newid], 1, d, d/(s0+1), d/(s1+1), s0, s1, 0, ...
                        'VariableNames', {'EndNodes', 'Weight', 'd_current', 'd_rel0', 'd_rel1', 'strain0', 'strain1', 'muscle'});

                    
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
  
    %% Remesh diamonds that get squashed
    %If an edge has a really high strain, check to see if diamond is
    %getting squashed
    for i = 1:numedges(jelly)
        %if jelly.Edges.strain0(i) > 0.15 
            a = jelly_temp.Edges.EndNodes(i,:);
            squashed = [];
            %is it a diamond? 
            co_neighs = intersect(neighbors(jelly_temp, a(1)), neighbors(jelly_temp, a(2)));
            if length(co_neighs) == 2 %can there be more than 2?
                for j = 1:length(co_neighs)
                    l1 = jelly_temp.Edges.d_current(findedge(jelly_temp, a(1), co_neighs(j))) + jelly_temp.Edges.d_current(findedge(jelly_temp, a(2), co_neighs(j)));
                    if length(l1) > 1
                        lim_reached = 1;
                        return
                    end
                    squashed = cat(1, squashed, jelly_temp.Edges.d_current(i)/l1);
                end
                %is it squashed?
                if max(squashed) > 0.9
                    %connect opposite corners
                    d = ((jelly.Nodes.x_coord(co_neighs(1)) - jelly.Nodes.x_coord(co_neighs(2)))^2 + (jelly.Nodes.y_coord(co_neighs(1)) - jelly.Nodes.y_coord(co_neighs(2)))^2)^(1/2);
                    s0 = mean([jelly_temp.Edges.strain0(findedge(jelly_temp, a(1), co_neighs(squashed==min(squashed)))), jelly_temp.Edges.strain0(findedge(jelly_temp, a(2), co_neighs(squashed==min(squashed))))]);
                    s1 = mean([jelly_temp.Edges.strain1(findedge(jelly_temp, a(1), co_neighs(squashed==min(squashed)))), jelly_temp.Edges.strain1(findedge(jelly_temp, a(2), co_neighs(squashed==min(squashed))))]);
                    %s_viscous = mean([jelly_temp.Edges.strainviscous(findedge(jelly_temp, a(1), co_neighs(squashed==min(squashed)))), jelly_temp.Edges.strainviscous(findedge(jelly_temp, a(2), co_neighs(squashed==min(squashed))))]);

                    props = table([co_neighs(1), co_neighs(2)], 1, d, d/(s0+1), d/(s1+1), s0, s1, 0, ...
                    'VariableNames', {'EndNodes', 'Weight', 'd_current', 'd_rel0', 'd_rel1', 'strain0', 'strain1', 'muscle'});

                    jelly_temp = addedge(jelly_temp, props); 

                    %update strain and d_relax on diamond edges
                    jelly_temp.Edges.strain0(findedge(jelly_temp, a(1), co_neighs(squashed==max(squashed)))) = mean([jelly_temp.Edges.strain0(findedge(jelly_temp, a(1), co_neighs(squashed==max(squashed)))), jelly_temp.Edges.strain0(findedge(jelly_temp, a(1), a(2)))]);
                    jelly_temp.Edges.strain0(findedge(jelly_temp, a(2), co_neighs(squashed==max(squashed)))) = mean([jelly_temp.Edges.strain0(findedge(jelly_temp, a(2), co_neighs(squashed==max(squashed)))), jelly_temp.Edges.strain0(findedge(jelly_temp, a(1), a(2)))]);

                    jelly_temp.Edges.strain1(findedge(jelly_temp, a(1), co_neighs(squashed==max(squashed)))) = mean([jelly_temp.Edges.strain1(findedge(jelly_temp, a(1), co_neighs(squashed==max(squashed)))), jelly_temp.Edges.strain1(findedge(jelly_temp, a(1), a(2)))]);
                    jelly_temp.Edges.strain1(findedge(jelly_temp, a(2), co_neighs(squashed==max(squashed)))) = mean([jelly_temp.Edges.strain1(findedge(jelly_temp, a(2), co_neighs(squashed==max(squashed)))), jelly_temp.Edges.strain1(findedge(jelly_temp, a(1), a(2)))]);

                    jelly_temp.Edges.d_rel0(findedge(jelly_temp, a(1), co_neighs(squashed==max(squashed)))) = jelly_temp.Edges.d_current(findedge(jelly_temp, a(1), co_neighs(squashed==max(squashed))))/(jelly_temp.Edges.strain0(findedge(jelly_temp, a(1), co_neighs(squashed==max(squashed)))) + 1);
                    jelly_temp.Edges.d_rel0(findedge(jelly_temp, a(2), co_neighs(squashed==max(squashed)))) = jelly_temp.Edges.d_current(findedge(jelly_temp, a(2), co_neighs(squashed==max(squashed))))/(jelly_temp.Edges.strain0(findedge(jelly_temp, a(2), co_neighs(squashed==max(squashed)))) + 1);
                    
                    jelly_temp.Edges.d_rel1(findedge(jelly_temp, a(1), co_neighs(squashed==max(squashed)))) = jelly_temp.Edges.d_current(findedge(jelly_temp, a(1), co_neighs(squashed==max(squashed))))/(jelly_temp.Edges.strain1(findedge(jelly_temp, a(1), co_neighs(squashed==max(squashed)))) + 1);
                    jelly_temp.Edges.d_rel1(findedge(jelly_temp, a(2), co_neighs(squashed==max(squashed)))) = jelly_temp.Edges.d_current(findedge(jelly_temp, a(2), co_neighs(squashed==max(squashed))))/(jelly_temp.Edges.strain1(findedge(jelly_temp, a(2), co_neighs(squashed==max(squashed)))) + 1);
                    
                    if jelly.Edges.muscle(i) ~= 0
                        jelly_temp.Edges.muscle(findedge(jelly_temp, a(1), co_neighs(squashed==max(squashed)))) = jelly.Edges.muscle(i);
                        jelly_temp.Edges.muscle(findedge(jelly_temp, a(2), co_neighs(squashed==max(squashed)))) = jelly.Edges.muscle(i);
                        
%                         jelly_temp.Nodes.inmus(co_neighs(squashed==max(squashed))) = mean([jelly.Nodes.inmus(a(1)), jelly.Nodes.inmus(a(2))]);
%                         jelly_temp.Nodes.outmus(co_neighs(squashed==max(squashed))) = mean([jelly.Nodes.outmus(a(1)), jelly.Nodes.outmus(a(2))]);
                    end
                    %%Now kill the old edge
                    jelly_temp = rmedge(jelly_temp, a(1), a(2));

                end
            end
        %end
    end
    jelly = jelly_temp;    

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

        if findedge(jelly, node1, node2) == 0 && length(intersect(neighbors(jelly, node1), neighbors(jelly, node2))) == 1
            
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
                %s_viscous = (jelly.Edges.strainviscous(prev_edge) + jelly.Edges.strainviscous(prev_edge2))/2;

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
            
            

