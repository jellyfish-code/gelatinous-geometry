%{
======================================================================
    Main function evolving jellyfish butterfly mesh.  
======================================================================

INPUT:
        elast0 (float):             Elasticity of spring (in Pascals).
        elast1 (float):             Elasticity of spring (in Pascals).
        vis (float):                Viscosity of dashpot (in Pascal*seconds).
        bulk_modulus (float):       Bulk modulus of jellyfish (in Pascals).
        area0 (float):              Initial area of jellyfish (in ?).
        muscle_strain (float):      Strain of jellyfish (in ?).
        contraction_rate (float):   Number of jellyfish contractions per minute.
        max_dR (float):             Maximum change in radius during contraction.
        dR_rate (float):            Increase in radius change with distance from anchored end (Figure S6c in paper).
        folder_save (string):       Name of folder in which simulation images is to be saved.
        datapath (string):          Directory in which folder specificed in folder_save can be found.

OUTPUT:
        None. 
%}

function visco_butterfly_SLM_newmus(elast0, elast1, vis, bulk_modulus, area0, muscle_strain, contraction_rate, max_dR, dR_rate, folder_save, datapath)
%% Set up parameters, everything in Pa(N/m^2) and s
    %Measured parameters
    contraction_duration = 0.8; %s
    contraction_strength = (elast0+elast1)*muscle_strain; %Pa
    
    %Graft parameters
    relax_duration = (60-contraction_rate*contraction_duration)/(contraction_rate + 1); %seconds
    vel = [];
    if contraction_rate > 75
        contraction_duration = 60/contraction_rate;
    end
    if relax_duration < 0
        relax_duration = 0;
    end

    %time
    time_step = 15; %minutes
    time_end = 2000; %hours
    time_steps = time_end*60/time_step;
    
    relax_param = (1-exp(-1*elast1/vis * (time_step * 60)));

    wing_width = 11;
    body_height = 2;
    body_start = [7, 6, 6, 7, 9];
    body_end = [9, 11, 11, 12, 11];
    
    %% Define the centers of all of the pieces
    wing1 = [6, 1]; %lower
    wing2 = [6, 1+wing_width];
    body = [6 + body_height, 1+wing_width/2];

    %% Writing images
    path0 = datapath;
    %     folder_save = ['KV_052620_off', num2str(offset), 'rate', num2str(contraction_rate)];

    if ~isempty(path0)
        Dr = dir([path0 '/' folder_save]);
        if isempty(Dr)
            S = mkdir(path0,folder_save);
            if ~S, disp('Fail to make folder!'); return;  
            end
        else
        end                        
        path1 =[path0, '/', folder_save];         % the place to store images
    else
        path0 = pwd; 
    end

    %% Set up jellyfish geometry as an array
    %This creates the initial offset array, assigning coordinates for each
    %node and connecting nodes. It also designates muscle and non-muscle
    %nodes.
    
    %% Start with uncut jellyfish, find the balance of elastic and pressure forces
    %%This function sets the "relaxed area" at 105% of the initial area,
    %%and allows the jellyfish F_elastic and F_pressure to find equilibrium
    %%over 120 hours
    [jelly_eq, ~, ~] = equilibrium_initial_KV(elast0, vis, bulk_modulus, 100, 1200, area0, contraction_strength);

    %% Make a matrix with the right shape 
    [jelly_initial, row_start, row_end] = butterfly_mesh();
    [muscle_top_outer, muscle_top_inner, muscle_lowleft_outer,  muscle_lowleft_inner, muscle_lowright_outer, muscle_lowright_inner] = butterfly_muscle();
    a = size(jelly_initial);
    j_eq_center = jelly_eq(6,6,:);
    j_temp = zeros(a);
    %cut jelly into butterfly shape
    %upper wing
    jelly_initial(1:6,1:11, 1) = jelly_eq(6:11, :, 1) - j_eq_center(1)+wing1(1);
    jelly_initial(1:6,1:11, 2) = jelly_eq(6:11, :, 2) - j_eq_center(2)+wing1(2);

    %lower wing
    jelly_initial(a(1)-5:a(1), 7:a(2), 1) = jelly_eq(1:6, :, 1) - j_eq_center(1)+wing2(1);
    jelly_initial(a(1)-5:a(1), 7:a(2), 2) = jelly_eq(1:6, :, 2) - j_eq_center(2)+wing2(2);

    %center
    for i = 1:length(body_start)
        j_temp(4+i, body_start(i):body_end(i), 1) = jelly_eq(3+i, body_start(i)-5:body_end(i)-5, 1) - j_eq_center(1)+body(1);
        j_temp(4+i, body_start(i):body_end(i), 2) = jelly_eq(3+i, body_start(i)-5:body_end(i)-5, 2)- j_eq_center(2)+body(2);
    end

    %% integrate the pieces
    for i = 1:length(body_start)
        jelly_initial(4+i,body_start(i):body_end(i),1) = (jelly_initial(4+i,body_start(i):body_end(i),1) + j_temp(4+i,body_start(i):body_end(i),1))/2;
        jelly_initial(4+i,body_start(i):body_end(i),2) = (jelly_initial(4+i,body_start(i):body_end(i),2) + j_temp(4+i,body_start(i):body_end(i),2))/2;
    end
    %% Convert array to graph
    jelly = convert_jelly_graph(jelly_initial, row_start, row_end);
    %%Add additional parameters
    jelly.Nodes.velocity = zeros(numnodes(jelly), 2);
    jelly.Nodes.F_net = zeros(numnodes(jelly), 2);
    jelly.Nodes.F_elastic = jelly.Nodes.F_net;
    jelly.Nodes.F_pressure = jelly.Nodes.F_net;
    jelly.Nodes.F_muscle = jelly.Nodes.F_net;
    jelly.Nodes.outmus = zeros(numnodes(jelly),1);
    jelly.Nodes.inmus = zeros(numnodes(jelly),1);
    jelly.Nodes.mus_num = zeros(numnodes(jelly),1);

    %% Define nodes that make up the edges (in order)
    edges = find_edges_butterfly(row_start, row_end);
    j_edges = zeros(numnodes(jelly),1);
    c = 1;
    for i = 1:length(edges)-1
        for j = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(j,1) == edges(1, i) && jelly.Nodes.Node_Names(j,2) == edges(2, i)
                j_edges(j) = c;
                c = c+1;
            end
        end
    end
    jelly.Nodes.edges = j_edges;

    %%Define a bunch of variables in the graph
    [jelly_i, row_start, row_end] = butterfly_mesh(); 
    jelly_but_i = convert_jelly_graph(jelly_i, row_start, row_end);
    jelly_but_i.Nodes.edges = j_edges;
    [j_area, ~] = area(jelly_but_i);
    area_relax = area0*j_area;
    jelly.Edges.d_rel0 = jelly_but_i.Edges.d_current; %This is a weird one. The relaxed length is the length
    %before equilibrium is found. So I'm just initializing another graft
    jelly.Edges.d_rel1 = jelly.Edges.d_rel0;
    jelly.Edges.strain0 = (jelly.Edges.d_current - jelly.Edges.d_rel0)./jelly.Edges.d_rel0;
    jelly.Edges.strain1 = (jelly.Edges.d_current - jelly.Edges.d_rel1)./jelly.Edges.d_rel1;

    m1 = zeros(length(muscle_top_outer),1);
    m2 = zeros(length(muscle_lowleft_outer),1);
    m3 = zeros(length(muscle_lowright_outer),1);
    m4 = zeros(length(muscle_top_inner),1);
    m5 = zeros(length(muscle_lowleft_inner),1);
    m6 = zeros(length(muscle_lowright_inner),1);

    outcount = 1;
    incount = 1;
    %% Define edges that make up the muscles
    for j = 1:length(muscle_top_outer)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_top_outer(2,j) && jelly.Nodes.Node_Names(i,2) == muscle_top_outer(1,j)
                m1(j,1) = i;
                jelly.Nodes.outmus(i) = outcount;
                outcount = outcount+1;
            end
        end
    end

    for j = 1:length(muscle_top_outer)-1
        a = findedge(jelly, m1(j), m1(j+1));
        jelly.Edges.muscle(a) = 1;
    end
    for j = 1:length(muscle_lowleft_outer)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_lowright_outer(2,j) && jelly.Nodes.Node_Names(i,2) == muscle_lowright_outer(1,j)
                m3(j) = i;
                jelly.Nodes.outmus(i) = outcount;
                outcount = outcount+1;
            end
        end
    end

    for j = 1:length(muscle_lowleft_outer)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_lowleft_outer(2,j) && jelly.Nodes.Node_Names(i,2) == muscle_lowleft_outer(1,j)
                m2(j) = i;
                jelly.Nodes.outmus(i) = outcount;
                outcount = outcount+1;
            end
        end
    end
    for j = 1:length(muscle_lowleft_outer)-1
        a = findedge(jelly, m2(j), m2(j+1));
        jelly.Edges.muscle(a) = 1;
        b = findedge(jelly, m3(j), m3(j+1));
        jelly.Edges.muscle(b) = 1;
    end

    for j = 1:length(muscle_top_inner)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_top_inner(2,j) && jelly.Nodes.Node_Names(i,2) == muscle_top_inner(1,j)
                m4(j) = i;
                jelly.Nodes.inmus(i) = incount;
                incount = incount+1;
            end
        end
    end
    for j = 1:length(muscle_top_inner)-1
        a = findedge(jelly, m4(j), m4(j+1));
        jelly.Edges.muscle(a) = 0.5;
    end
    for j = 1:length(muscle_lowleft_inner)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_lowright_inner(2,j) && jelly.Nodes.Node_Names(i,2) == muscle_lowright_inner(1,j)
                m6(j) = i;
                jelly.Nodes.inmus(i) = incount;
                incount = incount+1;
            end
        end
    end
    
    for j = 1:length(muscle_lowleft_inner)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_lowleft_inner(2,j) && jelly.Nodes.Node_Names(i,2) == muscle_lowleft_inner(1,j)
                m5(j) = i;
                jelly.Nodes.inmus(i) = incount;
                incount = incount+1;
            end
        end
    end

    for j = 1:length(muscle_lowleft_inner)-1
        a = findedge(jelly, m5(j), m5(j+1));
        jelly.Edges.muscle(a) = 0.5;
        b = findedge(jelly, m6(j), m6(j+1));
        jelly.Edges.muscle(b) = 0.5;
    end
    
    jelly.Nodes.mus_num(m1(1:5)) = 1;   jelly.Nodes.mus_num(m4(1:5)) = 1;
    jelly.Nodes.mus_num(m1(6)) = 1.5;   jelly.Nodes.mus_num(m4(6)) = 1.5;
    jelly.Nodes.mus_num(m1(7)) = 2;     jelly.Nodes.mus_num(m4(7:9)) = 2;
    jelly.Nodes.mus_num(m1(8)) = 2.5;   jelly.Nodes.mus_num(m4(10)) = 2.5;
    jelly.Nodes.mus_num(m1(9:13)) = 3;  jelly.Nodes.mus_num(m4(11:15)) = 3;
    jelly.Nodes.mus_num(m3) = 4;        jelly.Nodes.mus_num(m6) = 4;
    jelly.Nodes.mus_num(m2) = 5;        jelly.Nodes.mus_num(m5) = 5;
    
    muscle_length = sum(jelly.Edges.d_current(jelly.Edges.muscle == 1));
    [jelly, ~] = remesh_SLM_newmus(jelly, muscle_length);

    %% Image the graph
    if any(isfinite(jelly.Nodes.x_coord)-1) == 1 || any(isfinite(jelly.Nodes.y_coord)-1) == 1
        return
    end
    
    if max(jelly.Edges.d_current) > 10
        return
    end
    
    figure1 = plot(jelly, 'XData', jelly.Nodes.x_coord, 'YData', jelly.Nodes.y_coord, 'EdgeCData', jelly.Edges.strain0, 'LineWidth', 1, 'NodeLabel', {});
    %hold on
    %figure1 = quiver(jelly.Nodes.x_coord, jelly.Nodes.y_coord, jelly.Nodes.F_net(:,1), jelly.Nodes.F_net(:,2));

    hold off
    xlim([-1,14]);
    ylim([-1,14]);

    cd(path1);                                                            % write the image data
    saveas(figure1, '000.jpg')
    cd(path0);     
    pause(0.001)

   %%Setup takes about 4 seconds%%%
    %% Start the sim
    for time = 1:time_steps

        %Print the current time_step (in hours)
        hours = time/(60/time_step);
        
        %if hours == 1000
            
        %    pause(0.001)
        %    a_rat = area(jelly)/area_relax;
        %    jelly = sequential_reorg(jelly);
        %    area_relax = area(jelly)/a_rat;
        %end
            

        %% Update the current length of edges
        for i = 1:numedges(jelly)
            [node1, node2] = findedge(jelly, i);
            dx = jelly.Nodes.x_coord(node1) - jelly.Nodes.x_coord(node2);
            dy = jelly.Nodes.y_coord(node1) - jelly.Nodes.y_coord(node2);
            dist_current = (dx^2 + dy^2)^(1/2);
            jelly.Edges.d_current(i) = dist_current;
        end
        
        jelly.Edges.d_rel1 = -1*(jelly.Edges.d_current.*jelly.Edges.d_rel1)./((jelly.Edges.d_rel1 - jelly.Edges.d_current).*relax_param - jelly.Edges.d_rel1);
        if max(jelly.Edges.d_current) > 10 || max(jelly.Edges.strain0) > 1 %|| min(jelly.Edges.strain0) < -0.5
            cd(path1)
            writematrix(vel, 'velocity.xlsx');
            cd(path0)
         
            return
        end       
        %% Calculate new muscle coordinates from contraction
        %muscles are "synchronized", calculations for all muscle bands happen
        %simultaneously.
        [jelly, done] = contraction5(jelly, contraction_strength, muscle_strain, max_dR, dR_rate);
        if done == 0
            cd(path1)
            writematrix(vel, 'velocity.xlsx');
            cd(path0)
            return
        end
        
        jelly = SLM_elastic(jelly, elast0, elast1);

        %% pressure force
        jelly.Nodes.F_pressure = find_f_pressure(jelly, area_relax, bulk_modulus);  
        F_contract = jelly.Nodes.F_elastic + jelly.Nodes.F_pressure + jelly.Nodes.F_muscle;
        F_relax = jelly.Nodes.F_elastic + jelly.Nodes.F_pressure;
        
        %% Instead of separating out by contraction and relaxation phases, we are just finding the average F_net over the time step
        jelly.Nodes.F_net = (F_contract*contraction_rate*contraction_duration + F_relax*relax_duration*(contraction_rate+1))/60;

        %% Update the position of each node
        jelly.Nodes.velocity = jelly.Nodes.F_net./vis;
        contraction_displacement = jelly.Nodes.velocity*time_step*60;
        jelly.Nodes.x_coord = jelly.Nodes.x_coord + contraction_displacement(:,1);
        jelly.Nodes.y_coord = jelly.Nodes.y_coord + contraction_displacement(:,2);

       
        %% Find the current aspect ratio
        if any(isfinite(jelly.Nodes.x_coord)-1) == 1 || any(isfinite(jelly.Nodes.y_coord)-1) == 1
            cd(path1);                                                            % write the image data
            writematrix(vel, 'velocity.xlsx');
            cd(path0); 
            return
        end
       

        %% Remesh every 10 hours
        if mod(hours, 5) == 0 || min(jelly.Edges.d_current) < 0.35
            [jelly, lim_reached] = remesh_SLM_butterfly(jelly, muscle_length);
            if lim_reached == 1
                cd(path1);                                                            % write the image data
                writematrix(vel, 'velocity.xlsx');
                cd(path0); 
                return
            end
        end
    %     
        if mod(hours, 20) == 0
            %% image new relaxed jelly every 2 hours
            figure1 = plot(jelly, 'XData', jelly.Nodes.x_coord, 'YData', jelly.Nodes.y_coord, 'EdgeCData', jelly.Edges.strain0 + jelly.Edges.strain1, 'LineWidth', 1, 'NodeLabel', {});
            %hold on
            %figure1 = quiver(jelly.Nodes.x_coord, jelly.Nodes.y_coord, jelly.Nodes.F_net(:,1), jelly.Nodes.F_net(:,2));
            hold off
            xlim([-1,14]);
            ylim([-1,14]);
            %% Save images
            cd(path1);                                                            % write the image data
            i = floor(hours)*10;
            saveas(figure1, [num2str(i) '.jpg'])
            cd(path0);     
            pause(0.001)
            
        end
        if mod(hours, 5) == 0
            [~, edge_idx] = area(jelly);
            vel_max = max((jelly.Nodes.velocity(edge_idx,1).^2) + jelly.Nodes.velocity(edge_idx,2).^2).^(1/2);
            vel = cat(1, vel, [vel_max, hours]);
        end  

    end

    %% final image
    figure1 = plot(jelly, 'XData', jelly.Nodes.x_coord, 'YData', jelly.Nodes.y_coord);
    hold on
    figure1 = quiver(jelly.Nodes.x_coord, jelly.Nodes.y_coord, jelly.Nodes.F_net(:,1), jelly.Nodes.F_net(:,2));

    hold off
    xlim([-1,14]);
    ylim([-1,14]);
    
    %% Plot aspect ratio
    cd(path1);                                                            % write the image data
    writematrix(vel, 'velocity.xlsx');
    cd(path0);  
%  