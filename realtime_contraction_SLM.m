function realtime_contraction_SLM(elast0, elast1, vis, bulk_modulus, muscle_strain, folder_save)
%% Set up parameters, everything in Pa(N/m^2) and s
%elast0 = 0.05;
%elast1 = 0.05;
%vis = 750; %only used during set up
damping = 0.1; 
%bulk_modulus = 0.75;
area0 = 1.001;
%muscle_strain = 1.4;
contraction_rate = 20;
mass = 0.01; %grams per node
mus_length = []; %track muscle length over time

graft_type = 'offset';
offset = 2;
%folder_save = 'live_contraction_butterfly1';

%Measured parameters
contraction_duration = 0.8; %s
contraction_strength = (elast0+elast1)*muscle_strain; %Pa
%Graft parameters
relax_duration = (60-contraction_rate*contraction_duration)/(contraction_rate + 1); %seconds

if contraction_rate > 75
    contraction_duration = 60/contraction_rate;
end
if relax_duration < 0
    relax_duration = 0;
end

%time
time_step = 0.1; %seconds
time_end = 10; %seconds
time_steps = time_end/time_step;

%% Writing images
path0 = '/central/home/mgong/Documents/Model';
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
[jelly_eq, ~, d_uncut] = equilibrium_initial_KV(elast0, vis, bulk_modulus, 0, 1, area0, contraction_strength);

if contains(graft_type, 'butterfly') 
    wing_width = 11;
    body_height = 2;
    body_start = [7, 6, 6, 7, 9];
    body_end = [9, 11, 11, 12, 11];
    
    %% Define the centers of all of the pieces
    wing1 = [6, 1]; %lower
    wing2 = [6, 1+wing_width];
    body = [6 + body_height, 1+wing_width/2];

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

    muscle_length = sum(jelly.Edges.d_current(jelly.Edges.muscle == 1));
    [jelly, ~] = remesh_SLM(jelly, muscle_length);
    xlim([0,12]);
    ylim([-1,14]);
    if contains(graft_type, 'disconnected')
        %% keep jelly not quite connected
        %create copies of nodes that are at junctions
        %for butterfly, nodes 46, 58, 51, 63
        jelly_temp = jelly;
        new = numnodes(jelly);
        props = vertcat(jelly.Nodes(46,:), jelly.Nodes(58,:), jelly.Nodes(51,:), jelly.Nodes(63,:));
        jelly_temp = addnode(jelly_temp, props);
        %these are for the inside piece
        %connect new nodes to inside piece
        inside1 = findedge(jelly, 46, 52);
        inside2 = findedge(jelly, 46, 53);
        inside3 = findedge(jelly, 46, 47);       
        inside4 = findedge(jelly, 52, 58);
        inside5 = findedge(jelly, 53, 58);
        inside6 = findedge(jelly, 59, 58);
        inside7 = findedge(jelly, 57, 51);
        inside8 = findedge(jelly, 50, 51);
        inside9 = findedge(jelly, 57, 63);
        inside0 = findedge(jelly, 62, 63);
        props = vertcat(jelly.Edges(inside1,:), jelly.Edges(inside2,:), jelly.Edges(inside3,:), jelly.Edges(inside4,:), jelly.Edges(inside5,:), ...
            jelly.Edges(inside6,:), jelly.Edges(inside7,:), jelly.Edges(inside8,:), jelly.Edges(inside9,:), jelly.Edges(inside0,:));
        props.EndNodes = [new+1, 52; new+1, 53; new+1, 47; new+2, 52; new+2, 53; new+2, 59; new+3, 57; new+3, 50; new+4, 57; new+4, 62];
        jelly_temp = addedge(jelly_temp, props);
        
        %connect some more stuff up
        d1 = ((jelly.Nodes.x_coord(41) - jelly.Nodes.x_coord(46))^2 + (jelly.Nodes.y_coord(41) - jelly.Nodes.y_coord(46))^2)^(1/2);
        d2 = ((jelly.Nodes.x_coord(43) - jelly.Nodes.x_coord(51))^2 + (jelly.Nodes.y_coord(43) - jelly.Nodes.y_coord(51))^2)^(1/2);
        props = table([46,41; 58,66; 43,51; 63,68; new+1, 41; new+2, 66; new+3, 43; new+4, 68], [1; 1; 1; 1; 1; 1; 1; 1], ...
            [d1; d1; d2; d2; d1; d1; d2; d2], [d1; d1; d2; d2; d1; d1; d2; d2], [d1; d1; d2; d2; d1; d1; d2; d2], [0;0;0;0;0;0;0;0],...
            [0;0;0;0;0;0;0;0], [0;0;0;0;0;0;0;0],'VariableNames', {'EndNodes', 'Weight', 'd_current', 'd_rel0', 'd_rel1', 'strain0', 'strain1', 'muscle'});
        jelly_temp = addedge(jelly_temp, props);
        
        %outside piece no longer connect to inside piece
        jelly_temp = rmedge(jelly_temp, 46, 52);
        jelly_temp = rmedge(jelly_temp, 46, 53);
        jelly_temp = rmedge(jelly_temp, 46, 47);
        jelly_temp = rmedge(jelly_temp, 58, 52);
        jelly_temp = rmedge(jelly_temp, 58, 53);
        jelly_temp = rmedge(jelly_temp, 58, 59);
        jelly_temp = rmedge(jelly_temp, 57, 51);
        jelly_temp = rmedge(jelly_temp, 50, 51);
        jelly_temp = rmedge(jelly_temp, 57, 63);
        jelly_temp = rmedge(jelly_temp, 62, 63);

        %edges from 40-47, 65-59, 62-69, 50-44 are no longer connected.
        jelly_temp = rmedge(jelly_temp, 40, 47);
        jelly_temp = rmedge(jelly_temp, 59, 65);
        jelly_temp = rmedge(jelly_temp, 62, 69);
        jelly_temp = rmedge(jelly_temp, 44, 50);

        jelly = jelly_temp;
    end
elseif contains(graft_type, 'offset')
    %% Make a matrix with the right shape for offset graft
    [jelly_initial, row_start, row_end] = offset_mesh(offset); 
    [muscle_outer, muscle_inner] = offset_muscle(offset);
    a = size(jelly_initial);

    %% Then cut the jellyfish (update the x_ and y_coordinates) 
    jelly_initial(6:a(1), 1: a(1), :) = jelly_eq(6:a(1), 1: a(1), :);
    jelly_initial(1:6, 1+offset: a(1)+offset, 1) = jelly_eq(1:6, 1: a(1), 1) + offset*d_uncut;
    jelly_initial(1:6, 1+offset: a(1)+offset, 2) = jelly_eq(1:6, 1: a(1), 2);
    jelly_initial(6, 1+offset:6+offset, 1) = jelly_eq(6, 1:6, 1) + offset*d_uncut;
    jelly_initial(6, 7:a(1), 1) = jelly_eq(6, 7:a(1), 1);

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

    %% Define nodes that make up the edges (in order)
    edges = find_edges(muscle_outer);
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
    m1 = zeros(length(muscle_outer),1);
    m2 = zeros(length(muscle_outer),1);
    m3 = zeros(length(muscle_inner),1);
    m4 = zeros(length(muscle_inner),1);

    %%More parameters
    [jelly_i, row_start, row_end] = offset_mesh(offset); 
    jelly_off_i = convert_jelly_graph(jelly_i, row_start, row_end);
    jelly_off_i.Nodes.edges = j_edges;
    [j_area, ~] = area(jelly_off_i);
    area_relax = area0*j_area;
    jelly.Edges.d_rel0 = jelly_off_i.Edges.d_current; %This is a weird one. The relaxed length is the length
    %before equilibrium is found. So I'm just initializing another offset graft
    jelly.Edges.d_rel1 = jelly.Edges.d_current; %This is assumed to be fully relaxed
    jelly.Edges.strain0 = (jelly.Edges.d_current - jelly.Edges.d_rel0)./jelly.Edges.d_rel0;
    jelly.Edges.strain1 = (jelly.Edges.d_current - jelly.Edges.d_rel1)./jelly.Edges.d_rel1;
    
    outcount = 1;
    incount = 1;
    %% Define edges that make up the muscles
    for j = 1:length(muscle_outer)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_outer(1,j,1) && jelly.Nodes.Node_Names(i,2) == muscle_outer(2,j,1)
                m1(j) = i;
                jelly.Nodes.outmus(i) = outcount;
                outcount = outcount+1;
            end
        end
    end
    for j = 1:length(muscle_outer)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_outer(1,j,2) && jelly.Nodes.Node_Names(i,2) == muscle_outer(2,j,2)
                m2(j) = i;
                jelly.Nodes.outmus(i) = outcount;
                outcount = outcount+1;
            end
        end
    end
    for j = 1:length(muscle_outer)-1
        a = findedge(jelly, m1(j), m1(j+1));
        jelly.Edges.muscle(a) = 1;
        b = findedge(jelly, m2(j), m2(j+1));
        jelly.Edges.muscle(b) = 1;
    end
    for j = 1:length(muscle_inner)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_inner(1,j,1) && jelly.Nodes.Node_Names(i,2) == muscle_inner(2,j,1)
                m3(j) = i;
                jelly.Nodes.inmus(i) = incount;
                incount = incount+1;
            end
        end
    end
    for j = 1:length(muscle_inner)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_inner(1,j,2) && jelly.Nodes.Node_Names(i,2) == muscle_inner(2,j,2)
                m4(j) = i;
                jelly.Nodes.inmus(i) = incount;
                incount = incount+1;
            end
        end
    end
    for j = 1:length(muscle_inner)-1
        a = findedge(jelly, m3(j), m3(j+1));
        jelly.Edges.muscle(a) = 0.5;
        b = findedge(jelly, m4(j), m4(j+1));
        jelly.Edges.muscle(b) = 0.5;
    end

    muscle_length = sum(jelly.Edges.d_current(jelly.Edges.muscle == 1));
    [jelly, ~] = remesh_SLM(jelly, muscle_length);
            
    if max(jelly.Edges.d_current) > 10
        return
    end
    xlim([0, 11 + offset]);
    ylim([-1, 11]);
end
    
    %% Image the graph
    if any(isfinite(jelly.Nodes.x_coord)-1) == 1 || any(isfinite(jelly.Nodes.y_coord)-1) == 1
        return
    end
    
    if max(jelly.Edges.d_current) > 10
        return
    end
    
    figure1 = plot(jelly, 'XData', jelly.Nodes.x_coord, 'YData', jelly.Nodes.y_coord, 'EdgeCData', jelly.Edges.strain0, 'LineWidth', 0.8);

    cd(path1);                                                            % write the image data
    saveas(figure1, '000.jpg')
    cd(path0);     
    pause(0.001)

    c_count = 1;
    
    %% Start the sim
    for time = 1:time_steps
       
        %% Update the current length of edges
        for i = 1:numedges(jelly)
            [node1, node2] = findedge(jelly, i);
            dx = jelly.Nodes.x_coord(node1) - jelly.Nodes.x_coord(node2);
            dy = jelly.Nodes.y_coord(node1) - jelly.Nodes.y_coord(node2);
            dist_current = (dx^2 + dy^2)^(1/2);
            jelly.Edges.d_current(i) = dist_current;
        end
        
        muscle_length = sum(jelly.Edges.d_current(jelly.Edges.muscle == 1));
        mus_length = cat(1, mus_length, muscle_length);
        if max(jelly.Edges.d_current) > 10
            cd(path1);
            writematrix(mus_length, 'muscle_length.xlsx');
            cd(path0); 
        return
        end
        %% Calculate new muscle coordinates from contraction
        %muscles are "synchronized", calculations for all muscle bands happen
        %simultaneously.
        jelly.Nodes.F_muscle = contraction3(contraction_strength, jelly).*sin(c_count*pi/(contraction_duration/time_step));
        jelly = SLM_elastic(jelly, elast0, elast1);

        %% pressure force
        %jelly.Nodes.F_pressure = find_f_pressure(jelly, area_relax, bulk_modulus);  

        if c_count > 0
            %contracting
            jelly.Nodes.F_net = jelly.Nodes.F_elastic + jelly.Nodes.F_muscle - damping.*jelly.Nodes.velocity; %+ jelly.Nodes.F_pressure 
            c_count = c_count + 1;
            if c_count >= contraction_duration/time_step
                c_count = -1;
            end
        else
            %relaxing
            c_count = c_count - 1;
            jelly.Nodes.F_net = jelly.Nodes.F_elastic - damping.*jelly.Nodes.velocity; %+ jelly.Nodes.F_pressure 
            if c_count <= -1*(relax_duration/time_step)
                c_count = 1;
            end
        end

        %% Update the position of each node
        jelly.Nodes.velocity = jelly.Nodes.velocity + jelly.Nodes.F_net./mass*time_step;
        contraction_displacement = jelly.Nodes.velocity*time_step;
        jelly.Nodes.x_coord = jelly.Nodes.x_coord + contraction_displacement(:,1);
        jelly.Nodes.y_coord = jelly.Nodes.y_coord + contraction_displacement(:,2);

        figure1 = plot(jelly, 'XData', jelly.Nodes.x_coord, 'YData', jelly.Nodes.y_coord, 'EdgeCData', jelly.Edges.strain0 + jelly.Edges.strain1, 'LineWidth', 0.8);
%         hold on
%         figure1 = quiver(jelly.Nodes.x_coord, jelly.Nodes.y_coord, jelly.Nodes.F_net(:,1), jelly.Nodes.F_net(:,2));
%         hold off
        %% Save images
        cd(path1);                                                            % write the image data
        saveas(figure1, [num2str(time) '.jpg'])
        cd(path0);     
        pause(0.001)
            

    end
cd(path1);
writematrix(mus_length, 'muscle_length.xlsx');
cd(path0); 