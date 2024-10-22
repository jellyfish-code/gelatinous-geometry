%{
======================================================================
    Function that (?)
======================================================================

INPUT:
        jelly (?):
        contraction_strength (?): 
        muscle_strain (?): 
        max_dR (?): 
        dR_rate (?):

OUTPUT:
        jelly (?):
        done (?):
%}

function [jelly, done] = contraction5offset_old(jelly, contraction_strength, muscle_strain, max_dR, dR_rate)
%Reset F_muscle
F_muscle = zeros(numnodes(jelly), 2);
%Define the different piece of muscle
[trash, idx] = sortrows(jelly.Nodes, {'outmus'});
mus_idx = idx(trash.outmus~=0);
mus_idx = cat(1, mus_idx, mus_idx(1));
edge1 = 1;
mus_piece = [];
new_mus_piece = [];
mus_lengths = [];
mus_length_true = 0;
mus_anchor = [];

done = 1;
%% Need to find inflection point for butterfly
for i = 1:length(mus_idx)-1
    a = findedge(jelly, mus_idx(i), mus_idx(i+1));   
    new_mus_piece = cat(2, new_mus_piece, mus_idx(i));
    
    if a == 0 || jelly.Edges.muscle(a) == 0
        edge2 = i;
        check1 = outedges(jelly, mus_idx(edge1));
        check2 = outedges(jelly, mus_idx(edge2));
        check_2 = 0;
        anchor = 0;
        if length(check1) > 2
            anchor = anchor + 1;
        end
        if length(check2) > 2
            anchor = anchor + 1;
            check_2 = 1;
        end
        %make anchored end list first
        if anchor == 1 && check_2 == 1
            new_mus_piece = flip(new_mus_piece);
        end
        %newpiece = [mus_idx(edge1), mus_idx(edge2), mus_length_true];
        num = size(mus_piece);
        if num(1) > 0
            if num(2) > length(new_mus_piece)
                new_mus_piece = cat(2, new_mus_piece, zeros(1, num(2) - length(new_mus_piece)));
            elseif length(new_mus_piece) > num(2)
                mus_piece = cat(2, mus_piece, zeros(num(1), length(new_mus_piece)-num(2)));
            end
        end
        %update the table
        mus_piece = cat(1, mus_piece, new_mus_piece);
        mus_anchor = cat(1, mus_anchor, anchor);
        mus_lengths = cat(1, mus_lengths, mus_length_true);
        %initiate new muscle
        edge1 = i+1;
        new_mus_piece = [];
        mus_length_true = 0;
    else
        mus_length_true = mus_length_true + jelly.Edges.d_current(a);
    end
    
end

if isempty(mus_piece)
    mus_piece = cat(1, mus_piece, new_mus_piece);
    mus_anchor = cat(1, mus_anchor, 3);
    mus_lengths = cat(1, mus_lengths, mus_length_true);
end
    
%Each row of mus_piece now defines a muscle piece with edges and midpoint
%defined.
%For each piece, find the ends and the midpoint
m = size(mus_piece);
theta = zeros(m(1),1);
radius = zeros(m(1),1);
center = zeros(m(1), 2);
for i = 1:m(1)
    if mus_anchor(i) < 3
        a = find(mus_piece(i,:));
        cur_mus = mus_piece(i,a);
        mus_mid = ceil(length(cur_mus)/2);
        x1 = jelly.Nodes.x_coord(cur_mus(1));
        y1 = jelly.Nodes.y_coord(cur_mus(1));
        x2 = jelly.Nodes.x_coord(cur_mus(mus_mid));
        y2 = jelly.Nodes.y_coord(cur_mus(mus_mid));
        x3 = jelly.Nodes.x_coord(cur_mus(length(cur_mus)));
        y3 = jelly.Nodes.y_coord(cur_mus(length(cur_mus)));

    
        %Find the muscle length
        %x1, y1 = end1, x2 y2 = midpoint, x3 y3 = end2, Ox, Oy = center of circle
        syms ox oy
        eqn1 = (2*x1-2*x2)*ox + (2*y1 - 2*y2)*oy == x1^2 + y1^2 - x2^2 - y2^2;
        eqn2 = (2*x1-2*x3)*ox + (2*y1 - 2*y3)*oy == x1^2 + y1^2 - x3^2 - y3^2;

        [A, B] = equationsToMatrix([eqn1, eqn2]);
        X = linsolve(A, B);
        Ox = X(1); Oy = X(2); 


        %Can now define the radius
        radius(i) = ((x1 - Ox)^2 + (y1 - Oy)^2)^(1/2);
        center(i,:) = X;
    else
        a = find(mus_piece(i,:));
        cur_mus = mus_piece(i,a);
        center(i,1) = mean(jelly.Nodes.x_coord(cur_mus));
        center(i,2) = mean(jelly.Nodes.y_coord(cur_mus));
        radius(i) = mean(((jelly.Nodes.x_coord(cur_mus) - center(i,1)).^2 + (jelly.Nodes.y_coord(cur_mus) - center(i,2)).^2).^(1/2));
    end
    
    for j = 1:length(find(mus_piece(i,:)))
        X = jelly.Nodes.x_coord(mus_piece(i,j));
        Y = jelly.Nodes.y_coord(mus_piece(i,j));
        theta(i,j) = atan2(Y - center(i,2), X-center(i,1));
    end
end

%create table to store muscle info
mus_tab = table(mus_lengths, mus_anchor, mus_piece, radius, theta, center, 'VariableNames',["length", "anchored", "Nodes", "radius", "theta", "center"]);

%% Yay, all muscle info is extracted! Now to find the contracted valueeee
%new muscle contraction
%function F_muscle = ThreeDcontract(jelly, theta, radius, mus_idx)

%muscle_strain = 0.15;
contracted_nodes = zeros(height(mus_tab), length(mus_tab.Nodes(i,:)),2);

%%structure for each muscle piece as a table...?
%Find the new "contracted" position of muscle pieces
for i = 1:height(mus_tab)
    if mus_tab.anchored(i) == 1             
        
        %% which way is the muscle facing? pos theta or neg theta?
        edge = findedge(jelly, mus_tab.Nodes(i,1), mus_tab.Nodes(i,2));
        d_theta1 = jelly.Edges.d_current(edge)*muscle_strain/mus_tab.radius(i);
        if mus_tab.theta(i,2) < mus_tab.theta(i,1)
            if abs(mus_tab.theta(i,2) - mus_tab.theta(i,1)) < pi 
                d_theta1 = -1*d_theta1;
            end
        elseif mus_tab.theta(i,2) - mus_tab.theta(i,1) > pi
            d_theta1 = -1*d_theta1;
        end
        
        O = mus_tab.center(i,:);
        actual = find(mus_tab.Nodes(i,:));
        midpt = ceil(length(actual)/2);
        mid_theta = mus_tab.theta(i,midpt);
        cur_theta = mid_theta;
        
        %Midpoint goes in by dR
        new_R =mus_tab.radius(i) - min(max_dR, dR_rate*midpt);
        Ax = O(1) + new_R*(jelly.Nodes.x_coord(mus_tab.Nodes(i,midpt))-O(1))/mus_tab.radius(i);
        Ay = O(2) + new_R*(jelly.Nodes.y_coord(mus_tab.Nodes(i,midpt))-O(2))/mus_tab.radius(i);
        
        contracted_nodes(i,midpt,1) = Ax;
        contracted_nodes(i,midpt,2) = Ay;
        F_muscle(mus_tab.Nodes(i,midpt),1) = F_muscle(mus_tab.Nodes(i,midpt),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,midpt)));
        F_muscle(mus_tab.Nodes(i,midpt),2) = F_muscle(mus_tab.Nodes(i,midpt),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,midpt)));
        
        for j = midpt-1:-1:1
            edge = findedge(jelly, mus_tab.Nodes(i,j+1), mus_tab.Nodes(i,j));
            L = jelly.Edges.d_current(edge);
            %Define new mus, new R
            new_m =L*(1-muscle_strain);
            new_R =mus_tab.radius(i) - min(max_dR, dR_rate*j);
            d_theta = new_m/new_R;
            
            if d_theta1 > 0
                Ax = new_R*cos(cur_theta - abs(d_theta)) + O(1);
                Ay = new_R*sin(cur_theta - abs(d_theta)) + O(2); 
                
                contracted_nodes(i,j,1) = Ax;
                contracted_nodes(i,j,2) = Ay;
                F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));
                
                cur_theta = atan2(Ay - O(2), Ax - O(1));
                
            else
                Ax = new_R*cos(cur_theta + abs(d_theta)) + O(1);
                Ay = new_R*sin(cur_theta + abs(d_theta)) + O(2); 
                
                contracted_nodes(i,j,1) = Ax;
                contracted_nodes(i,j,2) = Ay;
                F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));
                
                cur_theta = atan2(Ay - O(2), Ax - O(1));
                
            end
            
        end
        cur_theta = mid_theta;
        
        for j = midpt+1:length(actual)
            edge = findedge(jelly, mus_tab.Nodes(i,j-1), mus_tab.Nodes(i,j));
            L = jelly.Edges.d_current(edge);
            %Define new mus, new R
            new_m =L*(1-muscle_strain);
            new_R =mus_tab.radius(i) - min(max_dR, dR_rate*j);
            d_theta = new_m/new_R;
            
            if d_theta1 > 0
                Ax = new_R*cos(cur_theta + abs(d_theta)) + O(1);
                Ay = new_R*sin(cur_theta + abs(d_theta)) + O(2); 
                
                contracted_nodes(i,j,1) = Ax;
                contracted_nodes(i,j,2) = Ay;
                F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));
                
                cur_theta = atan2(Ay - O(2), Ax - O(1));
                
            else
                Ax = new_R*cos(cur_theta - abs(d_theta)) + O(1);
                Ay = new_R*sin(cur_theta - abs(d_theta)) + O(2); 
                
                contracted_nodes(i,j,1) = Ax;
                contracted_nodes(i,j,2) = Ay;
                F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));
                
                cur_theta = atan2(Ay - O(2), Ax - O(1));
                
            end
        end
        
    elseif mus_tab.anchored(i) == 2

        O = mus_tab.center(i,:);
        actual = find(mus_tab.Nodes(i,:));

        %% first move anchors inward by mus_strain. first anchor:
        edge = findedge(jelly, mus_tab.Nodes(i,1), mus_tab.Nodes(i,2));
        d_theta1 = jelly.Edges.d_current(edge)*muscle_strain/mus_tab.radius(i);
        if mus_tab.theta(i,2) < mus_tab.theta(i,1)
            if abs(mus_tab.theta(i,2) - mus_tab.theta(i,1)) < pi 
                d_theta1 = -1*d_theta1;
            end
        elseif mus_tab.theta(i,2) - mus_tab.theta(i,1) > pi
            d_theta1 = -1*d_theta1;
        end

        midpt = ceil(length(actual)/2);
        mid_theta = mus_tab.theta(i,midpt);
        cur_theta = mid_theta;

        %Midpoint goes in by dR
        new_R =mus_tab.radius(i) - min(max_dR, dR_rate*midpt);
        Ax = O(1) + new_R*(jelly.Nodes.x_coord(mus_tab.Nodes(i,midpt))-O(1))/mus_tab.radius(i);
        Ay = O(2) + new_R*(jelly.Nodes.y_coord(mus_tab.Nodes(i,midpt))-O(2))/mus_tab.radius(i);

        contracted_nodes(i,midpt,1) = Ax;
        contracted_nodes(i,midpt,2) = Ay;
        F_muscle(mus_tab.Nodes(i,midpt),1) = F_muscle(mus_tab.Nodes(i,midpt),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,midpt)));
        F_muscle(mus_tab.Nodes(i,midpt),2) = F_muscle(mus_tab.Nodes(i,midpt),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,midpt)));

        for j = midpt-1:-1:1
            edge = findedge(jelly, mus_tab.Nodes(i,j+1), mus_tab.Nodes(i,j));
            L = jelly.Edges.d_current(edge);
            %Define new mus, new R
            new_m =L*(1-muscle_strain);
            new_R =mus_tab.radius(i) - min(max_dR, dR_rate*j);
            d_theta = new_m/new_R;

            if d_theta1 > 0
                Ax = new_R*cos(cur_theta - abs(d_theta)) + O(1);
                Ay = new_R*sin(cur_theta - abs(d_theta)) + O(2); 

                contracted_nodes(i,j,1) = Ax;
                contracted_nodes(i,j,2) = Ay;
                F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));

                cur_theta = atan2(Ay - O(2), Ax - O(1));
            else
                Ax = new_R*cos(cur_theta + abs(d_theta)) + O(1);
                Ay = new_R*sin(cur_theta + abs(d_theta)) + O(2); 

                contracted_nodes(i,j,1) = Ax;
                contracted_nodes(i,j,2) = Ay;
                F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));

                cur_theta = atan2(Ay - O(2), Ax - O(1));
            end

        end
        cur_theta = mid_theta;

        if length(actual) > 2
            for j = midpt+1:length(actual)
                edge = findedge(jelly, mus_tab.Nodes(i,j-1), mus_tab.Nodes(i,j));
                L = jelly.Edges.d_current(edge);
                %Define new mus, new R
                new_m =L*(1-muscle_strain);
                new_R =mus_tab.radius(i) - min(max_dR, dR_rate*(length(actual)-j));
                d_theta = new_m/new_R;

                if d_theta1 > 0
                    Ax = new_R*cos(cur_theta + abs(d_theta)) + O(1);
                    Ay = new_R*sin(cur_theta + abs(d_theta)) + O(2); 

                    contracted_nodes(i,j,1) = Ax;
                    contracted_nodes(i,j,2) = Ay;
                    F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                    F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));

                    cur_theta = atan2(Ay - O(2), Ax - O(1));

                else
                    Ax = new_R*cos(cur_theta - abs(d_theta)) + O(1);
                    Ay = new_R*sin(cur_theta - abs(d_theta)) + O(2); 

                    contracted_nodes(i,j,1) = Ax;
                    contracted_nodes(i,j,2) = Ay;
                    F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                    F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));

                    cur_theta = atan2(Ay - O(2), Ax - O(1));

                end
            end
        end
    elseif mus_tab.anchored(i) == 0
        %start in the middle
        %muscle_strain = 0.15;
        O = mus_tab.center(i,:);
        actual = find(mus_tab.Nodes(i,:));
        midpt = ceil(length(actual)/2);
        new_R = mus_tab.radius(i) - max_dR;
        Ax = O(1) + new_R*(jelly.Nodes.x_coord(mus_tab.Nodes(i,midpt))-O(1))/mus_tab.radius(i);
        Ay = O(2) + new_R*(jelly.Nodes.y_coord(mus_tab.Nodes(i,midpt))-O(2))/mus_tab.radius(i);
        
        contracted_nodes(i,midpt,1) = Ax;
        contracted_nodes(i,midpt,2) = Ay;
        
        
        for j = midpt:length(actual)-1
            edge = findedge(jelly, mus_tab.Nodes(i,j), mus_tab.Nodes(i,j+1));
            L = jelly.Edges.d_current(edge);
            
            %Define new mus, new R
            new_m =L*(1-muscle_strain);

            %define the center of circle Ox, Oy
            
            C1 = new_m^2 - new_R^2 + O(1)^2 + O(2)^2 - Ax^2 - Ay^2;
            C2 = O(2) - Ay;
            C3 = O(1) - Ax;
            
            p = [C2^2/C3^2+1, -1*C2*C1/C3^2+2*Ax*C2/C3-2*Ay, (C1/(2*C3))^2-Ax*C1/C3+Ax^2+Ay^2-new_m^2];
            if any(isinf(p)==1) || any(isnan(p)==1)
                done = 0;
            else
                y1 = roots(p);
                x1 = (C1 - 2*y1*C2)/(2*C3);
            end
            
            v1 = [x1(1)-Ax, y1(1)-Ay];
            v2 = [x1(2)-Ax, y1(2)-Ay];
            v0 = [jelly.Nodes.x_coord(mus_tab.Nodes(i,j+1)) - jelly.Nodes.x_coord(mus_tab.Nodes(i,j)), jelly.Nodes.y_coord(mus_tab.Nodes(i,j+1))-jelly.Nodes.y_coord(mus_tab.Nodes(i,j))];
            d1 = sqrt(sum((v0-v1).^2));
            d2 = sqrt(sum((v0-v2).^2));
            if d1 < d2
                Ax = x1(1);
                Ay = y1(1);
            else
                Ax = x1(2);
                Ay = y1(2);
            end
            
            F_muscle(mus_tab.Nodes(i,j+1),1) = F_muscle(mus_tab.Nodes(i,j+1),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j+1)));
            F_muscle(mus_tab.Nodes(i,j+1),2) = F_muscle(mus_tab.Nodes(i,j+1),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j+1)));
            %save the points
            contracted_nodes(i,j+1,1) = Ax;
            contracted_nodes(i,j+1,2) = Ay;
        end
        
        Ax = contracted_nodes(i,midpt,1);
        Ay = contracted_nodes(i,midpt,2);
        
        if length(actual) > 2
            for j = midpt:-1:2
                edge = findedge(jelly, mus_tab.Nodes(i,j), mus_tab.Nodes(i,j-1));
                L = jelly.Edges.d_current(edge);

                %Define new mus, new R
                new_m =L*(1-muscle_strain);

                %define the center of circle Ox, Oy

                C1 = new_m^2 - new_R^2 + O(1)^2 + O(2)^2 - Ax^2 - Ay^2;
                C2 = O(2) - Ay;
                C3 = O(1) - Ax;

                p = [C2^2/C3^2+1, -1*C2*C1/C3^2+2*Ax*C2/C3-2*Ay, (C1/(2*C3))^2-Ax*C1/C3+Ax^2+Ay^2-new_m^2];
                if any(isinf(p)==1) || any(isnan(p)==1)
                    done = 0;
                else
                    y1 = roots(p);
                    x1 = (C1 - 2*y1*C2)/(2*C3);
                end

                v1 = [x1(1)-Ax, y1(1)-Ay];
                v2 = [x1(2)-Ax, y1(2)-Ay];
                v0 = [jelly.Nodes.x_coord(mus_tab.Nodes(i,j-1)) - jelly.Nodes.x_coord(mus_tab.Nodes(i,j)), jelly.Nodes.y_coord(mus_tab.Nodes(i,j-1))-jelly.Nodes.y_coord(mus_tab.Nodes(i,j))];
                d1 = sqrt(sum((v0-v1).^2));
                d2 = sqrt(sum((v0-v2).^2));
                if d1 < d2
                    Ax = x1(1);
                    Ay = y1(1);
                else
                    Ax = x1(2);
                    Ay = y1(2);
                end

                F_muscle(mus_tab.Nodes(i,j-1),1) = F_muscle(mus_tab.Nodes(i,j-1),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j-1)));
                F_muscle(mus_tab.Nodes(i,j-1),2) = F_muscle(mus_tab.Nodes(i,j-1),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j-1)));
                %save the points
                contracted_nodes(i,j-1,1) = Ax;
                contracted_nodes(i,j-1,2) = Ay;
            end
        end
          
    elseif mus_tab.anchored(i) == 3
        O = mus_tab.center(i,:);
        actual = find(mus_tab.Nodes(i,:));
        % everything moves toward the center by a set amount
        new_R = mus_tab.radius(i) - max_dR;
        
        for j = 1:length(actual)
            Ax = O(1) + new_R*(jelly.Nodes.x_coord(mus_tab.Nodes(i,j))-O(1))/mus_tab.radius(i);
            Ay = O(2) + new_R*(jelly.Nodes.y_coord(mus_tab.Nodes(i,j))-O(2))/mus_tab.radius(i);
            
            F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
            F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));
        end
    end
end

%% inner muscle

%Define the different piece of muscle
[trash, idx] = sortrows(jelly.Nodes, {'inmus'});
mus_idx = idx(trash.inmus~=0);
mus_idx = cat(1, mus_idx, mus_idx(1));
edge1 = 1;
mus_piece = [];
new_mus_piece = [];
mus_lengths = [];
mus_length_true = 0;
mus_anchor = [];

%% Need to find inflection point for butterfly
for i = 1:length(mus_idx)-1
    a = findedge(jelly, mus_idx(i), mus_idx(i+1));   
    new_mus_piece = cat(2, new_mus_piece, mus_idx(i));
    
    % do we need to start a new muscle?
    %yes if i+1 is not connected to i, then i is in the current muscle but
    %a new muscle should be initiated   
    if a == 0 || jelly.Edges.muscle(a) == 0
        edge2 = i;
        check1 = outedges(jelly, mus_idx(edge1));
        check2 = outedges(jelly, mus_idx(edge2));
        check_2 = 0;
        anchor = 0;
        if length(check1) > 4
            anchor = anchor + 1;
        end
        if length(check2) > 4
            anchor = anchor + 1;
            check_2 = 1;
        end
        %make anchored end list first
        if anchor == 1 && check_2 == 1
            new_mus_piece = flip(new_mus_piece);
        end
        %newpiece = [mus_idx(edge1), mus_idx(edge2), mus_length_true];
        num = size(mus_piece);
        if num(1) > 0
            if num(2) > length(new_mus_piece)
                new_mus_piece = cat(2, new_mus_piece, zeros(1, num(2) - length(new_mus_piece)));
            elseif length(new_mus_piece) > num(2)
                mus_piece = cat(2, mus_piece, zeros(num(1), length(new_mus_piece)-num(2)));
            end
        end
        %update the table
        mus_piece = cat(1, mus_piece, new_mus_piece);
        mus_anchor = cat(1, mus_anchor, anchor);
        mus_lengths = cat(1, mus_lengths, mus_length_true);
        %initiate new muscle
        edge1 = i+1;
        new_mus_piece = [];
        mus_length_true = 0;
    else
        mus_length_true = mus_length_true + jelly.Edges.d_current(a);
    end
    
end
%Each row of mus_piece now defines a muscle piece with edges and midpoint
%defined.
%For each piece, find the ends and the midpoint
m = size(mus_piece);
theta = zeros(m(1),1);
radius = zeros(m(1),1);
center = zeros(m(1), 2);
for i = 1:m(1)
    if mus_anchor(i) < 3
        a = find(mus_piece(i,:));
        cur_mus = mus_piece(i,a);
        mus_mid = ceil(length(cur_mus)/2);
        x1 = jelly.Nodes.x_coord(cur_mus(1));
        y1 = jelly.Nodes.y_coord(cur_mus(1));
        x2 = jelly.Nodes.x_coord(cur_mus(mus_mid));
        y2 = jelly.Nodes.y_coord(cur_mus(mus_mid));
        x3 = jelly.Nodes.x_coord(cur_mus(length(cur_mus)));
        y3 = jelly.Nodes.y_coord(cur_mus(length(cur_mus)));

    
        %Find the muscle length
        %x1, y1 = end1, x2 y2 = midpoint, x3 y3 = end2, Ox, Oy = center of circle
        syms ox oy
        eqn1 = (2*x1-2*x2)*ox + (2*y1 - 2*y2)*oy == x1^2 + y1^2 - x2^2 - y2^2;
        eqn2 = (2*x1-2*x3)*ox + (2*y1 - 2*y3)*oy == x1^2 + y1^2 - x3^2 - y3^2;

        [A, B] = equationsToMatrix([eqn1, eqn2]);
        X = linsolve(A, B);
        Ox = X(1); Oy = X(2); 


        %Can now define the radius
        radius(i) = ((x1 - Ox)^2 + (y1 - Oy)^2)^(1/2);
        center(i,:) = X;
    else
        a = find(mus_piece(i,:));
        cur_mus = mus_piece(i,a);
        center(i,1) = mean(jelly.Nodes.x_coord(cur_mus));
        center(i,2) = mean(jelly.Nodes.y_coord(cur_mus));
        radius(i) = mean(((jelly.Nodes.x_coord(cur_mus) - center(i,1)).^2 + (jelly.Nodes.y_coord(cur_mus) - center(i,2)).^2).^(1/2));
    end
    
    for j = 1:length(find(mus_piece(i,:)))
        X = jelly.Nodes.x_coord(mus_piece(i,j));
        Y = jelly.Nodes.y_coord(mus_piece(i,j));
        theta(i,j) = atan2(Y - center(i,2), X-center(i,1));
    end
end

%create table to store muscle info
mus_tab = table(mus_lengths, mus_anchor, mus_piece, radius, theta, center, 'VariableNames',["length", "anchored", "Nodes", "radius", "theta", "center"]);

%% Yay, all muscle info is extracted! Now to find the contracted valueeee
%new muscle contraction
%function F_muscle = ThreeDcontract(jelly, theta, radius, mus_idx)

%muscle_strain = 0.15;
contracted_nodes = zeros(height(mus_tab), length(mus_tab.Nodes(i,:)),2);

%%structure for each muscle piece as a table...?
%Find the new "contracted" position of muscle pieces
for i = 1:height(mus_tab)
    if mus_tab.anchored(i) == 1
                %% which way is the muscle facing? pos theta or neg theta?
        edge = findedge(jelly, mus_tab.Nodes(i,1), mus_tab.Nodes(i,2));
        d_theta1 = jelly.Edges.d_current(edge)*muscle_strain/mus_tab.radius(i);
        if mus_tab.theta(i,2) < mus_tab.theta(i,1)
            if abs(mus_tab.theta(i,2) - mus_tab.theta(i,1)) < pi 
                d_theta1 = -1*d_theta1;
            end
        elseif mus_tab.theta(i,2) - mus_tab.theta(i,1) > pi
            d_theta1 = -1*d_theta1;
        end
        
        %prioritize dR
        %priority = 'dR';
        %muscle_strain = 0.15;
        O = mus_tab.center(i,:);
        actual = find(mus_tab.Nodes(i,:));
        midpt = ceil(length(actual)/2);
        mid_theta = mus_tab.theta(i,midpt);
        cur_theta = mid_theta;
        
        %Midpoint goes in by dR
        new_R =mus_tab.radius(i) - min(max_dR, dR_rate*midpt);
        Ax = O(1) + new_R*(jelly.Nodes.x_coord(mus_tab.Nodes(i,midpt))-O(1))/mus_tab.radius(i);
        Ay = O(2) + new_R*(jelly.Nodes.y_coord(mus_tab.Nodes(i,midpt))-O(2))/mus_tab.radius(i);
        
        contracted_nodes(i,midpt,1) = Ax;
        contracted_nodes(i,midpt,2) = Ay;
        F_muscle(mus_tab.Nodes(i,midpt),1) = F_muscle(mus_tab.Nodes(i,midpt),1) + 0.5*contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,midpt)));
        F_muscle(mus_tab.Nodes(i,midpt),2) = F_muscle(mus_tab.Nodes(i,midpt),2) + 0.5*contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,midpt)));
        
        for j = midpt-1:-1:1
            edge = findedge(jelly, mus_tab.Nodes(i,j+1), mus_tab.Nodes(i,j));
            L = jelly.Edges.d_current(edge);
            %Define new mus, new R
            new_m =L*(1-muscle_strain);
            new_R =mus_tab.radius(i) - min(max_dR, dR_rate*j);
            d_theta = new_m/new_R;
            
            if d_theta1 > 0
                Ax = new_R*cos(cur_theta - abs(d_theta)) + O(1);
                Ay = new_R*sin(cur_theta - abs(d_theta)) + O(2); 
                
                contracted_nodes(i,j,1) = Ax;
                contracted_nodes(i,j,2) = Ay;
                F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + 0.5*contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + 0.5*contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));
                
                cur_theta = atan2(Ay - O(2), Ax - O(1));
                
            else
                Ax = new_R*cos(cur_theta + abs(d_theta)) + O(1);
                Ay = new_R*sin(cur_theta + abs(d_theta)) + O(2); 
                
                contracted_nodes(i,j,1) = Ax;
                contracted_nodes(i,j,2) = Ay;
                F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + 0.5*contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + 0.5*contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));
                
                cur_theta = atan2(Ay - O(2), Ax - O(1));
                
            end
            
        end
        cur_theta = mid_theta;
        
        for j = midpt+1:length(actual)
            edge = findedge(jelly, mus_tab.Nodes(i,j-1), mus_tab.Nodes(i,j));
            L = jelly.Edges.d_current(edge);
            %Define new mus, new R
            new_m =L*(1-muscle_strain);
            new_R =mus_tab.radius(i) - min(max_dR, dR_rate*j);
            d_theta = new_m/new_R;
            
            if d_theta1 > 0
                Ax = new_R*cos(cur_theta + abs(d_theta)) + O(1);
                Ay = new_R*sin(cur_theta + abs(d_theta)) + O(2); 
                
                contracted_nodes(i,j,1) = Ax;
                contracted_nodes(i,j,2) = Ay;
                F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + 0.5*contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + 0.5*contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));
                
                cur_theta = atan2(Ay - O(2), Ax - O(1));
                
            else
                Ax = new_R*cos(cur_theta - abs(d_theta)) + O(1);
                Ay = new_R*sin(cur_theta - abs(d_theta)) + O(2); 
                
                contracted_nodes(i,j,1) = Ax;
                contracted_nodes(i,j,2) = Ay;
                F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + 0.5*contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + 0.5*contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));
                
                cur_theta = atan2(Ay - O(2), Ax - O(1));
                
            end
        end

    elseif mus_tab.anchored(i) == 2

        %priority = 'dR';
        %muscle_strain = 0.15;
        O = mus_tab.center(i,:);
        actual = find(mus_tab.Nodes(i,:));

        %% first move anchors inward by mus_strain. first anchor:
        edge = findedge(jelly, mus_tab.Nodes(i,1), mus_tab.Nodes(i,2));
        d_theta1 = jelly.Edges.d_current(edge)*muscle_strain/mus_tab.radius(i);
        if mus_tab.theta(i,2) < mus_tab.theta(i,1)
            if abs(mus_tab.theta(i,2) - mus_tab.theta(i,1)) < pi 
                d_theta1 = -1*d_theta1;
            end
        elseif mus_tab.theta(i,2) - mus_tab.theta(i,1) > pi
            d_theta1 = -1*d_theta1;
        end

        midpt = ceil(length(actual)/2);
        mid_theta = mus_tab.theta(i,midpt);
        cur_theta = mid_theta;

        %Midpoint goes in by dR
        new_R =mus_tab.radius(i) - min(max_dR, dR_rate*midpt);
        Ax = O(1) + new_R*(jelly.Nodes.x_coord(mus_tab.Nodes(i,midpt))-O(1))/mus_tab.radius(i);
        Ay = O(2) + new_R*(jelly.Nodes.y_coord(mus_tab.Nodes(i,midpt))-O(2))/mus_tab.radius(i);

        contracted_nodes(i,midpt,1) = Ax;
        contracted_nodes(i,midpt,2) = Ay;
        F_muscle(mus_tab.Nodes(i,midpt),1) = F_muscle(mus_tab.Nodes(i,midpt),1) + 0.5*contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,midpt)));
        F_muscle(mus_tab.Nodes(i,midpt),2) = F_muscle(mus_tab.Nodes(i,midpt),2) + 0.5*contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,midpt)));

        for j = midpt-1:-1:1
            edge = findedge(jelly, mus_tab.Nodes(i,j+1), mus_tab.Nodes(i,j));
            L = jelly.Edges.d_current(edge);
            %Define new mus, new R
            new_m =L*(1-muscle_strain);
            new_R =mus_tab.radius(i) - min(max_dR, dR_rate*j);
            d_theta = new_m/new_R;

            if d_theta1 > 0
                Ax = new_R*cos(cur_theta - abs(d_theta)) + O(1);
                Ay = new_R*sin(cur_theta - abs(d_theta)) + O(2); 

                contracted_nodes(i,j,1) = Ax;
                contracted_nodes(i,j,2) = Ay;
                F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + 0.5*contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + 0.5*contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));

                cur_theta = atan2(Ay - O(2), Ax - O(1));
            else
                Ax = new_R*cos(cur_theta + abs(d_theta)) + O(1);
                Ay = new_R*sin(cur_theta + abs(d_theta)) + O(2); 

                contracted_nodes(i,j,1) = Ax;
                contracted_nodes(i,j,2) = Ay;
                F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + 0.5*contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + 0.5*contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));

                cur_theta = atan2(Ay - O(2), Ax - O(1));
            end

        end
        cur_theta = mid_theta;

        if length(actual) > 2
            for j = midpt+1:length(actual)
                edge = findedge(jelly, mus_tab.Nodes(i,j-1), mus_tab.Nodes(i,j));
                L = jelly.Edges.d_current(edge);
                %Define new mus, new R
                new_m =L*(1-muscle_strain);
                new_R =mus_tab.radius(i) - min(max_dR, dR_rate*(length(actual)-j));
                d_theta = new_m/new_R;

                if d_theta1 > 0
                    Ax = new_R*cos(cur_theta + abs(d_theta)) + O(1);
                    Ay = new_R*sin(cur_theta + abs(d_theta)) + O(2); 

                    contracted_nodes(i,j,1) = Ax;
                    contracted_nodes(i,j,2) = Ay;
                    F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + 0.5*contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                    F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + 0.5*contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));

                    cur_theta = atan2(Ay - O(2), Ax - O(1));

                else
                    Ax = new_R*cos(cur_theta - abs(d_theta)) + O(1);
                    Ay = new_R*sin(cur_theta - abs(d_theta)) + O(2); 

                    contracted_nodes(i,j,1) = Ax;
                    contracted_nodes(i,j,2) = Ay;
                    F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + 0.5*contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
                    F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + 0.5*contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));

                    cur_theta = atan2(Ay - O(2), Ax - O(1));

                end
            end
        end
    elseif mus_tab.anchored(i) == 0
        %start in the middle
        %muscle_strain = 0.15;
        O = mus_tab.center(i,:);
        actual = find(mus_tab.Nodes(i,:));
        midpt = ceil(length(actual)/2);
        new_R = mus_tab.radius(i) - max_dR;
        Ax = O(1) + new_R*(jelly.Nodes.x_coord(mus_tab.Nodes(i,midpt))-O(1))/mus_tab.radius(i);
        Ay = O(2) + new_R*(jelly.Nodes.y_coord(mus_tab.Nodes(i,midpt))-O(2))/mus_tab.radius(i);
        
        contracted_nodes(i,midpt,1) = Ax;
        contracted_nodes(i,midpt,2) = Ay;
        
        
        for j = midpt:length(actual)-1
            edge = findedge(jelly, mus_tab.Nodes(i,j), mus_tab.Nodes(i,j+1));
            L = jelly.Edges.d_current(edge);
            
            %Define new mus, new R
            new_m =L*(1-muscle_strain);

            %define the center of circle Ox, Oy
            
            C1 = new_m^2 - new_R^2 + O(1)^2 + O(2)^2 - Ax^2 - Ay^2;
            C2 = O(2) - Ay;
            C3 = O(1) - Ax;
            
            p = [C2^2/C3^2+1, -1*C2*C1/C3^2+2*Ax*C2/C3-2*Ay, (C1/(2*C3))^2-Ax*C1/C3+Ax^2+Ay^2-new_m^2];
            if any(isinf(p)==1) || any(isnan(p)==1)
                done = 0;
            else
                y1 = roots(p);
                x1 = (C1 - 2*y1*C2)/(2*C3);
            end
            
            v1 = [x1(1)-Ax, y1(1)-Ay];
            v2 = [x1(2)-Ax, y1(2)-Ay];
            v0 = [jelly.Nodes.x_coord(mus_tab.Nodes(i,j+1)) - jelly.Nodes.x_coord(mus_tab.Nodes(i,j)), jelly.Nodes.y_coord(mus_tab.Nodes(i,j+1))-jelly.Nodes.y_coord(mus_tab.Nodes(i,j))];
            d1 = sqrt(sum((v0-v1).^2));
            d2 = sqrt(sum((v0-v2).^2));
            if d1 < d2
                Ax = x1(1);
                Ay = y1(1);
            else
                Ax = x1(2);
                Ay = y1(2);
            end
            
            F_muscle(mus_tab.Nodes(i,j+1),1) = F_muscle(mus_tab.Nodes(i,j+1),1) + 0.5*contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j+1)));
            F_muscle(mus_tab.Nodes(i,j+1),2) = F_muscle(mus_tab.Nodes(i,j+1),2) + 0.5*contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j+1)));
            %save the points
            contracted_nodes(i,j+1,1) = Ax;
            contracted_nodes(i,j+1,2) = Ay;
        end
        
        Ax = contracted_nodes(i,midpt,1);
        Ay = contracted_nodes(i,midpt,2);
        
        if length(actual) > 2
            for j = midpt:-1:2
                edge = findedge(jelly, mus_tab.Nodes(i,j), mus_tab.Nodes(i,j-1));
                L = jelly.Edges.d_current(edge);

                %Define new mus, new R
                new_m =L*(1-muscle_strain);

                %define the center of circle Ox, Oy

                C1 = new_m^2 - new_R^2 + O(1)^2 + O(2)^2 - Ax^2 - Ay^2;
                C2 = O(2) - Ay;
                C3 = O(1) - Ax;

                p = [C2^2/C3^2+1, -1*C2*C1/C3^2+2*Ax*C2/C3-2*Ay, (C1/(2*C3))^2-Ax*C1/C3+Ax^2+Ay^2-new_m^2];
                if any(isinf(p)==1) || any(isnan(p)==1)
                    done = 0;
                else
                    y1 = roots(p);
                    x1 = (C1 - 2*y1*C2)/(2*C3);
                end

                v1 = [x1(1)-Ax, y1(1)-Ay];
                v2 = [x1(2)-Ax, y1(2)-Ay];
                v0 = [jelly.Nodes.x_coord(mus_tab.Nodes(i,j-1)) - jelly.Nodes.x_coord(mus_tab.Nodes(i,j)), jelly.Nodes.y_coord(mus_tab.Nodes(i,j-1))-jelly.Nodes.y_coord(mus_tab.Nodes(i,j))];
                d1 = sqrt(sum((v0-v1).^2));
                d2 = sqrt(sum((v0-v2).^2));
                if d1 < d2
                    Ax = x1(1);
                    Ay = y1(1);
                else
                    Ax = x1(2);
                    Ay = y1(2);
                end

                F_muscle(mus_tab.Nodes(i,j-1),1) = F_muscle(mus_tab.Nodes(i,j-1),1) + 0.5*contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j-1)));
                F_muscle(mus_tab.Nodes(i,j-1),2) = F_muscle(mus_tab.Nodes(i,j-1),2) + 0.5*contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j-1)));
                %save the points
                contracted_nodes(i,j-1,1) = Ax;
                contracted_nodes(i,j-1,2) = Ay;
            end
        end
          
    elseif mus_tab.anchored(i) == 3
        O = mus_tab.center(i,:);
        actual = find(mus_tab.Nodes(i,:));
        % everything moves toward the center by a set amount
        new_R = mus_tab.radius(i) - max_dR;
        
        for j = 1:length(actual)
            Ax = O(1) + new_R*(jelly.Nodes.x_coord(mus_tab.Nodes(i,j))-O(1))/mus_tab.radius(i);
            Ay = O(2) + new_R*(jelly.Nodes.y_coord(mus_tab.Nodes(i,j))-O(2))/mus_tab.radius(i);
            
            F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
            F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));
        end
        
    end
end
jelly.Nodes.stress_muscle = F_muscle;
