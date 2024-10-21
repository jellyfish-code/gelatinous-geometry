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

function [jelly, done] = contraction4offset(jelly, contraction_strength, muscle_strain, max_dR, dR_rate)
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
% prev_slope = 1;
% prev_inf = 0;
% prev_infcos = 0;
%prev_tan = 0;
% prev_cot = 0;
done = 1;
%% Need to find inflection point for butterfly
for i = 1:length(mus_idx)-1
    a = findedge(jelly, mus_idx(i), mus_idx(i+1));   
    %slope = (jelly.Nodes.y_coord(mus_idx(i+1)) - jelly.Nodes.y_coord(mus_idx(i)))/(jelly.Nodes.x_coord(mus_idx(i+1)) - jelly.Nodes.x_coord(mus_idx(i)));
%     inf = (slope - prev_slope)/abs(slope-prev_slope);
%     inf_cos = (acot(1/slope) - acot(1/prev_slope))/abs(acot(1/slope) - acot(1/prev_slope));
    %cur_tan = atan(slope);
%     cur_cot = acot(1/slope);
    %node i gets added to the current muscle 
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
%         prev_inf = 0;
        %prev_tan =0;
        
        
    %yes if inflection changed, then i is in the current muscle and also in
    %the next muscle
%     elseif prev_inf~=0 && inf~=prev_inf && inf_cos~=prev_infcos
%     elseif prev_tan~=0 && abs(abs(cur_tan-prev_tan)-pi/2) < pi/2 - 0.4 %&& abs(cur_cot-prev_cot) > 0.3
%         edge2 = i;
%         check1 = outedges(jelly, mus_idx(edge1));
%         check2 = outedges(jelly, mus_idx(edge2));
%         check_2 = 0;
%         anchor = 0;
%         if length(check1) > 2
%             anchor = anchor + 1;
%         end
%         if length(check2) > 2
%             anchor = anchor + 1;
%             check_2 = 1;
%         end
%         %make anchored end list first
%         if anchor == 1 && check_2 == 1
%             new_mus_piece = flip(new_mus_piece);
%         end
%         %newpiece = [mus_idx(edge1), mus_idx(edge2), mus_length_true];
%         num = size(mus_piece);
%         if num(1) > 0
%             if num(2) > length(new_mus_piece)
%                 new_mus_piece = cat(2, new_mus_piece, zeros(1,num(2) - length(new_mus_piece)));
%             elseif length(new_mus_piece) > num(2)
%                 mus_piece = cat(2, mus_piece, zeros(num(1), length(new_mus_piece)-num(2)));
%             end
%         end
%         %update the table
%         mus_piece = cat(1, mus_piece, new_mus_piece);
%         mus_anchor = cat(1, mus_anchor, anchor);
%         mus_lengths = cat(1, mus_lengths, mus_length_true);
%         %initiate new muscle
%         edge1 = i;
%         new_mus_piece = [];
%         mus_length_true = 0;
% %         prev_inf = 0;
%         new_mus_piece = cat(2, new_mus_piece, mus_idx(i));
%         mus_length_true = mus_length_true + jelly.Edges.d_current(a);
%         prev_tan=0;
    else
%         if length(new_mus_piece) > 1
%             prev_inf = inf;
%             prev_infcos = inf_cos;
%         end
%         prev_slope = slope;
%         prev_cot = cur_cot;
        %prev_tan = cur_tan;
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
    for j = 1:length(find(mus_piece(i,:)))
        X = jelly.Nodes.x_coord(mus_piece(i,j));
        Y = jelly.Nodes.y_coord(mus_piece(i,j));
        theta(i,j) = atan2(Y - Oy, X-Ox);
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
        O = mus_tab.center(i,:);
        actual = find(mus_tab.Nodes(i,:));
        
        %% rotate by d_current*musclestrain first
        % what theta is an arclength = d_current*musstrain about the
        % center? shift all nodes that theta 
        edge = findedge(jelly, mus_tab.Nodes(i,1), mus_tab.Nodes(i,2));
        d_theta = jelly.Edges.d_current(edge)*muscle_strain/mus_tab.radius(i);
        if mus_tab.theta(i,2) < mus_tab.theta(i,1)
            if abs(mus_tab.theta(i,2) - mus_tab.theta(i,1)) < pi 
                d_theta = -1*d_theta;
            end
        elseif mus_tab.theta(i,2) - mus_tab.theta(i,1) > pi
            d_theta = -1*d_theta;
        end
        

        Ax = cos(mus_tab.theta(i,1)+d_theta)*mus_tab.radius(i) + O(1);
        Ay = sin(mus_tab.theta(i,1)+d_theta)*mus_tab.radius(i) + O(2);

        contracted_nodes(i,1,1) = Ax;
        contracted_nodes(i,1,2) = Ay;
     
        F_muscle(mus_tab.Nodes(i,1),1) = F_muscle(mus_tab.Nodes(i,1),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,1)));
        F_muscle(mus_tab.Nodes(i,1),2) = F_muscle(mus_tab.Nodes(i,1),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,1)));
        
        for j = 1:length(actual)-1

            edge = findedge(jelly, mus_tab.Nodes(i,j), mus_tab.Nodes(i,j+1));
            L = jelly.Edges.d_current(edge);
            
            %Define new mus, new R
            new_R =mus_tab.radius(i) - min(max_dR, dR_rate*j);
            new_m =L*(1-muscle_strain);

            %define the center of circle Ox, Oy
            
            C1 = new_m^2 - new_R^2 + O(1)^2 + O(2)^2 - Ax^2 - Ay^2;
            C2 = O(2) - Ay;
            C3 = O(1) - Ax;
            
            p = [C2^2/C3^2+1, -1*C2*C1/C3^2+2*Ax*C2/C3-2*Ay, (C1/(2*C3))^2-Ax*C1/C3+Ax^2+Ay^2-new_m^2];
            
            if any(isinf(p)==1) || any(isnan(p)==1)
                done = 0;
                return
            else
                y1 = roots(p);
                x1 = (C1 - 2*y1*C2)/(2*C3);
            end
            
            %% Check!
%             check = [new_m, sqrt((x1(1)-Ax)^2 + (y1(1)-Ay)^2)]
%             check = [new_m, sqrt((x1(2)-Ax)^2 + (y1(2)-Ay)^2)]
%             check = [new_R, sqrt((x1(1)-O(1))^2 + (y1(1)-O(2))^2)]
%             check = [new_R, sqrt((x1(2)-O(1))^2 + (y1(2)-O(2))^2)]
            %test which point is the correct one with theta
            if isreal(y1) && isreal(x1)
                t1 = mod(atan2(y1(1) - O(2), x1(1) - O(1)),2*pi)-mod(atan2(contracted_nodes(i,j,2) - O(2), contracted_nodes(i,j,1)-O(1)), 2*pi);
                if abs(t1) > pi
                    t1 = mod(atan2(y1(1) - O(2), x1(1) - O(1)),2*pi) - atan2(contracted_nodes(i,j,2) - O(2), contracted_nodes(i,j,1)-O(1));
                end
                t2 = mod(atan2(y1(2) - O(2), x1(2) - O(1)),2*pi)-mod(atan2(contracted_nodes(i,j,2) - O(2), contracted_nodes(i,j,1)-O(1)), 2*pi);
                if abs(t2) > pi
                    t2 = mod(atan2(y1(2) - O(2), x1(2) - O(1)),2*pi) - atan2(contracted_nodes(i,j,2) - O(2), contracted_nodes(i,j,1)-O(1));
                end
            else
                done=0;
                return
            end
            if d_theta > 0
                if t1 > 0
                    Ax = x1(1); Ay = y1(1);
                elseif t2 > 0
                    Ax = x1(2); Ay = y1(2);
                end
            else
                if t1 < 0
                    Ax = x1(1); Ay = y1(1);
                elseif t2 < 0
                    Ax = x1(2); Ay = y1(2);
                end
            end
            
       
%             v1 = [x1(1)-Ax, y1(1)-Ay];
%             v2 = [x1(2)-Ax, y1(2)-Ay];
%             v0 = [mus_temp(1,j+1) - mus_temp(1,j), mus_temp(2,j+1)-mus_temp(2,j)];
%             d1 = sqrt(sum((v0-v1).^2));
%             d2 = sqrt(sum((v0-v2).^2));
%             if d1 < d2
%                 Ax = x1(1);
%                 Ay = y1(1);
%             else
%                 Ax = x1(2);
%                 Ay = y1(2);
%             end
            
            F_muscle(mus_tab.Nodes(i,j+1),1) = F_muscle(mus_tab.Nodes(i,j+1),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j+1)));
            F_muscle(mus_tab.Nodes(i,j+1),2) = F_muscle(mus_tab.Nodes(i,j+1),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j+1)));
            %save the points
            contracted_nodes(i,j+1,1) = Ax;
            contracted_nodes(i,j+1,2) = Ay;
            
            
            %save points as Ax, Ay
        end


    elseif mus_tab.anchored(i) == 2
        %prioritize dR
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
        Ax = cos(mus_tab.theta(i,1)+d_theta1)*mus_tab.radius(i) + O(1);
        Ay = sin(mus_tab.theta(i,1)+d_theta1)*mus_tab.radius(i) + O(2);
        contracted_nodes(i,1,1) = Ax;
        contracted_nodes(i,1,2) = Ay;
        F_muscle(mus_tab.Nodes(i,1),1) = F_muscle(mus_tab.Nodes(i,1),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,1)));
        F_muscle(mus_tab.Nodes(i,1),2) = F_muscle(mus_tab.Nodes(i,1),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,1)));
        anc1_t = atan2(Ay - O(2), Ax-O(1));
        
        %second anchor
        edge = findedge(jelly, mus_tab.Nodes(i,length(actual)), mus_tab.Nodes(i,length(actual)-1));
        d_theta = jelly.Edges.d_current(edge)*muscle_strain/mus_tab.radius(i);
        if mus_tab.theta(i,length(actual)-1) < mus_tab.theta(i,length(actual))
            if abs(mus_tab.theta(i,length(actual)-1) - mus_tab.theta(i,length(actual))) < pi 
                d_theta = -1*d_theta;
            end
        elseif mus_tab.theta(i,length(actual)-1) - mus_tab.theta(i,length(actual)) > pi
            d_theta = -1*d_theta;
        end        
        Ax = cos(mus_tab.theta(i,length(actual))+d_theta)*mus_tab.radius(i) + O(1);
        Ay = sin(mus_tab.theta(i,length(actual))+d_theta)*mus_tab.radius(i) + O(2);
        contracted_nodes(i,length(actual),1) = Ax;
        contracted_nodes(i,length(actual),2) = Ay;
        F_muscle(mus_tab.Nodes(i,length(actual)),1) = F_muscle(mus_tab.Nodes(i,length(actual)),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,length(actual))));
        F_muscle(mus_tab.Nodes(i,length(actual)),2) = F_muscle(mus_tab.Nodes(i,length(actual)),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,length(actual))));
        anc2_t = atan2(Ay - O(2), Ax - O(1));
        
        %%now for nodes between the anchors
        d_theta = (mod(anc2_t,2*pi) - mod(anc1_t,2*pi))/(length(actual)-1);
        if d_theta1 < 0
            d_theta = -1*abs(d_theta);
        end
        for j = 2:length(actual)-1
            %each node moves toward the center by dR*distance from anchors
            d = min(length(actual)-j, j-1);

            Ax = (mus_tab.radius(i) - d*dR_rate)*cos(anc1_t + d_theta*(j-1)) + O(1);
            Ay = (mus_tab.radius(i) - d*dR_rate)*sin(anc1_t + d_theta*(j-1)) + O(2);
%         %vector is center - current
%         vx = (O(1) - jelly.Nodes.x_coord(mus_tab.Nodes(i,j)))/mus_tab.radius(i);
%         vy = (O(2) - jelly.Nodes.y_coord(mus_tab.Nodes(i,j)))/mus_tab.radius(i);
% 
%         Ax = jelly.Nodes.x_coord(mus_tab.Nodes(i,j)) + d*dR_rate*vx;
%         Ay = jelly.Nodes.y_coord(mus_tab.Nodes(i,j)) + d*dR_rate*vy;

        F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
        F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));
        %save the points
        contracted_nodes(i,j,1) = Ax;
        contracted_nodes(i,j,2) = Ay;
        end
        
        %check strain
%         for j = 1:length(actual)-1
%             edge = findedge(jelly, mus_tab.Nodes(i,j), mus_tab.Nodes(i,j+1));
%             L = jelly.Edges.d_current(edge);
%             new_m =L*(1-muscle_strain);
% 
%             dx = contracted_nodes(i,j+1,1) - contracted_nodes(i,j,1);
%             dy = contracted_nodes(i,j+1,2) - contracted_nodes(i,j,2);
% 
% %             if sqrt((dx^2 + dy^2)) < new_m
% %                 %priority = 'strain';
% %             end
%         end
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
                return
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
                    return
                else
                    y1 = roots(p);
                    x1 = (C1 - 2*y1*C2)/(2*C3);
                end

                if isreal(y1) && isreal(x1)
                
                    v1 = [x1(1)-Ax, y1(1)-Ay];
                    v2 = [x1(2)-Ax, y1(2)-Ay];
                else
                    done=0;
                    return
                end
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
% prev_slope = 1;
% prev_inf = 0;
% prev_infcos = 0;
%prev_tan = 0;
% prev_cot = 0;
%% Need to find inflection point for butterfly
for i = 1:length(mus_idx)-1
    a = findedge(jelly, mus_idx(i), mus_idx(i+1));   
    %slope = (jelly.Nodes.y_coord(mus_idx(i+1)) - jelly.Nodes.y_coord(mus_idx(i)))/(jelly.Nodes.x_coord(mus_idx(i+1)) - jelly.Nodes.x_coord(mus_idx(i)));
%     inf = (slope - prev_slope)/abs(slope-prev_slope);
%     inf_cos = (acot(1/slope) - acot(1/prev_slope))/abs(acot(1/slope) - acot(1/prev_slope));
    %cur_tan = atan(slope);
%     cur_cot = acot(1/slope);
    %node i gets added to the current muscle 
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
%         prev_inf = 0;
        %prev_tan =0;
        
        
    %yes if inflection changed, then i is in the current muscle and also in
    %the next muscle
%     elseif prev_inf~=0 && inf~=prev_inf && inf_cos~=prev_infcos
%     elseif prev_tan~=0 && abs(abs(cur_tan-prev_tan)-pi/2) < pi/2 - 0.4 %&& abs(cur_cot-prev_cot) > 0.3
%         edge2 = i;
%         check1 = outedges(jelly, mus_idx(edge1));
%         check2 = outedges(jelly, mus_idx(edge2));
%         check_2 = 0;
%         anchor = 0;
%         if length(check1) > 4
%             anchor = anchor + 1;
%         end
%         if length(check2) > 4
%             anchor = anchor + 1;
%             check_2 = 1;
%         end
%         %make anchored end list first
%         if anchor == 1 && check_2 == 1
%             new_mus_piece = flip(new_mus_piece);
%         end
%         %newpiece = [mus_idx(edge1), mus_idx(edge2), mus_length_true];
%         num = size(mus_piece);
%         if num(1) > 0
%             if num(2) > length(new_mus_piece)
%                 new_mus_piece = cat(2, new_mus_piece, zeros(1,num(2) - length(new_mus_piece)));
%             elseif length(new_mus_piece) > num(2)
%                 mus_piece = cat(2, mus_piece, zeros(num(1), length(new_mus_piece)-num(2)));
%             end
%         end
%         %update the table
%         mus_piece = cat(1, mus_piece, new_mus_piece);
%         mus_anchor = cat(1, mus_anchor, anchor);
%         mus_lengths = cat(1, mus_lengths, mus_length_true);
%         %initiate new muscle
%         edge1 = i;
%         new_mus_piece = [];
%         mus_length_true = 0;
% %         prev_inf = 0;
%         new_mus_piece = cat(2, new_mus_piece, mus_idx(i));
%         mus_length_true = mus_length_true + jelly.Edges.d_current(a);
%         prev_tan=0;
    else
%         if length(new_mus_piece) > 1
%             prev_inf = inf;
%             prev_infcos = inf_cos;
%         end
%         prev_slope = slope;
%         prev_cot = cur_cot;
        %prev_tan = cur_tan;
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
    for j = 1:length(find(mus_piece(i,:)))
        X = jelly.Nodes.x_coord(mus_piece(i,j));
        Y = jelly.Nodes.y_coord(mus_piece(i,j));
        theta(i,j) = atan2(Y - Oy, X-Ox);
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
        O = mus_tab.center(i,:);
        actual = find(mus_tab.Nodes(i,:));
        
        %% rotate by d_current*musclestrain first
        % what theta is an arclength = d_current*musstrain about the
        % center? shift all nodes that theta 
        edge = findedge(jelly, mus_tab.Nodes(i,1), mus_tab.Nodes(i,2));
        d_theta = jelly.Edges.d_current(edge)*muscle_strain/mus_tab.radius(i);
        if mus_tab.theta(i,2) < mus_tab.theta(i,1)
            if abs(mus_tab.theta(i,2) - mus_tab.theta(i,1)) < pi 
                d_theta = -1*d_theta;
            end
        elseif mus_tab.theta(i,2) - mus_tab.theta(i,1) > pi
            d_theta = -1*d_theta;
        end
        

        Ax = cos(mus_tab.theta(i,1)+d_theta)*mus_tab.radius(i) + O(1);
        Ay = sin(mus_tab.theta(i,1)+d_theta)*mus_tab.radius(i) + O(2);

        contracted_nodes(i,1,1) = Ax;
        contracted_nodes(i,1,2) = Ay;
     
        F_muscle(mus_tab.Nodes(i,1),1) = F_muscle(mus_tab.Nodes(i,1),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,1)));
        F_muscle(mus_tab.Nodes(i,1),2) = F_muscle(mus_tab.Nodes(i,1),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,1)));
        
        for j = 1:length(actual)-1

            edge = findedge(jelly, mus_tab.Nodes(i,j), mus_tab.Nodes(i,j+1));
            L = jelly.Edges.d_current(edge);
            
            %Define new mus, new R
            new_R =mus_tab.radius(i) - min(max_dR, dR_rate*j);
            new_m =L*(1-muscle_strain);

            %define the center of circle Ox, Oy
            
            C1 = new_m^2 - new_R^2 + O(1)^2 + O(2)^2 - Ax^2 - Ay^2;
            C2 = O(2) - Ay;
            C3 = O(1) - Ax;
            
            p = [C2^2/C3^2+1, -1*C2*C1/C3^2+2*Ax*C2/C3-2*Ay, (C1/(2*C3))^2-Ax*C1/C3+Ax^2+Ay^2-new_m^2];
            
            if any(isinf(p)==1) || any(isnan(p)==1)
                done = 0;
                return
            else
                y1 = roots(p);
                x1 = (C1 - 2*y1*C2)/(2*C3);
            end
            
            if isreal(y1) && isreal(x1)
            
            %% Check!
%             check = [new_m, sqrt((x1(1)-Ax)^2 + (y1(1)-Ay)^2)]
%             check = [new_m, sqrt((x1(2)-Ax)^2 + (y1(2)-Ay)^2)]
%             check = [new_R, sqrt((x1(1)-O(1))^2 + (y1(1)-O(2))^2)]
%             check = [new_R, sqrt((x1(2)-O(1))^2 + (y1(2)-O(2))^2)]
            %test which point is the correct one with theta
                t1 = mod(atan2(y1(1) - O(2), x1(1) - O(1)),2*pi)-mod(atan2(contracted_nodes(i,j,2) - O(2), contracted_nodes(i,j,1)-O(1)), 2*pi);
                if abs(t1) > pi
                    t1 = mod(atan2(y1(1) - O(2), x1(1) - O(1)),2*pi) - atan2(contracted_nodes(i,j,2) - O(2), contracted_nodes(i,j,1)-O(1));
                end
                t2 = mod(atan2(y1(2) - O(2), x1(2) - O(1)),2*pi)-mod(atan2(contracted_nodes(i,j,2) - O(2), contracted_nodes(i,j,1)-O(1)), 2*pi);
                if abs(t2) > pi
                    t2 = mod(atan2(y1(2) - O(2), x1(2) - O(1)),2*pi) - atan2(contracted_nodes(i,j,2) - O(2), contracted_nodes(i,j,1)-O(1));
                end
            else
                done = 0;
                return
            end
            if d_theta > 0
                if t1 > 0
                    Ax = x1(1); Ay = y1(1);
                elseif t2 > 0
                    Ax = x1(2); Ay = y1(2);
                end
            else
                if t1 < 0
                    Ax = x1(1); Ay = y1(1);
                elseif t2 < 0
                    Ax = x1(2); Ay = y1(2);
                end
            end
            
       
%             v1 = [x1(1)-Ax, y1(1)-Ay];
%             v2 = [x1(2)-Ax, y1(2)-Ay];
%             v0 = [mus_temp(1,j+1) - mus_temp(1,j), mus_temp(2,j+1)-mus_temp(2,j)];
%             d1 = sqrt(sum((v0-v1).^2));
%             d2 = sqrt(sum((v0-v2).^2));
%             if d1 < d2
%                 Ax = x1(1);
%                 Ay = y1(1);
%             else
%                 Ax = x1(2);
%                 Ay = y1(2);
%             end
            
            F_muscle(mus_tab.Nodes(i,j+1),1) = F_muscle(mus_tab.Nodes(i,j+1),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j+1)));
            F_muscle(mus_tab.Nodes(i,j+1),2) = F_muscle(mus_tab.Nodes(i,j+1),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j+1)));
            %save the points
            contracted_nodes(i,j+1,1) = Ax;
            contracted_nodes(i,j+1,2) = Ay;
            
            
            %save points as Ax, Ay
        end


    elseif mus_tab.anchored(i) == 2
        %prioritize dR
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
        Ax = cos(mus_tab.theta(i,1)+d_theta1)*mus_tab.radius(i) + O(1);
        Ay = sin(mus_tab.theta(i,1)+d_theta1)*mus_tab.radius(i) + O(2);
        contracted_nodes(i,1,1) = Ax;
        contracted_nodes(i,1,2) = Ay;
        F_muscle(mus_tab.Nodes(i,1),1) = F_muscle(mus_tab.Nodes(i,1),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,1)));
        F_muscle(mus_tab.Nodes(i,1),2) = F_muscle(mus_tab.Nodes(i,1),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,1)));
        anc1_t = atan2(Ay - O(2), Ax-O(1));
        
        %second anchor
        edge = findedge(jelly, mus_tab.Nodes(i,length(actual)), mus_tab.Nodes(i,length(actual)-1));
        d_theta = jelly.Edges.d_current(edge)*muscle_strain/mus_tab.radius(i);
        if mus_tab.theta(i,length(actual)-1) < mus_tab.theta(i,length(actual))
            if abs(mus_tab.theta(i,length(actual)-1) - mus_tab.theta(i,length(actual))) < pi 
                d_theta = -1*d_theta;
            end
        elseif mus_tab.theta(i,length(actual)-1) - mus_tab.theta(i,length(actual)) > pi
            d_theta = -1*d_theta;
        end        
        Ax = cos(mus_tab.theta(i,length(actual))+d_theta)*mus_tab.radius(i) + O(1);
        Ay = sin(mus_tab.theta(i,length(actual))+d_theta)*mus_tab.radius(i) + O(2);
        contracted_nodes(i,length(actual),1) = Ax;
        contracted_nodes(i,length(actual),2) = Ay;
        F_muscle(mus_tab.Nodes(i,length(actual)),1) = F_muscle(mus_tab.Nodes(i,length(actual)),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,length(actual))));
        F_muscle(mus_tab.Nodes(i,length(actual)),2) = F_muscle(mus_tab.Nodes(i,length(actual)),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,length(actual))));
        anc2_t = atan2(Ay - O(2), Ax - O(1));
        
        %%now for nodes between the anchors
        d_theta = (mod(anc2_t,2*pi) - mod(anc1_t,2*pi))/(length(actual)-1);
        if d_theta1 < 0
            d_theta = -1*abs(d_theta);
        end
        for j = 2:length(actual)-1
            %each node moves toward the center by dR*distance from anchors
            d = min(length(actual)-j, j-1);

            Ax = (mus_tab.radius(i) - d*dR_rate)*cos(anc1_t + d_theta*(j-1)) + O(1);
            Ay = (mus_tab.radius(i) - d*dR_rate)*sin(anc1_t + d_theta*(j-1)) + O(2);
%         %vector is center - current
%         vx = (O(1) - jelly.Nodes.x_coord(mus_tab.Nodes(i,j)))/mus_tab.radius(i);
%         vy = (O(2) - jelly.Nodes.y_coord(mus_tab.Nodes(i,j)))/mus_tab.radius(i);
% 
%         Ax = jelly.Nodes.x_coord(mus_tab.Nodes(i,j)) + d*dR_rate*vx;
%         Ay = jelly.Nodes.y_coord(mus_tab.Nodes(i,j)) + d*dR_rate*vy;

        F_muscle(mus_tab.Nodes(i,j),1) = F_muscle(mus_tab.Nodes(i,j),1) + contraction_strength*(Ax-jelly.Nodes.x_coord(mus_tab.Nodes(i,j)));
        F_muscle(mus_tab.Nodes(i,j),2) = F_muscle(mus_tab.Nodes(i,j),2) + contraction_strength*(Ay-jelly.Nodes.y_coord(mus_tab.Nodes(i,j)));
        %save the points
        contracted_nodes(i,j,1) = Ax;
        contracted_nodes(i,j,2) = Ay;
        end
        
        %check strain
%         for j = 1:length(actual)-1
%             edge = findedge(jelly, mus_tab.Nodes(i,j), mus_tab.Nodes(i,j+1));
%             L = jelly.Edges.d_current(edge);
%             new_m =L*(1-muscle_strain);
% 
%             dx = contracted_nodes(i,j+1,1) - contracted_nodes(i,j,1);
%             dy = contracted_nodes(i,j+1,2) - contracted_nodes(i,j,2);
% 
% %             if sqrt((dx^2 + dy^2)) < new_m
% %                 %priority = 'strain';
% %             end
%         end
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
                return
            else
                y1 = roots(p);
                x1 = (C1 - 2*y1*C2)/(2*C3);
            end
            
            if isreal(y1) && isreal(x1)
                v1 = [x1(1)-Ax, y1(1)-Ay];
                v2 = [x1(2)-Ax, y1(2)-Ay];     
            else
                done=0;
                return
            end
            
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
                    return
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
          
        
    end
end
jelly.Nodes.F_muscle = F_muscle;
