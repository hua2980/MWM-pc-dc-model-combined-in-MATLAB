% 211209 WH modified from Charlotte(WANG Shiqi) and Wenyue's R codes
% chose pc-based or dc-based model by weight: always choose dc-based model
% in the end 

% takes weights, PC coordinates & parameters and returns new weights and performance measures (single trial)
function [weights_pc,weights_dc, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants,...
    latency, speed_std, speed_ps, mean_angle,time_step] = ...
    run_trial (curr_model, weights0_pc,weights0_dc, Wmult, sigma_pc, sigma_dc,sigma_ac, PC_x, PC_y, DC,Vdecay, ac_const, beta, etdecay, ...
    lrate, discf, noise, platform_x, platform_y, starting_x, starting_y, speed, wall_pun)
% speed in meters per second (!) taken as input - latency computed accordingly
% *Add parameters for DC
    
    % FIXED PARAMETRES OF THE EXPERMENT

    pool_diameter = 1.4; %Maze diameter (\phi) in metres (m)
    platform_radius = 0.06; %Platform radius

    N_dc = 10; % *Population of distance cells
    N_pc = 211; %Population of place cells
    N_ac = 36; %Population of action cells
    which_pc = 0; % *add which_pc
    which_dc = 0;
    
    count = 0; %initialize count of loop
    
    dist = 0;
    wall_zone = 0;
    quadrants = zeros(4,1); %Percentage spent on each quadrant 

    weights_pc = weights0_pc; %Initialize modifiable weights 
    weights_dc = weights0_dc; % *add initialize dc weights    


    el_tr_pc = zeros(N_pc, N_ac); %Initialize eligibility traces matrix
    el_tr_dc = zeros(N_dc, N_ac); % *add dc eligibility traces matrix

    %Initialize trajectories
    track_x = starting_x; %Current position of trajectory is equal to the 
                          %starting location of the animal
    track_y = starting_y;
    vel_x = 0;
    vel_y = 0;

    % NAVIGATION LOOP
    while ((track_x(end) - platform_x)^2 + (track_y(end) - platform_y)^2 > ... 
        platform_radius^2)

        % *add dc
        weights_pc = weights_pc*(1-noise) + rand(N_pc, N_ac)*Wmult*noise;
        weights_dc = weights_dc*(1-noise) + rand(N_dc, N_ac)*Wmult*noise;

        dist_to_wall = pool_diameter/2 - sqrt(track_x(end)^2+track_y(end)^2);

        %Calculate PC activation 
        % *add DC activation
        PC_activation = zeros(1,N_pc);
        for i = 1:N_pc
            PC_activation(i) = exp(-((track_x(end) - PC_x(i))^2 + ...
                (track_y(end)- PC_y(i))^2)/(2*sigma_pc^2));
        end
        DC_activation = zeros(1,N_dc);
        for i = 1:N_dc
            DC_activation(i) = exp(-(dist_to_wall - DC(i))^2/(2*sigma_dc^2));
        end

        %Calculate AC activation (i.e. value of the action, Q)
        if (length(track_x) > 1)
            prevQ_pc = AC_activation_pc(which_pc); %Displays the Q value before movement
            prevQ_dc = AC_activation_dc(which_dc); % *Add dc 
        end

        AC_activation_pc = PC_activation * weights_pc;
        AC_activation_dc = DC_activation * weights_dc;

        % *add dc
        %Make an action- softmax equation
        ACsel_pc = AC_activation_pc.^beta;
        ACsel_pc = max(0,ACsel_pc)./sum(ACsel_pc); 
        ACsel_dc= AC_activation_dc.^beta;
        ACsel_dc= max(0,ACsel_dc)./sum(ACsel_dc);

        % *add  #which = sample(c(0:35), 1, replace = TRUE, prob = ACsel)
        ASrand_pc = rand;
        which_pc = 1; 
        ASsum_pc = ACsel_pc(1);

        while (which_pc < N_ac && ASsum_pc < ASrand_pc)
            which_pc = which_pc + 1;
            ASsum_pc = ASsum_pc + ACsel_pc(which_pc);
        end

        % *add dc
        ASrand_dc = rand;
        which_dc = 1; 
        ASsum_dc = ACsel_dc(1);

        while (which_dc < N_ac && ASsum_dc < ASrand_dc)
            which_dc = which_dc + 1;
            ASsum_dc = ASsum_dc + ACsel_dc(which_dc);
        end

        %Eligibility traces
        el_tr_pc = el_tr_pc * etdecay;
        el_tr_dc = el_tr_dc * etdecay;

        for j = 1:N_ac
            itmp_pc = min(abs(j-which_pc),N_ac-abs(j-which_pc));
            itmp_dc = min(abs(j-which_dc),N_ac-abs(j-which_dc));
            actgaus_pc = exp(-(itmp_pc*itmp_pc)/(2*sigma_ac*sigma_ac));
            actgaus_dc = exp(-(itmp_dc*itmp_dc)/(2*sigma_ac*sigma_ac));
            el_tr_pc(:,j) = el_tr_pc(:,j) + actgaus_pc*AC_activation_pc(j)*(PC_activation');
            el_tr_dc(:,j) = el_tr_dc(:,j) + actgaus_dc*AC_activation_dc(j)*(DC_activation'); % *r_i^pc should be column matrix
        end 

        % *add direction
        if (track_x(end)>=0)
            central.angle = pi + atan(track_y(end)/track_x(end));
        else
            central.angle = atan(track_y(end)/track_x(end));
        end
        moving.dir = central.angle + which_dc/N_ac*2*pi;

        % *add dc
        if (curr_model ==1)
            vel_x = [vel_x (vel_x(end)+ac_const*cos(which_pc/N_ac*2*pi))*Vdecay];
            vel_y = [vel_y (vel_y(end)+ac_const*sin(which_pc/N_ac*2*pi))*Vdecay];
        else
            vel_x = [vel_x, (vel_x(end)+ac_const*cos(moving.dir))*Vdecay];
            vel_y = [vel_y, (vel_y(end)+ac_const*sin(moving.dir))*Vdecay];
        end

        %velocity per time step (not second)
        track_x = [track_x track_x(end)+vel_x(end)];
        track_y = [track_y track_y(end)+vel_y(end)];

        %Check if not out of bounds, reset location & speed if so
        if (track_x(end)^2 + track_y(end)^2 > (pool_diameter/2)^2)
            ratio = (track_x(end)^2 + track_y(end)^2)/((pool_diameter/2)^2);
            track_x(end) = track_x(end)/sqrt(ratio);
            track_y(end) = track_y(end)/sqrt(ratio);
            vel_x(end) = track_x(end) - track_x(end-1);
            vel_y(end) = track_y(end) - track_y(end-1);
        end

        if (length(track_x) > 2)
           if ((track_x(end) - platform_x)^2 + (track_y(end) - platform_y)^2 < ... 
        platform_radius^2)
               rew = 10; %found platform - reward
           elseif (track_x(end)^2+track_y(end)^2 > (0.99*pool_diameter/2)^2)        %%%%%%%%%%%%%why 0.99 maybe try to set the wall reward constant at 0? or -10?
               rew = -wall_pun; %hit wall - punishment
           else
               rew = 0; %didn't find - no reward
           end
           
           % *add DC
           currQ_pc = AC_activation_pc(which_pc);
           currQ_dc = AC_activation_dc(which_dc);
           %disp(['Qs: ',num2str(currQ),' - ',num2str(prevQ)]);
           tderr_pc = rew + discf*currQ_pc - prevQ_pc; %temporal difference error
           tderr_dc = rew + discf*currQ_dc - prevQ_dc; %temporal difference error
           weights_pc = max(0,weights_pc + lrate*tderr_pc*el_tr_pc);
           weights_dc = max(0,weights_dc + lrate*tderr_dc*el_tr_dc);
    %         for i = 1:N_ac
    %             for j = 1:N_pc
    %                 weights(j,i) = max(weights(j,i) + lrate*tderr*el_tr(j,i),0);
    %             end
    %         end
        end
        %disp([num2str(track_x(end)),',',num2str(track_y(end))]); %display the 
        %coordinates of each action cell within the loop
        laststep = sqrt((track_x(end)-track_x(end-1))^2 + (track_y(end)-track_y...
            (end-1))^2);
        dist = dist + laststep;


        if (track_x(end)^2 + track_y(end)^2 > 0.8*(pool_diameter/2)^2)
            wall_zone = wall_zone + 1;
        elseif (track_x(end) > 0 && track_y(end) > 0)
            quadrants(1) = quadrants(1) + 1;
        elseif (track_x(end) > 0 && track_y(end) < 0)
            quadrants(2) = quadrants(2) + 1;
        elseif (track_x(end) < 0 && track_y(end) < 0)
            quadrants(3) = quadrants(3) + 1;
        else
            quadrants(4) = quadrants(4) + 1;
        end

        if (length(track_x) > 100) % evaluate latency after 100+ steps
            speed_ts = mean((vel_x(2:end).^2+vel_y(2:end).^2).^0.5); % speed in meters/time step
            latency = (length(track_x)-1) * speed_ts / speed; % convert to seconds
            if (latency > 60) % if more than a minute, stop
                break;
            end
        end
        count = count + 1;
    end
    
    



    latency = length(track_x)-1;  % latency in time steps
    wall_zone = wall_zone/latency;
    quadrants = quadrants/latency;
    speed_ts = mean((vel_x(2:end).^2+vel_y(2:end).^2).^0.5); % speed in meters/time step

    % speed per action step from Hanbing
    speed_ps = (vel_x(2:end).^2+vel_y(2:end).^2).^0.5; 

    %time step
    time_step = speed_ts / speed; 

    %mean turning angle from Barbara
    vel = [vel_x' vel_y'];
    angle = [];
    for steps = 2:(length(vel_x)-1)
        A= vel(steps,:);
        B= vel(steps+1,:);
        angle = [angle acos(dot(A,B)/(norm(A)*norm(B)))]; % radian result
    end
    angle = angle * 180 / pi;
    % angle = min(angle, pi-angle) * speed_ts / speed;
    mean_angle = mean(angle);

    % speed standard deviation from Hua
    speed_std = std((vel_x(2:end).^2+vel_y(2:end).^2).^0.5,1); 
    speed_std = speed_std / time_step;

    latency = latency * time_step; % latency in seconds
end