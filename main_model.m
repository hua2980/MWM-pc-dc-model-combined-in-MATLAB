
% 211209 Wanghua modified from Charlotte(WANG Shiqi) and Wenyue's R code
% *add for combined model

%%
% Set parameters
tic

N_dc = 10; %Population of distance sensory cell [?]
N_pc = 211; %Population of place cells [100..300]
N_ac = 36; %Population of action cells [25..50]

plot_trajectories = 1; %yes - 1, no - 0
plot_cognitive_maps =  1; %yes - 1, no - 0
plot_dc_maps = 1; % yes-1, no-0
pln = plot_trajectories + plot_cognitive_maps + plot_dc_maps;
% if pln > 0.5
%     f=figure;
% end
Nruns = 1; % how many runs to run if not plotting anything

pool_diameter = 1.4; % #Maze diameter (\phi) in metres (m)
platform_radius = 0.06; % #Platform radius (m)
sigma_dc = 0.1; %place cell sigma (standard deviation), in meters [0.05..0.2]
sigma_pc = 0.1; %place cell sigma (standard deviation), in meters [0.05..0.2]
sigma_ac = 2;  % #action cell sigma (standard deviation), in action cells [1..3]

etdecay = 0.83; %Eligibility trace decay (lambda) [0.75..0.95] LESS THAN GAMMA!
alpha = 0.01; %Learning rate (\alpha) [0.005..0.02]
beta = 7;  %Exploration-exploitation factor (\beta) [0.5..12]
gamma = 0.85; %Discount factor (\gamma) [0.75..0.95]

Vdecay = 0.82; %velocity decay [0.75..0.95]
ac_const = 0.02; %acceleration const [0.01..0.03]
Wnoise =  0.0004; %Weight noise [0.0001..0.0007]
Wmult = 0.1; %Weight multiplier [0.05..0.15]
hitwall =  0.2; %punishment for hitting the wall [0..1]
speed =  0.175; %mouse speed (m/s) [0.1..0.25]

Ntrials = 10; %#number of trials per day
Ndays = 5; %#number of days
track_x_sum=cell(1, Ntrials*Ndays);
track_y_sum=cell(1, Ntrials*Ndays);

Npsets = 1;

cumm_weights_mod = zeros(2, Ntrials*Ndays); % for record the weight changes in the two model, hua 20211217
cumm_which_mod = zeros(1, Ntrials*Ndays); % for record the model changes, hua 20211217

% #performance measures to compute: latency, distance, time in target quadrant, opposite quadrant, and wall zone
% AM: performance measures to comput: speed per step
if (pln > 0.5) %if any plots
    clf
    PMs = zeros(8,Ndays,Ntrials);
    AMs = zeros(Ndays, Ntrials);
else
    PMs = zeros(8,Ndays,Ntrials,Nruns,Npsets);
    AMs = zeros(Ndays,Ntrials,Nruns); % #multiple runs
end

%for the real model
%Platform coordinates:
platform_x = cos(-pi/4)*pool_diameter/4; %x coordinate
platform_y = sin(-pi/4)*pool_diameter/4; %y coordinate
  
%Starting locations of the modeled animal (4 different ones):
strad = pool_diameter/2*0.85; %15% of maze radius to the wall
starting_xs = strad * [cos(pi/6) cos(pi/3) cos(7*pi/6) cos(4*pi/3)]; %x coordinates
starting_ys = strad * [sin(pi/6) sin(pi/3) sin(7*pi/6) sin(4*pi/3)]; %y coordinates
th = 0:pi/50:2*pi; %for plotting circles :)

if (pln > 0.5)
    %Generate initial weights
    weights_pc = rand(N_pc,N_ac)*Wmult;
    weights_dc = rand(N_dc,N_ac)*Wmult;
    weights_mod= rand(2,1)*Wmult;
    
    %Generate distance cells
    DC = zeros(1,N_dc);

    %Generate place cells
    PC_x = zeros(1,N_pc); %1xN_pc matrix containing the x = 0 coordinate for each place cell
    PC_y = zeros(1,N_pc); %1xN_pc matrix containing the y = 0 coordinate for each place cell

    for i = 1:N_dc %#For each distance cell:
        DC(i) = rand*(pool_diameter); %#Random positions of distance cells
    end
    
    for i = 1:N_pc %For each place cell:
        PC_x(i) = (rand - 0.5)*pool_diameter; %Random positions of place cells
        PC_y(i) = (rand - 0.5)*pool_diameter;
        while (PC_x(i)^2 + PC_y(i)^2 > (pool_diameter/2)^2) %Checks for out of bounds
            PC_x(i) = (rand - 0.5)*pool_diameter;
            PC_y(i) = (rand - 0.5)*pool_diameter;
        end
    end

%%
record = 0; % for test, created by hua 20211217 can be deleted later

    for day = 1 : Ndays
        for trial = 1:Ntrials
                        
            %# fix or variable platform(not match betaval code)
            %whichplatform = 1; %# fix
            whichplatform=1-randi([0,1],1)*2; %# variable
            
            idx = randi(4); %randomly choose one of 4 starting locations
            starting_x = starting_xs(idx);
            starting_y = starting_ys(idx);
            
            % add noise to weights_mod
            weights_mod = weights_mod*(1-Wnoise) + rand(2, 1)*Wmult*Wnoise;
            [curr_model, mod_prob] = model_selection (weights_mod, beta); % model selection
            prev_modQ = weights_mod; % output prev_modQ
            trial
            record = record + 1;
            cumm_which_mod(:,record) = curr_model; % record the choice
            
          
            %#run trial
            [wres_pc,wres_dc, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants,...
                latency, speed_std, speed_ps, mean_angle,time_step]=...
                run_trial(curr_model, weights_pc, weights_dc, Wmult, sigma_pc, ...
                sigma_dc, sigma_ac, PC_x, PC_y, DC, Vdecay, ac_const, beta, etdecay, ...
                alpha, gamma, Wnoise, whichplatform*platform_x, whichplatform*platform_y, ...
                starting_x, starting_y, speed, hitwall);   
            
            weights_pc = wres_pc;
            weights_dc = wres_dc;
            
            % compute reward:
            id = (1.5 - 0.5 * curr_model); %convert -1/1 into rows 1/2; pc-based (1) is row 1, dc-based (-1) is row 2
            rewd0=zeros(2,1);
            rew = 0;
            if ((track_x(end) - whichplatform*platform_x)^2 + (track_y(end) - whichplatform*platform_y)^2 < ... 
               platform_radius^2)
               rew = rew+10; %found platform - reward
            end
            if wall_zone > 0
               rew = rew-hitwall; %hit wall - punishment
            end
            rewd0(id)= rewd0(id) + rew; % get reward after decision (say, rewd = 1)    
            % update weights_mod
            [wres_mod, ~] = update_weight (rewd0, weights_mod, gamma, prev_modQ, alpha);
            wres_mod
    
            weights_mod = wres_mod;
            cumm_weights_mod(:,record) = wres_mod;
            
      
            PMs(1,day,trial) = latency; 
            PMs(2,day,trial) = dist; 
            if whichplatform == 1
                PMs(3,day,trial) = quadrants(4,1)*100; %#target quadrant
                PMs(4,day,trial) = quadrants(2,1)*100; %#opposite quadrant
            else
                PMs(3,day,trial) = quadrants(2,1)*100; %#target quadrant
                PMs(4,day,trial) = quadrants(4,1)*100; %#opposite quadrant
            end
            PMs(5,day,trial) = wall_zone*100; %#wall zone
            PMs(6,day,trial) = speed_std*100; %# speed_std
            PMs(7,day,trial) = mean_angle;
            PMs(8,day,trial) = time_step;
            
            if (day == Ndays && trial == Ntrials)
            if (plot_trajectories)
            %subplot(Ndays, Ntrials*pln,(day-1)*Ntrials*pln+(trial-1)*pln+1); 
            subplot(1, pln, 1); 

            hold on
            %plot the trajectory
            for i = 1:(length(track_x))-1
                line(track_x(i:i+1),track_y(i:i+1),'Color',[i/length(track_x),0,1-i/length(track_x)]);
            end
            %plot the maze and platform
            plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),'k');
            plot(whichplatform*platform_x+platform_radius*cos(th),whichplatform*platform_y+platform_radius*sin(th),'k');
            end

            if (plot_dc_maps)
            %subplot(Ndays, Ntrials*pln,(day-1)*Ntrials*pln+(trial-1)*pln + 2);
            subplot(1, pln, 2);

            hold on
            %plot the cognitive map (for DC)
            for x = -(pool_diameter/2):(pool_diameter/6):(pool_diameter/2)
                for y = -(pool_diameter/2):(pool_diameter/6):(pool_diameter/2)
                    if (x^2 + y^2 <= (pool_diameter/2)^2)
                        x2 = x;
                        y2 = y;
                        dist_to_wall = pool_diameter/2-sqrt(x2^2+y2^2);

                        for k = 1:N_ac
                            DC_activation = zeros(1,N_dc);
                            for i = 1:N_dc
                                DC_activation(i) = exp(-(dist_to_wall -DC(i))^2/(2*sigma_dc^2));
                            end
                    %Calculate AC activation (i.e. value of the action)
                            AC_activation = zeros(1,N_ac);
                            for i = 1:N_ac
                                for j = 1:N_dc
                                    AC_activation(i) = AC_activation(i) + ...
                                        DC_activation(j)*weights_dc(j,i);
                                end
                            end
                                if x >0
                                    central_angle = atan(y/x);
                                else
                                    central_angle = pi+atan(y/x);
                                end
                            moving_dir = pi+central_angle+k/N_ac*2*pi;
                            x2 = x + (AC_activation(k)/10)*cos(moving_dir);
                            y2 = y + (AC_activation(k)/10)*sin(moving_dir);
                            hold on;
                            line([x x2],[y y2],'Color',[k/N_ac 0 1-k/N_ac]);
                            %xlabel('day'+ day + ', trial' + trial),ylabel('DC cognitive map, 2% alteration rates');
                        end
                    end
                end
            end
            %plot the maze and platform
            plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),'k');
            plot(whichplatform*platform_x+platform_radius*cos(th),whichplatform*platform_y+platform_radius*sin(th),'k'); 
            end
            
            if (plot_cognitive_maps)
            %subplot(Ndays, Ntrials*pln,(day-1)*Ntrials*pln+trial*pln);
            subplot(1, pln, 3);
            
            hold on
            %plot the cognitive map (for PC)
            for x = -(pool_diameter/2):(pool_diameter/6):(pool_diameter/2)
                for y = -(pool_diameter/2):(pool_diameter/6):(pool_diameter/2)
                    if (x^2 + y^2 <= (pool_diameter/2)^2)
                        for k = 1:N_ac
                            PC_activation = zeros(1,N_pc);
                            for i = 1:N_pc
                                PC_activation(i) = exp(-((x - PC_x(i))^2 + (y...
                                - PC_y(i))^2)/(2*sigma_pc^2));
                            end
                    %Calculate AC activation (i.e. value of the action)
                            AC_activation = zeros(1,N_ac);
                            for i = 1:N_ac
                                for j = 1:N_pc
                                    AC_activation(i) = AC_activation(i) + ...
                                        PC_activation(j)*weights_pc(j,i);
                                end
                            end
                            x2 = x + (AC_activation(k)/10)*cos(k/N_ac*2*pi);
                            y2 = y + (AC_activation(k)/10)*sin(k/N_ac*2*pi);
                            hold on;
                            line([x x2],[y y2],'Color',[k/N_ac 0 1-k/N_ac]);
                        end
                    end
                end
            end
            %plot the maze and platform
            plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),'k');
            plot(whichplatform*platform_x+platform_radius*cos(th),whichplatform*platform_y+platform_radius*sin(th),'k'); 
            end
            end
            
            
            toc
         end
    end
    

else
    %run multiple times without plotting
    disp('else');
    for rep = 1:Nruns
        display(num2str(rep));
        weights_pc = rand(N_pc,N_ac)*Wmult;
        weights_dc = rand(N_dc,N_ac)*Wmult;

        DC = zeros(0,N_dc);
        PC_x = zeros(0,N_pc);
        PC_y = zeros(0,N_pc);
    
        for i = 1:N_dc %#For each distance cell:
        DC(i) = rand*(pool_diameter); %#Random positions of distance cells
        end
        for i = 1:N_pc %For each place cell:
        PC_x(i) = (rand - 0.5)*pool_diameter; %Random positions of place cells
        PC_y(i) = (rand - 0.5)*pool_diameter;
            while (PC_x(i)^2 + PC_y(i)^2 > (pool_diameter/2)^2) %Checks for out of bounds
                PC_x(i) = (rand - 0.5)*pool_diameter;
                PC_y(i) = (rand - 0.5)*pool_diameter;
            end
        end
        
        for day = 1 : Ndays
            for trial = 1:Ntrials
                %# fix or variable platform(not match betaval code)
                %whichplatform = 1 %# fix
                whichplatform=randi([0,1],1); %# variable

                idx = randi(4); %randomly choose one of 4 starting locations
                starting_x = starting_xs(idx);
                starting_y = starting_ys(idx);

                %No betas!

                %#run trial
                [wres_pc,wres_dc, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants,...
                    latency, speed_std, speed_ps, mean_angle,time_step]=...
                    run_trial(curr_model, weights_pc, weights_dc, Wmult, sigma_pc, ...
                    sigma_dc, sigma_ac, PC_x, PC_y, DC, Vdecay, ac_const, beta, etdecay, ...
                    alpha, gamma, Wnoise, whichplatform*platform_x, whichplatform*platform_y, ...
                    starting_x, starting_y, speed, hitwall);        
                weights_mod = wres_mod;
                weights_pc = wres_pc;
                weights_dc = wres_dc;
            
                PMs(1,day,trial,rep) = latency; 
                PMs(2,day,trial,rep) = dist; 
                if whichplatform == 1
                    PMs(3,day,trial,rep) = quadrants(4,1)*100; %#target quadrant
                    PMs(4,day,trial,rep) = quadrants(2,1)*100; %#opposite quadrant
                else
                    PMs(3,day,trial,rep) = quadrants(2,1)*100; %#target quadrant
                    PMs(4,day,trial,rep) = quadrants(4,1)*100; %#opposite quadrant
                end
                PMs(5,day,trial,rep) = wall_zone*100; %#wall zone
                PMs(6,day,trial,rep) = speed_std*100; %# speed_std
                PMs(7,day,trial,rep) = mean_angle;
                PMs(8,day,trial,rep) = time_step;
            end
        end
        toc
    end
end