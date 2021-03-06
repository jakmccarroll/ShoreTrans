%% ST_OPTIMIZER
% !!!!! -> This function is horrendous, needs a complete rewrite !!!!
% Update May 2021
% Use: To optimise (minimise) the numer of profile translation distance calculated ...
    % ... to find balance the volume (minimise the error).
% Called by: ST_MAIN

function [Xi_winner, Xi_ind, outPROFS, OPT] = ST_OPTIMIZER(x0, z0, dV_input, OPT)

outPROFS = [];
if OPT.shortOutput == 0
    disp(['Finding OPTIMISED estimate of TRANSLATION DISTANCE'])
end

%% FIRST-PASS ESTIMATE
% SLR TRANSLATION -> BRUUN ESTIMATE
% R = dS * (W/h)
dS = OPT.dS   ;
h  = OPT.toeCrest_level - OPT.DoC   ;
W  = x0(OPT.DoC_ind) - OPT.ToCr_x   ;
Xi_est_a = round(-dS*(W/h))         ;

% VOLUME CHANGE TRANSLATION ESTIMATE
% R = dV/h
Xi_est_b = round(dV_input/h)   ;

% INITIAL ESTIMATE
Xi_est = Xi_est_a + Xi_est_b   ;

if -Xi_est >= OPT.toeCrest_ind   % catches too big initial estimates, and sets to max possible for prof extent
% if Xi_est < -200
    disp('CATCH 1: -Xi_est >= OPT.toeCrest_ind !!!!!!!!!!!!!')
    Xi_est = -(OPT.toeCrest_ind - 1);
end

% First-pass translation
[outPROFS1, OPT] = ST_TRANSLATOR(x0, z0, Xi_est, dV_input, OPT)  ;
outPROFS = [outPROFS; outPROFS1]     ;

%% LOOP TO OPTIMIZE
dV_gettingClose = 25  ;
dV_error = outPROFS1.dV_error ;
if abs(dV_error) > 0
    optimized = 0;
elseif dV_error == 0
    optimized = 1;
end

dV_divider = round(h) ;

n = 1;
while optimized == 0

    if abs(dV_error) > dV_gettingClose && OPT.wall_z_min_REACHED == 0
        Xi_est = Xi_est + round(-(dV_error/dV_divider)) ; % CAN CHANGE THE DENOMINATOR TO TWEAK
        
        if ~ismember(Xi_est, [outPROFS.Xi]')
            
%             if -Xi_est >= OPT.toeCrest_ind   % catches too big initial estimates, and sets to max possible for prof extent
%                 % if Xi_est < -200
%                 disp('CATCH 2: Xi_est is less than limit set !!!!!!!!!!!!!')
%                 Xi_est = -(OPT.toeCrest_ind - 1);
%             end
            
            [outPROFS1, OPT] = ST_TRANSLATOR(x0, z0, Xi_est, dV_input, OPT);

            outPROFS = [outPROFS; outPROFS1]; %#ok<*AGROW>
        end
        dV_error = outPROFS1.dV_error;
        dV_divider = dV_divider + 1;
%         dV_divider = max([dV_divider, h]);

    elseif abs(dV_error) <= dV_gettingClose && OPT.wall_z_min_REACHED ~= 1
        for i = [-2,-1,1,2]
            Xi_est1 = Xi_est + i; 
            if ~ismember(Xi_est1, [outPROFS.Xi]')
                [outPROFS1, OPT] = ST_TRANSLATOR(x0, z0, Xi_est1, dV_input, OPT);
                outPROFS = [outPROFS; outPROFS1];
            end
        end
        optimized = 1; % END THE LOOP
        OPT.opt_outcome = 1;

    elseif OPT.wall_z_min_REACHED == 1  % WALL Z MIN REACHED
        while OPT.wall_z_min_REACHED == 1
            
            Xi_est = Xi_est + 1 ;
            [outPROFS1, OPT] = ST_TRANSLATOR(x0, z0, Xi_est, dV_input, OPT);
            outPROFS = [outPROFS; outPROFS1];
            disp('Z MIN REACHED in OPTIMIZER');
            disp(['OPT.wall_z_min_REACHED = ' num2str(OPT.wall_z_min_REACHED) ]);
        end
        optimized = 1; % END THE LOOP 
        OPT.opt_outcome = 2;
        
    elseif dV_error < 0 && OPT.wall_z_min_REACHED == 1    
        Xi_est = Xi_est + 5 ;
        [outPROFS1, OPT] = ST_TRANSLATOR(x0, z0, Xi_est, dV_input, OPT);
        outPROFS = [outPROFS; outPROFS1];

    end

    % CATCH ENDLESS LOOPS
    n = n + 1;
    if n > 100
        for i = [-5,-4,-3,-2,-1,1,2,3,4,5]
            Xi_est1 = Xi_est + i; 
            if ~ismember(Xi_est1, [outPROFS.Xi]')
                [outPROFS1, OPT] = ST_TRANSLATOR(x0, z0, Xi_est1, dV_input, OPT);
                outPROFS = [outPROFS; outPROFS1];
            end
        end
        optimized = 1;
%         OPT.opt_outcome = 0;
        warning('Optimization reached max iterations');
    end
end


%% GET THE WINNER

if OPT.wall_z_min_check == 0
    [dV_error, Xi_ind] = min(abs([outPROFS.dV_error]));
    Xi_all = [outPROFS.Xi]';
    Xi_winner = Xi_all(Xi_ind);
    
    if dV_error <= h
        OPT.opt_outcome = 1;
    elseif dV_error > h
        OPT.opt_outcome = 0;
        warning('Large dV_error, check results');
    end
        
elseif OPT.wall_z_min_check == 1

    dV_ERRORS = [outPROFS.dV_error]'  ;
    below_z_min = [outPROFS.below_z_min]'  ;
    ind1 = find(below_z_min == 0)  ;
    
    [dV_error, ind2] = min(abs(dV_ERRORS(ind1))  )   ;
    Xi_ind = ind1(ind2)   ;
    
    Xi_all = [outPROFS.Xi]'     ;
    Xi_winner = Xi_all(Xi_ind)  ;

end




%%










%%


