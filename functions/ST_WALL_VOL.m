function [z_temp, WALL] = ST_WALL_VOL(x0, z0, z_temp, z_noWall, OPT, Xi)
%% CALC BEHIND WALL VOL and APPLY TO OFFSHORE
%
% The difference volume (dV_behindWall) is the area...
    % ... ONSHORE of WALL, and ...
    % ... BELOW the initial Toe/Crest level (post-SLR)
    % dV_behindWall = V_ToCr_Block - V_noWall;
    % ... where V_noWall is the up-on translated profile, ignoring the wall
%
% [z_out, z_noWall_out]= ST_WALL_VOL...
            % ... (x0, z0, z_temp, z_noWall, OPT, Xi)
%
% INPUTS
    % x0 = input profile x-values
    % z0 = input profile z-values
    % z_temp = working profile -> translation applied (by other function,
        % e.g., TRANSLATOR) in front of wall, surface behind wall lifted to z0
    % z_noWall = working profile --. translated, with zone behind wall
        % allowed to erode 
    % OPT -> structure with model settings (see OPT =  ST_OPT_defaults)  
    % Xi -> translation distance z0 translated, to get to z_temp (and z_noWall)
%
% OUTPUTS 
    % z_out = output profile with Beuzen rule applied
    % WALL = structure with all outputs
%

%% Extract variables from OPT structure
ToCr_level  = OPT.toeCrest_level;
ToCr_level2 = OPT.toeCrest_level2;
ToCr_ind    = OPT.toeCrest_ind;
DoC         = OPT.DoC;
DoC_ind     = OPT.DoC_ind;

wall_ind      = OPT.wall_ind;
wall_level    = OPT.wall_level;
redist_ratio  = OPT.redist_ratio;
wall_z_min    = OPT.wall_z_min;
z_min_initial = OPT.z_min_initial;
wall_no_erode_behind = OPT.wall_no_erode_behind;

rock        = OPT.rockSwitch;
z_rock      = OPT.rockLayer;
dx          = x0(2) - x0(1);

dS          = OPT.dS;

%% BEHIND WALL EROSION CHECK
% added 16/6/2020 -> needs to put in "OPT.wall_no_erode_behind"
if wall_no_erode_behind == 1
    ind = find(x0 < x0(wall_ind) & z_temp < z0);
    z_temp(ind) = z0(ind);
end

%% VOL under NO WALL surface, from x=0 -> x = WALL_IND = CLIFF (z_rock < maxErodeLevel)

% 1. BOTTOM SURFACE -> z_noWall (capped at ToCr2)
z_noWall_cap = z_noWall; % puts a cap on z_noWall at toe/crest level post-SLR

z_min = min(z0, repmat(ToCr_level2, size(z0))); % a surface that is the 
    % min value of either z0 or ToCr_level2

% ADDED/CHANGED THIS STEP TO RAISE CAP WHEN BACK OF PROFILE IS FLAT (TORCROSS)
%  CHECK THIS DOESN'T BREAK PPT 
% z_noWall_cap(z_noWall > top_level2) = top_level2;
% z_noWall_cap(1:wall_ind + Xi) =  ToCr_level2; % Draft5 to MG
z_noWall_cap(1:wall_ind + Xi) =  ToCr_level2; % checked, not changed, Feb 2021, MG Rev1

ind = find(z_noWall_cap>ToCr_level2 & x0<x0(wall_ind)); % extra check to get all pts above ToCr2_level
z_noWall_cap(ind) = ToCr_level2;

V_noWall = trapz(x0(1:wall_ind), z_noWall_cap(1:wall_ind) + DoC) ; 
    % vol above DoC to avoid neg values cancelling

% 2. TOP SURFACE -> TOE/CREST-BLOCK
% VOL under "TOE/CREST-BLOCK" a rectangular block with top surface = ToCr_level + dS
    % + dS added Feb 2021
% take volume from x=0 -> x=wall_ind
z_ToCr_Block = ToCr_level2 .* ones(size(x0)) ; % checked, not changed, Feb 2021, MG Rev1


if wall_ind + Xi <= 0 
    error('wall_ind + Xi < 0. Cannot translate back far enough! Increase profile onshore length.')
end
ind = find(z0 <= ToCr_level2 & x0 > x0(wall_ind + Xi) & x0<x0(wall_ind)) ; % this checks that if ...
    % section of the top surface within the extent of the translation distance (Xi) are 
    % lower than ToCr2, they are reduced to z0. 
    % e.g., avoids over counting for a Brunn-type profile.
z_ToCr_Block(ind) = z0(ind); 
V_ToCr_Block = trapz(x0(1:wall_ind), z_ToCr_Block(1:wall_ind)  + DoC);

% 3. TOP SURFACE - BOTTOM SURFACE = dV BEHIND WALL
% The difference volume (dV_behindWall) is the area...
    % ... ONSHORE of WALL, and ...
    % ... BELOW the initial Toe/Crest level (post-SLR)

dV_behindWall = V_ToCr_Block - V_noWall ;

%% PUT OUTPUTS IN WALL STRUCTURE
WALL=[];
WALL.z_noWall = z_noWall;
WALL.z_noWall_cap = z_noWall_cap;
WALL.z_ToCr_Block = z_ToCr_Block;
WALL.dV_behindWall = dV_behindWall;

%% REDISTRIBUTE LOST VOL OFFSHORE

if dV_behindWall > 0 && z_min_initial == 0 
    
    % find the first point that is offshore of the wall 
    % (and above the rock layer, if on)
    if rock == 1
        ind_st = find(x0 >= x0(wall_ind) & z_temp > z_rock, 1, 'first') ;
    elseif rock == 0
        ind_st = find(x0>=x0(wall_ind), 1, 'first') ;
    end
%     ind_en=find(z_temp>z0 & x0>ind_IEL & z_temp>z_rock, 1, 'first');
    len_wall2DoC = round( x0(DoC_ind) - x0(ind_st) )  ; 
        % distance from 'First Erodible Pt' to DoC
    
    if len_wall2DoC > 0
        len_wall2redistEndPt = round(redist_ratio .* len_wall2DoC)   ;

        ind_en  = round( x0(ind_st) + x0(len_wall2redistEndPt ) )   ;
        ind_len = round( x0(ind_en) - x0(ind_st) + 1  )  ;
        IND = ind_st : ind_en  ;

        if length(IND) > 1
            % linear interp VOL REDIST from WALL to cross-over point 
            % using area of a triangle 
            % A = xz / 2 --> z=2A/x 
            % (here, z=elevation to erode at the wall, tapering to 0 at redist_ratio offshore)
            WALL.V          = dV_behindWall ;
            WALL.x          = ind_st:ind_en ;
            WALL.IND        = IND ;
            WALL.x_len      = ind_len ;
            WALL.z0         = (2.*WALL.V) ./ WALL.x_len ;
            WALL.z_end      = 0  ; % taper to zero change
            WALL.z_off      = interp1([WALL.x(1) WALL.x(end)],[WALL.z0, WALL.z_end],WALL.x)'   ;

            % REDISTRIBUTE the hypothetical erosion to offshore of the wall
            z_temp(IND) = z_temp(IND) - WALL.z_off   ;
        end
        
    else
        disp('Fully depleted, no wall redistribution possible');
        % z_temp(1:ind_st) = z_rock(1:ind_st);  % RAISE to ROCK  
            % removed/skipped this line, can't see how it helps, and it appears
            % to be in error as it applies z_rock to profile where NE-layer is switched off
        IND = [];
    end
    
elseif  dV_behindWall <= 0 && z_min_initial == 0
    IND=[];

elseif  z_min_initial == 1
    z_temp  =  z0;
    
end


%% ROCK CHECK (new on 18/5/2020)
% much simpler check that raises any wall vol redistibution, ...
% ... that drops it below the rock, back up to rock level
if rock == 1 && z_min_initial == 0
    % if ROCK = ON and max erosion had not already been reached...
    ind = find(z_temp < z_rock );
    z_temp(ind) = z_rock(ind);
end


%% ROCK CHECK and REDISTRIBUTE (OLD ITERATIVE VERSION -> REPEAT UNTIL DONE, SLOW, lots of issues)
% DROPPED ON 18/5/2019
% iteratively check if any profile points drop below the rock level
% raise these points above the rock, redistibrute, then repeat...

% try %!!! NEED TO FIX, this section will now skip 

% if ~isempty(IND) && rock==1
% 
%     ind_rock = find(z_temp(IND) < z_rock(IND)) ;
%     IND_rock = IND(ind_rock)' ;
% 
%     ind_noRock = find(z_temp(IND) >= z_rock(IND)) ;
%     IND_noRock = IND(ind_noRock)'  ;
% 
%     if length(IND_rock)>1
%         while length(IND_rock)>1
% 
%             % find ROCK-BLOCKED VOL
%             vol_temp = trapz(x0(IND_rock), z_temp(IND_rock) - DoC ) ;
%             vol_rock = trapz(x0(IND_rock), z_rock(IND_rock) - DoC ) ;
%             dV_temp  = vol_rock - vol_temp ;
% 
%             % RAISE to ROCK
%             z_temp(IND_rock) = z_rock(IND_rock) ;
% 
%             % linear interp ROCK VOL REDIST from WALL to cross-over point 
%             % EXCLUDING ROCK points (IND_noRock)
%             % A = xz / 2 --> z=2A/x 
%             WALL.V = dV_temp ; 
%             WALL.x = IND_noRock ;
%             WALL.IND = WALL.x ;
%             WALL.x_len = round(length(WALL.IND) .* dx) ;
%             
%             if WALL.x_len > 1 % catch error that occurs if profile is fully depleted (18/5/2020)
%                 WALL.z0 = (2.*WALL.V) ./ WALL.x_len   ;
%                 WALL.z_end = 0 ; % taper to zero change
%                 WALL.z_off = interp1([1 WALL.x_len],[WALL.z0, WALL.z_end],[1 : WALL.x_len])' ;
% 
%                 % REDIST
%                 z_temp(IND_noRock) = z_temp(IND_noRock) - WALL.z_off ; % ERROR ON THIS LINE
%                     % WALL.z_off length = 81
%                         % length WALL.z_off = length(WALL.IND) .* dx
%                     % IND_noRock length = 82
%                         % length IND_noRock = find pts where z_temp >= z_rock
%                     % SOLUTION -> needed to ROUND (length(WALL.IND) .* dx)
% 
%                 % RE-INDEX
%                 ind_rock = find(z_temp(IND_noRock) < z_rock(IND_noRock));
%                     % note using IND_noRock as new rock areas may be exposed
%                 IND_rock = IND(ind_rock);
% 
%                 ind_noRock = find(z_temp(IND_noRock) >= z_rock(IND_noRock));
%                 IND_noRock = IND(ind_noRock);
%             else
%                 disp(['Profile fully depleted to DoC,' 10 ...
%                     'no further wall redistribution possible'])
%                 IND_rock = 0;
%             end
% 
%         end
%     end
% end


%%

% catch
%     warning('WALL FUNCTION REDISTRIBUTE FAIL!!! Need to fix');
% %     warning(['ST_WALL_VOL -> "ROCK CHECK and REDISTRIBUTE" ' 10 ...
% %         'section ***FAILED TO RUN***. This is an undiagnosed bug which' 10 ...
% %         ' can be ignored if there are no isolated rocks in active profile']);
% end

%% CATCH ANY EROSION BELOW DoC, RAISE BED BACK TO DoC
if  z_min_initial == 0 % if max erosion had not already been reached...
    ind = find(z_temp < DoC & x0 <= x0(DoC_ind) );
    z_temp(ind) = DoC;
end

%% TAG PROFILES THAT ERODE BELOW OPT.WALL_Z_MIN

if  z_min_initial == 0 % if max erosion had not already been reached...
    z_offshore_wall = z_temp(OPT.wall_ind + 2); % added + 1 to go offshore of wall (9/6/2020)
    if z_offshore_wall < OPT.wall_z_min
        WALL.below_z_min = 1;
    else
        WALL.below_z_min = 0;
    end

elseif z_min_initial == 1
    WALL.below_z_min = 1;
    
    
end















%%



