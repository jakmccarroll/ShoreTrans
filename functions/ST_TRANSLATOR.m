%% ST_TRANSLATOR
% 27/8/2020, Updated May 2021
% Use: Performs the SHORETRANS profile translation
% Called by: ST_OPTIMIZER or ST_MAIN (if OPTIMIZER is switched off)
% INPUTS
    % Xi = translation distance

function [outProf, OPT] = ST_TRANSLATOR(x0, z0, Xi, dV_input, OPT)

if OPT.shortOutput == 0
    disp(['Running TRANSLATOR with estimated distance = ' num2str(Xi) ' m']);
end

%% UNPACK VARS FROM OPT
z0_temp = OPT.z0_temp;  % z0_temp = z0 + dune accretion + dune HORZ GROWTH + trend_x + (Interp over rocks)
z1 = OPT.z_raise; % z1 = z0_temp + dS; 
% X_range = OPT.X_range; % Translation range (X_range)
% X_len   = OPT.X_len; 

% NE layer
z_rock = OPT.rockLayer;
z0_rock = OPT.z0_rock; % Z0_ROCK (MAX OF UNION OF ROCK AND Z0)

% ToCr
ToCr_level  = OPT.toeCrest_level;
ToCr_level2 = OPT.toeCrest_level2; % post-SLR level
ToCr_ind    = OPT.toeCrest_ind;
ToCr_x      = OPT.ToCr_x;

% wall
wall       = OPT.wallSwitch;
wall_level = OPT.wall_level;
wall_ind   = OPT.wall_ind;
wall_x     = OPT.wall_x;

% SL and DoC
dS = OPT.dS;
S_initial = OPT.S_initial;
S_final   = OPT.S_final;
DoC       = OPT.DoC;
DoC_ind   = OPT.DoC_ind;
DoC2      = OPT.DoC2;
DoC2_ind   = OPT.DoC2_ind;

% Dune
%

%% PERFORM PROFILE TRANSLATION (LOOP - each x-shore position)
% Loop through the profile translation routine for each position in X_range
outProf=[];

% for n = 1 : X_len % n = shift onshore (in metres as dx=1)
% n = 1
%   	disp(num2str(n))

%% PROFTRANS 1: Set iteration of translation distance (Xi)
% Xi               = X_range(n);
%     outProf(n,1).z0    = z0;
outProf.Xi    = Xi;
z_temp = z1; % z_temp is the working profile the the loop iterations

% If wall = ON, set an initial profile to determine the...
% ... onshore-of-wall erosion volume
if wall==1
    z_noWall = z_temp;
end

%% PROFTRANS 2 - SHIFT PROFILE 
% active_ind -> the section that SHIFTS UP (by SLR) and ITERATES ONSHORE (by Xi)
active_ind = [ToCr_ind + Xi : DoC_ind - 1]'   ;

if active_ind(1) < 0
    disp('active_ind is NEGATIVE !!!!!!!!!!!!!!!!')
end

z_temp(active_ind) = z1(active_ind - Xi)  ;  % shift prof on/off
                %... (permitting erosion below rock + beyond wall)

if wall==1 && OPT.wall_overwash == 0
    % if wall=ON, translate the noWall profile (as per z_temp)...
    z_noWall(active_ind) = z_temp(active_ind);

    % ... then shift up z_temp where it is onshore of the wall
    % ... (but this caused an issue (can't remember) and had to switch off
    % ... update (9/1/2020). When switched on it prevents the slump erosion at top of cliff (e.g. PPT N DUNES)
    % ....upd(9/1/2020). Workaround was to put this at the start of the SLUMP section
    % ...upd(31/1/2020). Turned it BACK ON. Going in circles, need a solution that is tested on...
    % WALL prof + SLUMP at top of cliffs prof
%         z_noWall2 = z_temp; % need this for SLUMP TOP FIX (2/2/2020)
    % Update 10/4/2020 - Slumping occurs onshore of wall, had to switch off slumping (DFIG10, ConcMod)
    z_temp(x0<=x0(wall_ind) ) = z0_temp(x0 <= x0(wall_ind) );
end

if OPT.rockSwitch == 1 
    % if rock=ON, shift up sections of profile that fall below the rock layer
    rock_ind         = find(z_temp < z_rock);
    z_temp(rock_ind) = z_rock(rock_ind); 

end

if OPT.rollover == 0 % if rollover=OFF, prevent elevation increase behind the crest/toe
    % (21/5/2020) removed "  && wall==0 " 
    % ... this allows for accretion between the wall and the crest
    % ... it's intended for sheet piling / rock armour (e.g., Torcross) that can be buried,...
    % ... but with a taller wall behind, that can't be overwashed.
    noRoll_ind = find(z_temp > z0_temp & x0 <= x0(ToCr_ind));
    z_temp(noRoll_ind) = z0_temp(noRoll_ind);

elseif OPT.rollover==2 % && wall==0 % ROLLOVER = 2 -> barrier maintains height
    noRoll_ind = find(z_temp>z0_temp & x0<x0(ToCr_ind));
    z_temp(noRoll_ind) = OPT.toeCrest_level;

end

%% PROFTRANS 3 -  OFFSHORE SMOOTHING -> LINEAR INTERP FROM OFF_PTS TO DoC
% (contains necessary ASSUMPTIONS)

W = x0(DoC_ind) - x0(ToCr_ind);
st_min = round(.10 * W);
% st_min = 100;

% START (st) point for smoothing
if Xi <= -st_min
    st = DoC_ind + Xi - 1 ; % Xi is negative, so "st" is onshore of DoC
    % if profile is RECEEDING, smooth the translation distance...
elseif Xi > -st_min && Xi <= 0   
    st = DoC_ind - st_min; 
    % need the "50" to avoid a big step when there high SLR with low translation distance
    
% elseif Xi > 0 && Xi < 50
%     % if profile is PROGRADING, smooth (1x) the translation distance
% %         st = DoC_ind - 2*Xi;   % Xi is positive, so "st" is onshore of DoC 
%     en = DoC_ind + 50         ;     % Xi is positive, so "st" is OFFSHORE of DoC 
%     
elseif Xi >  0
    en = DoC_ind + Xi;
end
    
% END (en) point for smoothing
if Xi <= 0 
    en = DoC_ind  ; 
elseif Xi > 0 
%         en = st;
    st = DoC_ind - 1 ;
end

OPT.smooth_st = st  ;
OPT.smooth_en = en  ;

% Apply smoothing
if Xi <= 0 
    z_temp(st+1:en) = interp1([x0(st), x0(en)], [z_temp(st), z0(en)], x0(st+1:en));
elseif Xi > 0
    z_temp(st+1:en) = interp1([x0(st), x0(en)], [z_temp(st), z0(en)], x0(st+1:en));
end

% If wall=ON, copy the smoothing to the "no wall" profile
if wall == 1
    z_noWall(st+1:en) = z_temp(st+1:en);
end 

%% PROFTRANS 4 - ONSHORE CREST RAISE UP (when profile shifting OFFSHORE)
% As the profile marches offshore, this maintains the crest of the barrier
% top_ind and top_level refer to barrier crest in this instance
if Xi>0
    st = ToCr_ind + 1;
    en = ToCr_ind + Xi;
    z_temp(st:en) = ToCr_level + dS;

    if OPT.duneHGrow == 1 % Maintain DUNE HORZ GROWTH CREST HT (1/8/2020) 
        ind1 = OPT.duneHGrow_ind_crest ;
        ind2 = OPT.duneHGrow_ind_offPt + Xi;
        z_temp(ind1:ind2) = OPT.duneHGrow_crest + S_initial;
    end

end

%% PROFTRANS 5 -  DUNE SLUMPING (if ROLLOVER = OFF || SLUMPCHECK=ON)
% if rollover is OFF... and duneSlump is ON
% ... the dunes get eaten into and this is the angle of ...
% ... repose of the dune post-translation duneface.

z_noSlump=z_temp;
if (OPT.rollover==0 || OPT.slumpCheck ==1) && OPT.duneSlump==1

    [z_temp, SLUMP, z_maxSlump]=...
    ST_SLUMP(x0, z0_temp, z_temp, OPT);

    z_slump = z_temp;
    OPT.z_maxSlump = z_maxSlump;   

end

%% PROFTRANS 6 - ROLLOVER - BACK DUNE SLOPE (roll_backSlope) (!!!NEEDS SEPARATE FUNCTION!!!)
% NEEDS TO GO IN A FUNCTION!!!
% creates a (low-gradient) slope behind a barrier that is rolling back
if OPT.rollover > 0 %
    z_noWash = z_temp;

    washON=1;
    ind_wash_0 = ToCr_ind + Xi; % initial back of barrier pt
    ind_wash   = ind_wash_0; % this iterates to march onshore

    z_back = z_temp(ind_wash); % back of barrier elevation
    z_step = OPT.dX .* tand(OPT.roll_backSlope); % step down this much each dx at back of barrier
    z_step_n = z_step; % this iterates to decrease by z_step each step

    while washON==1
        ind_wash = ind_wash - OPT.dX; % iterate onshore
        z_temp(ind_wash) = z_back - z_step_n; % step down the profile
        z_step_n = z_step_n + z_step; % increase step size for next iteration

        if z_temp(ind_wash) <= z_noWash(ind_wash)
            z_temp(ind_wash) = z_noWash(ind_wash); 
                % if overwash has dipped below existing profile, bring it back up
            washON=0; % break the loop
        end

    end
end

%% PROF TRANS 7 - WALL -> CALC BEHIND-THE-WALL VOL and REDISTRIBUTE OFFSHORE
% this is the (hypothetical) volume eroded from "onshore of the wall" (or cliff)
% that has to be shifted offshore of the wall to the real profile

if wall == 1
    z_noRedist = z_temp;
    [z_temp, WALL]=...
        ST_WALL_VOL(x0, z0_temp, z_temp, z_noWall, OPT, Xi);
    outProf.WALL = WALL;
    outProf.below_z_min = WALL.below_z_min;    
end

%% PROF TRANS (7a) -> LAST ROCK CHECK
% Added 9/6/2020, fix for issue where bed drops below rNE layer
if OPT.rockSwitch == 1
   ind = find(z_temp <= z_rock);
   z_temp(ind) = z_rock(ind);
end

%% PROF TRANS 8 -  VOL CALCS, RECESSION CALCS
% Final z-profile
outProf.z_final = z_temp;
outProf.z_NE    = z_rock;

% Max of inital prof and NE-layer
outProf.z0_rock = z0_rock;

% Interim z-profiles (1-WALL)
if wall==1
    outProf.z_noWall   = z_noWall;
    outProf.z_noRedist = z_noRedist;
else
    outProf.z_noWall   = nan .* zeros(size(z0));
    outProf.z_noRedist = nan .* zeros(size(z0));   

end

% Interim z-profiles (2-SLUMP, duneSlope)
if OPT.rollover==0
    outProf.z_noSlump = z_noSlump;
else
    outProf.z_noSlump = nan .* zeros(size(z0));
end

% Interim z-profiles (3 - OVERWASH, backSlope)
if OPT.rollover==1
    outProf.z_noWash = z_noWash;
else
    outProf.z_noWash = nan .* zeros(size(z0));           
end

% Interim -> if SLUMP on 
% NEED TO FIX THIS
% INTERIM Z_PROFS SHOULD BE INCLUDING SOMETHING (E.G. Z_SLUMP)
% ... NOT EXCLUDING (Z_NoSLUMP)
if exist('z_slump')
    outProf.z_slump = z_slump;
end

% shorelines
 % removed the "last +1" (15/6/2020)
 % this method for shoreline finding is inconsistent with "ST_getContour" in "ST_timeSeries_v2", 
 %... producing an offset error in R_sh_total of -1 (8/4/2021)
ind1 = find(OPT.z0_rock > S_initial, 1, 'last')  ; % 20/2/2020 changed from 'first' to 'last'...
ind2 = find(outProf.z_final>S_final, 1, 'last')  ; % ... to account for profs with lagoons (<0 at back)
outProf.shore_x1 = x0(ind1)   ;                         % ... added "+1" so it will still be first point after going under SL
outProf.shore_ind1 = ind1     ;                    
outProf.shore_x2 = x0(ind2)   ; 
outProf.shore_ind2 = ind2     ; 
outProf.R_sh = outProf.shore_x2 - outProf.shore_x1;

% ToeCrest Lines
    % Added 15/6/2020. Inconsistency in using 'last + 1' for shoreline, but 'last' for ToCr
    % need to make method same throughout (15/6/2020)
ind1 = ToCr_ind;
ind2 = find(outProf.z_final >= ToCr_level2, 1, 'last');

if isempty(ind2) % for profiles where the toeCrest can't rise (e.g., behind a wall)
	ind2 = find(outProf.z_final >= ToCr_level, 1, 'last');
        % if new ToeCrest level not met, revert to initial ToeCrest level
    
    if isempty(ind2) % if still lower, take highest point in new profile
        max2 = max(z_temp);
        ind2 = find(z_temp >= max2, 1, 'last');
    end
end
ToCr2_ind = ind2;

outProf.ToCr_x1 = ToCr_x;                        
outProf.ToCr_ind1 = ind1;
outProf.ToCr_x2 = x0(ind2); 
outProf.ToCr_ind2 = ind2; 
outProf.R_ToCr = outProf.ToCr_x2 - outProf.ToCr_x1;

% VOLUME calcs (18/2/2020 -> changed z0 to z0_rock for vol calcs.
    % ... needed to handle offshore reefs
    % ... z0_rock also accounts for the z0_temp changes (DuneAccrete and DeepOn)
    % (20/2/2020) changed all DoC to DoC2
V0 = trapz(x0(1:DoC2_ind), z0_rock(1:DoC2_ind) - DoC2) ; % offset by DoC so no NEG VOLS
V2 = trapz(x0(1:DoC2_ind), z_temp(1:DoC2_ind) - DoC2)  ; % (19/2/2020) changed to (- DoC) previous was (+ DoC)
outProf.dV  = V2 - V0;
%     dV_n(n,1) = outProf.dV; 
outProf.dV_error = outProf.dV - dV_input  ;

% DUNE VOLUME
V0_DUNE = trapz(x0(1:ToCr_ind), z0_rock(1:ToCr_ind)); 
V2_DUNE = trapz(x0(1:ToCr_ind), z_temp(1:ToCr_ind)); 
outProf.dV_DUNE = V2_DUNE - V0_DUNE;

% SHOREFACE VOL
V0_SF = trapz(x0(ToCr_ind:DoC2_ind), z0_rock(ToCr_ind:DoC2_ind) - DoC2); % offset by DoC so no NEG VOLS 
V2_SF = trapz(x0(ToCr_ind:DoC2_ind), z_temp(ToCr_ind:DoC2_ind) - DoC2); 
outProf.dV_SF = V2_SF - V0_SF;

% BeachWidth
outProf.BchWidth_1 = round(outProf.shore_x1 - outProf.ToCr_x1 - 1); 
    % "-1" to avoid a totally depleted wall from having width = 1m
outProf.BchWidth_2 = round(outProf.shore_x2 - outProf.ToCr_x2 - 1);
outProf.BchW_ratio = outProf.BchWidth_2 / outProf.BchWidth_1;

% Wall Beach Width (wall position to shoreline)
outProf.WallBchWidth_1 = round(outProf.shore_x1 - OPT.wall_x ); 
    % "-1" to avoid a totally depleted wall from having width = 1m
outProf.WallBchWidth_2 = round(outProf.shore_x2 - OPT.wall_x );
outProf.WallBchW_ratio = outProf.WallBchWidth_2 / outProf.WallBchWidth_1;

% SHOREFACE WIDTH and HEIGHT
outProf.W_SF1 = x0(DoC_ind) - x0(outProf.ToCr_ind1);
outProf.W_SF2 = x0(DoC_ind) - x0(outProf.ToCr_ind2);

outProf.h_SF1 = z0(outProf.ToCr_ind1) - z0(DoC_ind);
outProf.h_SF2 = z0(outProf.ToCr_ind2) - z_temp(DoC_ind);

%% CHECK IF CONDITION FOR EXCESS EROSION OFFSHORE OF WALL HAS BEEN MET

if OPT.wall_z_min_check == 1
    z_off_wall = z_temp(OPT.wall_ind + 2);
    
    if z_off_wall <= OPT.wall_z_min
        disp('wall REACHED');
        OPT.wall_z_min_REACHED = 1;
%         OPT.BREAKER = 1;

    else
        OPT.wall_z_min_REACHED = 0;
    end
    
else
        OPT.wall_z_min_REACHED = 0;
%         OPT.BREAKER = 0;
end


%% SHORT OUTPUT (SAVES MEMORY)
if OPT.shortOutput == 1
    try outProf = rmfield(outProf,'WALL'); end
    try outProf = rmfield(outProf,'z0_rock'); end
    try outProf = rmfield(outProf,'z_noWall'); end
    try outProf = rmfield(outProf,'z_noRedist'); end
    try outProf = rmfield(outProf,'z_noSlump'); end
    try outProf = rmfield(outProf,'z_noWash'); end
    try outProf = rmfield(outProf,'z_slump'); end


end






%%














%%
