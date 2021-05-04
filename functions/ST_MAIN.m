%% ST_MAIN
% May 2021, Jak McCarroll
% Use: Run SHORETRANS shoreface translation model
%
% INPUTS
    % x0, z0 -> profile coords, x-positive is offshore
    % dV_input -> sediment budget change (m3/m)
    % varargin (see SETTINGS below)
% OUTPUTS
    % outProf -> all outputs to a single structure,...
    % OPT -> returns initial settings, with modifications
   
function [outProf, outPROFS, OPT] = ST_MAIN(x0, z0, dV_input, varargin)

%% SETTINGS (VARARGIN, use SETPROPERTY from openEarthTools)
% Translation increment/distance
% x = positive offshore
%%
OPT.dX    = 1;   % increment to translate profile off-/onshore (must be positive)
OPT.X_off = 0;    % max distance to translate profile OFFSHORE (progradation, default = 0)
OPT.X_on  = -100;  % max distance to translate profile ONSHORE (recession, default = -100)
OPT.X_dir = -1;   % trend profile movement -> -1=Receeding; +1=Prograding

% SLR options
OPT.dS = 1; % sea-level rise magnitude (default = 1 m)
    % set to 0 if no SLR
OPT.S_initial = 0; % initial sea-level (default = 0 m)
    % This will be non-zero when doing a second-stage of translation
OPT.DoC  = -10; % (UPPER) Depth of closure (default = 10 m)
OPT.DoC2 = []; % (LOWER) Depth of Closure 2 (default = if empty, set to DoC * 2)

% Profile shape options
    % provide either a z-level or index (relative to x0,z0)
    % for the position of the dune toe OR barrier crest
OPT.toeCrest_level=[];
OPT.toeCrest_ind=[];
    % If both _level and _ind are given, these are used.
    % If no input is given, the toe/crest ...
    % ... will be taken as the max elevation of the profile.   
    % wallLevel and toeCrestLevel can be equal
    % wall can be behind the toe (e.g. cliff buried behind dune face),...
    % ... OR wall can be in front of barrier (e.g. Torcross, where ...
    % ... sheet piling wall is offshore of top of berm/crest when moderately accreted)
OPT.toeCrest_fixedLevel = 0; % 1 = toeCrest elevation does NOT rise with SL ...
    
    
OPT.rollover = 0; % use rollover (1 = ON, 0 = OFF, 2=ON (but don't keep up with SLR)
    % (ROLLOVER ON), barrier crest rolls back, maintaining height,...
    % ... and/or rising with sea-level. 
    % (ROLLOVER OFF), erosion will eat back into dune/barrier,...
    % ... crest height not mainteined or raised with SLR
    % (ROLLOVER = 2), barrier rolls back but is restricted to current absolute height,
    % ... will NOT maintain height relative to SLR

OPT.roll_backSlope = 4; % if rollover is ON... (default = 2deg)
    % ... this is the angle (in degrees) of the onshore slope, ...
    % ... behind the the transformed barrier crest.
    
OPT.duneSlump = 1; % use DUNE SLUMPTING (1=ON, 0=OFF)
OPT.duneSlope = 30; % if rollover is OFF... (default = 30deg)
    % ... the dunes get eaten into and this is the angle of ...
    % ... repose of the dune post-translation duneface.
OPT.slumpCap   = 10; % when applying slumping, ignore pts above this level ...
    % ... this avoids slumping of steep cliffs / dunes that are well above the active profile.

OPT.duneAccrete   = 0;   % use DUNE ACCRETION (1=ON, 0=OFF)
OPT.duneAcc_Vol   = 100;   % VOLUME to accrete dunes
OPT.duneAcc_Dist  = 100; % DISTANCE (behind ToCr) over which to apply the DUNE VOL CHANGE

OPT.duneHGrow = 0;       % (1=ON, 0=OFF)  FORCE dunes to GROW HORIZONTALLY (take sediment from active shoreface)
OPT.duneHGrow_m = 0;     % METRES of HORZ DUNE GROWTH
% OPT.duneHGrow_V = 0;   % VOLUME of DUNE HORIZONTAL GROWTH
OPT.duneHGrow_crest = 0; % CREST of dune (ht above present SL) to grow from

% Rock layer options
OPT.rockSwitch=0; % rock switch (1 = ON, 0 = OFF [default])
    % RockLayer ON -> a check is made to ensure the translated profile,...
    % ... stays above the rock layer (see next input, OPT.Rock_z)
    
OPT.rockLayer = z0 - 100; % CHANGE 31/1/2020, create DUMMY ROCK LAYER
    % ... to solve issues (e.g. in _SLUMP function) requiring rocks to run
% OPT.rockLayer=[]; % rock layer (vector, same size as z0)
    % If RockSwitch is ON, RockLayer is the rock profile

% Wall-redistribution ("imaginary erosion" zone behind wall/cliff)
    % using Beuzen2018
OPT.wallSwitch = 0; % erode "imaginary" zone behind a wall (ON=1,OFF=0 [default])
    % If a wall (or cliff) exists, erode into the area
    % ... behind the wall, then transfer ...
    % ... the eroded volume offshore of the wall.
    % ... NO EROSION CAN OCCUR ONSHORE of the WALL.
    % ROLLOVER=ON does not work with WALL=ON
% If wallSwitch = ON, provide either...
    % 1. wallLevel (z-level) -> wall starts when profile first dips below this point)
    % 2. wall_ind  (x0,z0 - index) -> wall starts at this position
    % If no input, max elevation of the profile will be used.
OPT.wall_level = nan;  
OPT.wall_ind   = nan; 
OPT.wall_x     = nan;
OPT.wall_overwash = 0; % Wall overwash (21/5/2020) 
    % if wall_overwash = 1, sediment is allowed to pile behind the wall 
    % ... for accretion / positive sediment budget cases
OPT.wall_z_min_check = 0; % 1 = check to see if pt offshore of wall has eroded below WALL_Z_MIN
    % ... translation loop will break once that pt is reached
OPT.wall_z_min = -2; % Cap on how much erosion can occur offshore the wall
    % Once erosion immediately off the wall reaches WALL_Z_MIN (+ dS),...
    % ... these profiles will get tagged and be excluded from further analysis
OPT.wall_no_erode_behind = 0; % 0 = allows erosion behind wall; 1 = no erosion behind wall permitted

OPT.redist_ratio = 0.33; % if wallRedist = ON...
    % ... this is the cross-shore portion of the profile that the...
    % ... "behind the wall" volume is transfered to,...
    % ... starting offshore of the wallLevel (or wall_ind) point.

OPT.slumpCheck = 0; % (1=ON, 0=OFF0 runs across the final profile ...
    % ... to fix any vertical "sand cliffs", ...
    % ... this is required if prograding out from a vertical wall.
    % PROBLEM (19/2/2020) - SLUMPING can occur BEHIND WALLS
        % ... see DFIG3 script, and switch on duneSlump to investigate

OPT.deepOn     = 0; % Onshore tranport from <DoC1 (0=OFF [default], 1=ON)  
OPT.deepOn_V   = 0; % volume of deep onshore transport
                % Positive = Offshore; Negative = Onshore.
   
OPT.shortOutput = 1; % 0-FULL output (lots of memory, for debugging); 1 = SHORT Output (saves memory)
OPT.optimizer   = 1; % 1 = use OPTIMIZER function to find translation distance (FAST), 
    % 0 = translate profile across a fixed range

%% SETPROPERTY -> available with openEarthTools (Deltares)
[OPT, OPT.Set, OPT.Default] = setproperty(OPT, varargin); 
 
%% Translation increment (dX) and range (X_range)
% x is positive offshore
dx = x0(2) - x0(1); % little-dx = spacing of input profile
dX = OPT.dX;        % big-dX    = step-size to march profile on/offshore
X_dir = OPT.X_dir;
if X_dir == -1 % RECEEDING trend (default)
    X_range = [OPT.X_off: -dX : OPT.X_on]'; 
elseif X_dir == 1 % PROGRADING trend
    X_range = [OPT.X_on: dX : OPT.X_off]';  
else
    error(['X_dir not set correctly. X_dir must be -1 (receeding profile) or +1 (prograding profile),' ...
        ' dX sets the increment at which to calculate profile translation.']);
end
X_len = length(X_range);
OPT.X_len = X_len;
OPT.X_range = X_range;

% WARN IF X_on > 0
if OPT.X_on > 0 
    warning(['Onshore translation distance (X_on) has been set as a POSITIVE VALUE.'...
        10 'POSITIVE is OFFSHORE, therefore onshore translation must be a NEGATIVE VALUE']);
end

%% SLR settings
dS=OPT.dS; % sea-level rise magnitude (default = 1 m)
S_initial = OPT.S_initial;  % initial sea-level 
S_final   = S_initial + dS; % final sea-level
OPT.S_final = S_final;

%% DoC index (Depth of Closure)
% DoC - UPPER
DoC = OPT.DoC;
DoC_ind=find(z0>=DoC, 1, 'last') + 1  ; 
% DoC_ind = find(z0 < DoC, 1, 'first')   ; 
    % changed 13/5/2020, from 'first < DoC' to fix issue with full depletion 
    % in front of wall detecting a DoC immediately offshore of wall

% DoC2 - LOWER
if isempty(OPT.DoC2) 
    OPT.DoC2 = DoC - 1; %  round(2*DoC); % DEEPER DoC
    % 27/2/2020 - changed default DoC2 to just below DoC to avoid issues with prof not extending deep enough
end
DoC2     = OPT.DoC2;
DoC2_ind = find(z0<DoC2, 1, 'first');

if isempty(DoC2_ind)
    warning('WARNING: Profile does not reach specified DoC2 depth. Using min(z0) as DoC2.');
    [~, DoC2_ind] = min(z0);
       
end
OPT.DoC2_ind = DoC2_ind;

%% Rollover & Slope settings
rollover       = OPT.rollover;
backSlope      = OPT.roll_backSlope;
duneSlump      = OPT.duneSlump;
duneSlope      = OPT.duneSlope;
slumpCap       = OPT.slumpCap;
slumpCheck     = OPT.slumpCheck;

%% Toe / Crest Position (ToCr)
% In varargin, set OPT.toeCrest_level OR OPT.toeCrest_ind.
ToCr_level = OPT.toeCrest_level;
ToCr_ind   = OPT.toeCrest_ind;

% get ToCr_ind if ToCr_level used as input
if isempty(ToCr_ind) && ~isempty(ToCr_level)
    ToCr_ind = find(z0 >= ToCr_level, 1 , 'last');

% get ToCr_level if ToCr_ind used as input    
elseif ~isempty(ToCr_ind) && isempty(ToCr_level)       
    ToCr_level = z0(ToCr_ind);
    
% set ToCr_level and ToCr_ind if no inputs given   
elseif (isempty(ToCr_ind) && isempty(ToCr_level)) || isnan(ToCr_level)
	ToCr_level = max(z0);
    ToCr_ind = find(z0 == ToCr_level, 1, 'last');
    warning(['Both OPT.toeCrest_level and OPT.toeCrest_ind are empty.' ...
        ' Profile max elevation used as dune toe / barrier crest.' ]);
end

ToCr_level2 = ToCr_level + dS; % toe/crest level after SLR
% if rollover == 1
%     ToCr_level2 = ToCr_level + dS; % toe/crest level after SLR
% else 
%     ToCr_level2 = ToCr_level;
% end

OPT.ToCr_level2 = ToCr_level2;
ToCr_x     = x0(ToCr_ind); % x-location of ToCr
OPT.ToCr_x = ToCr_x;

%% Rock layer settings
rock   = OPT.rockSwitch;
z_rock = OPT.rockLayer;

% FIX NEGATIVE DUNE VOL ISSUES, IF INPUT PROF IS BELOW DUNE LEVEL
% Occurs when (z0 onshore of ToCr_ind is below z_rock)
% Amended (9/6/2020) -> offshore zone (>toeCrest_ind) allows hypothetical "bed below rocks"
    % ...            -> onshore zone  (<=toeCrest_ind) raises the initial bed to rock level (CHANGES Z0!!!)
if rock == 1 
    ind = find(x0 <= x0(ToCr_ind) & z0 < z_rock);
    z0(ind) = z_rock(ind);     % Raise ONSHORE Z0 to Z_ROCK ...
end

%% Dune Accretion / Horizontal Growth Settings
duneAccrete   = OPT.duneAccrete;   % use DUNE ACCRETION (1=ON, 0=OFF)
duneAcc_Vol   = OPT.duneAcc_Vol;   % VOLUME to accrete dunes
duneAcc_Dist  = OPT.duneAcc_Dist;  % DISTANCE (behind ToCr) over which to apply the DUNE VOL CHANGE
% disp(['Dune dist = ' num2str(duneAcc_Dist)])

duneHGrow        = OPT.duneHGrow;       % (1 = ON, 2 = OFF)
duneHGrow_m      = OPT.duneHGrow_m;     % horz distance (m) to grow dune
duneHGrow_crest  = OPT.duneHGrow_crest; % crest of dune (ht above present SL) to grow from

%% Deep Onshore Transport Settings
deepOn     =  OPT.deepOn; 
deepOn_V   =  OPT.deepOn_V;

%% Wall Redistribution (as per Beuzen2018)
% ---------- THIS CAN BE MOVED TO A SEPARATE FUNCTION !! -------
% For profiles with a wall or cliff, this sets the location of the
% wall. Potential erosion is calculated behind the wall 
% (using Atkinson2018 translation, assuming wall is not there),
% then the erosion volume is redistributed offshore of the wall.
wall       = OPT.wallSwitch;
wall_level = OPT.wall_level;
wall_ind   = OPT.wall_ind;
wall_x     = OPT.wall_x;

if wall == 1
    % if using wall_x, get wall_ind and wall_level
    if ~isnan(wall_x)
        [~,wall_ind]   = min(abs(x0 - wall_x));
        wall_level     = z0(wall_ind);
        OPT.wall_level  = wall_level;
        OPT.wall_ind    = wall_ind;
    % if wall_level is input, get wall_ind, then wall_x   
    elseif ~isnan(wall_level)
        wall_ind = find(z0<wall_level, 1, 'first') ;
        wall_x   = x0(wall_ind) ;
        OPT.wall_x = wall_x;
        OPT.wall_ind    = wall_ind;
    end
end

% Wall erosion limit
wall_z_min = OPT.wall_z_min  ;
if isnan(wall_ind) && wall == 1
    error('wall_ind is empty, check your wall settings!')
elseif wall == 1 && z0(wall_ind + 1) < wall_z_min
    OPT.z_min_initial = 1 ; % tag if INPUT PROF (z0) is already at MAX EROSION
else 
    OPT.z_min_initial = 0 ;
end

% Wall overwash (21/5/2020)
wall_overwash = OPT.wall_overwash; 
% if wall_overwash = 1, sediment is allowed to pile behind the wall 
    % ... for accretion / positive sediment budget cases

%% OLD WALL REDIST
% if wall == 1
%     % if wall_ind is input, get wallLevel
%     if isempty(wall_level) && ~isempty(wall_ind)
%         wall_level = z0(wall_ind);
%     
%     % if wallLevel is input, get wall_ind     
%     elseif ~isempty(wall_level) && isempty(wall_ind)
%         wall_ind = find(z0<wall_level, 1, 'first');
%     
%     % if both wallLevel and wall_ind are empty, use max height of profile
%     elseif isempty(wall_level) && isempty(wall_ind)
%         wall_level = max(z0);
%         wall_ind  = find(z0==wall_level, 1, 'last');
%         warning(['Both wall_ind and wall_level are empty.' ...
%             ' Profile max elevation used as wall location.' ]);
%         
%     end
% end

%% DUNE ACCRETION - if on, modify Z_RAISE and Z0
if duneAccrete==1
    % ACCRETE THE ACTIVE DUNES
    z0_temp = z0; % z0_temp is z0 modified by processes before the main loop,...
        % includes DUNE ACCRETION and DEEP ONSHORE TRANSPORT
    [~, ind_duneBack] = min(abs (x0 -  ( ToCr_x - duneAcc_Dist) ) ) ; % back point where the dune accretion starts
    dz_duneAcc = duneAcc_Vol / duneAcc_Dist;  % dz to accrete each dune cell
    z0_temp(ind_duneBack:ToCr_ind - 1) = z0(ind_duneBack:ToCr_ind - 1) + dz_duneAcc;
    
    % (DON'T!!) ERODE THE ACTIVE SHOREFACE -> allow translation to do the work (see notes below)
    
    active_len = x0(DoC_ind) - x0(ToCr_ind);
    dz_duneEro = - duneAcc_Vol / active_len;
%     z0_temp(ToCr_ind:DoC_ind) = z0(ToCr_ind:DoC_ind) + dz_duneEro;
    z0_temp(ToCr_ind:DoC_ind) = z0(ToCr_ind:DoC_ind);
        %  !!! (12/6/2020) leave z0_temp unaltered (don't add dune Erosion)...
        % ... instead rely on the algorithm to determine the correct (balancing) translation distance

    OPT.z_duneAcc       = z0_temp;
    OPT.dz_duneAcc      = dz_duneAcc;
    OPT.dz_duneEro      = dz_duneEro;
    
else
    z0_temp = z0;
end

%% DUNE HORIZONTAL GROWTH - if on, modify Z_RAISE and Z0
if duneHGrow == 1
    % HORIZ ACTIVE DUNES
    % includes DUNE ACCRETION, DUNE HORZ GROWTH, DEEP ONSHORE TRANSPORT
    
    ind1 = find(z0_temp >= S_initial + duneHGrow_crest , 1, 'last'); % dune crest index
    ind2 = find(x0 >= x0(ind1) + duneHGrow_m,  1, 'first'); % offshore pt to grow out to
    z0_temp(ind1:ind2) = S_initial + duneHGrow_crest ;
    
    % (DON'T!!) ERODE THE ACTIVE SHOREFACE -> allow translation to do the work (see notes below)

    OPT.duneHGrow_ind_crest       = ind1;  % dune crest index
    OPT.duneHGrow_ind_offPt       = ind2;  % offshore pt to grow out to
end

%% DEEP ONSHORE/OFFSHORE TRANSPORT (negative deepOn_V = onshore transport)
% see EMBf_VARB_TESTING_V2 for SINE CURVE AREA PROOFS

% SINE FUNCTION specifying WAVELENGTH (L) and AMPLITUDE (a)
% y = a .* sin(x.* (pi/ (L/2) ) ); 

% VOLUME OF HALF A SINE CURVE = AMPLITUDE (FULL WAVELENGTH/PI)
% V_(half wave sine) = a (L/pi)

% GIVEN A TARGET VOLUME, AND KNOWN HALF-WAVELENGTH, FIND AMPLITUDE
% a = pi*V / L

if deepOn == 1
	% ERODE(negative deepOn_V) / ACCRETE(positive) -> THE LOWER SHOREFACE
    lower_len = x0(DoC2_ind) - x0(DoC_ind+1);
%     dz_lowerEro = - deepOn_V / lower_len;
%     z0_temp(DoC_ind+1 : DoC2_ind) = z0_temp(DoC_ind+1 : DoC2_ind) + dz_lowerEro;

    bot_ind   = [DoC_ind+1 : DoC2_ind];
    x_bot     = x0(bot_ind);
    x_bot0    = x_bot - x_bot(1); % x-values for the bottom half of the profile, offset to start at 0
    a = (pi * -deepOn_V) / (2*lower_len); % amplitude
    dz_lowerEro = - a .* sin(x_bot0 .* (pi / (lower_len)));
    z0_temp(bot_ind) = z0_temp(bot_ind) + dz_lowerEro;

	% ACCRETE(neg)/ERODE(pos) THE UPPER SHOREFACE    
    active_len = x0(DoC_ind) - x0(ToCr_ind);
    dz_upperAcc =  -deepOn_V / active_len;
%     z0_temp(ToCr_ind:DoC_ind) = z0_temp(ToCr_ind:DoC_ind) + dz_upperAcc;
    z0_temp(ToCr_ind:DoC_ind) = z0_temp(ToCr_ind:DoC_ind);
    
    %  !!! (12/6/2020) leave z0_temp unaltered (don't add dune deepOn_V)...
    % ... instead rely on the algorithm to determine the correct (balancing) translation distance
        
end
OPT.z0_temp = z0_temp; % z0_temp = z0 + dune accretion + dune HORZ GROWTH + trend_x

%% Z0_ROCK (MAX OF UNION OF ROCK AND Z0) -> added for OFFSHORE REEFS (18/2/2020)
% use this for the final calcs, V0 = pre any changes (including the z0_temp changes above)

if rock == 1 
    z0_rock = max(z0, z_rock); % Take Z0_ROCK as max
elseif rock == 0
    z0_rock = z0; % added 9/6/2020, solved issue with NE layer causing issues, even when switched off
end
OPT.z0_rock = z0_rock;

%% Z0_ROCK -> Linear interp Z0 across ROCK OUTCROPS
if rock == 1
%     disp('interping across rocks')
    ind = (z0_rock <= z_rock);
    
    % INDEX OUTCROP SECTIONS
    ind_st = []; % last NON-OUTCROP section BEFORE an outcrop
    ind_en = []; % first NON-OUCROP section AFTER
    for i = 2:length(ind)-1
        if ind(i)==0 && ind(i+1) == 1
            ind_st = [ind_st; i];
        elseif ind(i)==0 && ind(i-1) == 1
            ind_en = [ind_en; i];
        end
    end    
    if length(ind_en) < length(ind_st) % add end index if outcrop extends...
        ind_en = [ind_en; ind_st(end)]; % ... to end of profile
    end
    
    % LINEAR INTERP ACROSS OUTCROP SECTIONS
    for i = 1:length(ind_st)
        x_temp = [x0(ind_st(i)) x0(ind_en(i)) ];
        z_temp = [z0(ind_st(i)) z0(ind_en(i)) ];
        ind = [ind_st(i):ind_en(i)];
        z0_temp(ind) = interp1(x_temp, z_temp, x0(ind));
    end
    clear x_temp z_temp
    
end

%% RAISE input profile by SLR (z0 -> z1)
z1 = z0_temp; % z1 = INTERIM1 -> SLR RAISED PROFILE
z1(ToCr_ind:DoC_ind-1) = z0_temp(ToCr_ind:DoC_ind-1) + dS;  % 31/8/2020 changed from DoC_ind to ind - 1
OPT.z_raise = z1;


%% UPDATE OPT STRUCTURE
OPT.S_final = S_final;
OPT.DoC_ind = DoC_ind;
OPT.toeCrest_level  = ToCr_level;
OPT.toeCrest_level2 = ToCr_level2;
OPT.toeCrest_ind    = ToCr_ind;
% OPT.wall_level      = wall_level;
% OPT.wall_ind        = wall_ind;

%% OPTIMISED PROFILE TRANSLATION 
% call OPTIMISER - finds x-translation distance with minimum dV_error
% outProf=[];


if OPT.optimizer == 1
    [Xi_winner, Xi_ind, outPROFS, OPT] = ST_OPTIMIZER(x0, z0, dV_input, OPT);
    OPT.Xi_winner = Xi_winner;
    OPT.Xi_ind    = Xi_ind;
    outProf = outPROFS(Xi_ind);
    OPT.n_guesses = length(outPROFS);

    if OPT.shortOutput == 0
        disp(['Final estimate:' 10 ...
            'Translation distance = ' num2str(Xi_winner) ' m' 10 ...
            'dV_error = ' num2str(outProf.dV_error, '%0.2f') ' m^3/m' 10 ...
            'n guesses = ' num2str(OPT.n_guesses)]);
    end
end

%% RANGE PROFILE TRANSLATION (NON-OPTIMISED)

if OPT.optimizer == 0
    outPROFS=[];
    for n = 1 : X_len % n = shift onshore (in metres as dx=1)

        Xi = X_range(n);
        [outProf1, OPT] = ST_TRANSLATOR(x0, z0, Xi, dV_input, OPT);
        outPROFS = [outPROFS; outProf1];
        
    end
    outProf = [];
end





%%


