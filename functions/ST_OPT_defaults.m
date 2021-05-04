%% function [OPT] = EMBf_profTrans_OPT_defaults
% RJM, 31/1/2020
% Get the OPT structure with default settings


function [OPT] = ST_OPT_defaults(x0,z0)

if nargin==0
    x0=nan;
    z0=nan;
end

%% SETTINGS (VARARGIN, use SETPROPERTY from openEarthTools)
% Translation increment/distance
    % x = positive offshore

OPT.dX    = 1;   % increment to translate profile off-/onshore (must be positive)
OPT.X_off = 0;    % max distance to translate profile OFFSHORE (progradation, default = 0)
OPT.X_on  = -100;  % max distance to translate profile ONSHORE (recession, default = -100)
% OPT.X_dir = -1;   % trend profile movement -> -1=Receeding; +1=Prograding

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
    
OPT.rollover=0; % use rollover (1 = ON, 0 = OFF, 2=ON (but don't keep up with SLR)
    % (ROLLOVER ON), barrier crest rolls back, maintaining height,...
    % ... and/or rising with sea-level. 
    % (ROLLOVER OFF), erosion will eat back into dune/barrier,...
    % ... crest height not mainteined or raised with SLR
    % (ROLLOVER = 2), barrier rolls back but is restricted to current absolute height,
    % ... will NOT maintain height relative to SLR
    % ROLLOVER=ON does not work with WALL=ON

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
    
% Rock layer options
OPT.rockSwitch=0; % rock switch (1 = ON, 0 = OFF [default])
    % RockLayer ON -> a check is made to ensure the translated profile,...
    % ... stays above the rock layer (see next input, OPT.Rock_z)
    
OPT.rockLayer = nan .* z0; % change 10/6/2020, change from -100 to nans
    % rockLayer =  z0 - 100; % CHANGE 31/1/2020, create DUMMY ROCK LAYER
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

OPT.slumpCheck = 0; % (1 = ON, 0 = OFF0 runs across the final profile ...
    % ... to fix any vertical "sand cliffs", ...
    % ... this is required if prograding out from a vertical wall.
    % PROBLEM (19/2/2020) - SLUMPING can occur BEHIND WALLS
        % ... see DFIG3 script, and switch on duneSlump to investigate

OPT.deepOn     = 0; % Onshore tranport from <DoC1 (0=OFF [default], 1=ON)        
OPT.deepOn_V   = 100; % volume of deep onshore transport
   
OPT.shortOutput = 1; % 0-FULL output (lots of memory, for debugging); 1 = SHORT Output (saves memory)
OPT.optimizer   = 1; % 1 = use OPTIMIZER function to find translation distance (FAST), 
    % 0 = translate profile across a fixed range



















%%