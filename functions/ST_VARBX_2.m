%% ST_VARBX_2 --> NSW STYLE 
% 3/4/2021
% STORM CUT --> as per the method in Kinsela 2017
% Take beachface angle between [0 -> 2m], shift onshore
%
% INPUTS
% x0 = input x-values for profile
% x_temp = input z-values for profile
% vx.TC_level = Toe / Crest level ; e.g.  = 2.5;
% vx.dX = distance to shift beachface onshore (metres)
% vx.Z1 = 2;    % TOP of beachface
% vx.Z2 = 0;    % BOTTOM of BEAHCFACE (take angle from TOP to MID and shift onshore)
% vx.Z3 = -1;   % PIVOT pt (join translated MID to original PIVOT pt)
% vx.Z4 = -6;   % BOTTOM of STORM BAR pt (lost volume shifted as a sine curve offshore to 
% vx.duneSlope = slump slope for dune

function [vx] = ST_VARBX_2(x0, z0, vx, OPT)

% unpack
% TC_level = vx.TC_level;
TC_level = OPT.toeCrest_level;
dX = vx.dX; % dX --> positive offshore (usually a NEGATIVE value for storm cut)
duneSlope = OPT.duneSlope;

% Index ToeCrest and "vx" elevations
TC_ind = find(z0 >= TC_level, 1, 'last') + 1; % toe/crest
vx.i1 = find(z0 <= vx.Z1 & x0 >= x0(TC_ind), 1, 'first'); % Top of beachface
vx.i2 = find(z0 <= vx.Z2 & x0 >= x0(TC_ind), 1, 'first'); % Bottom of beachface
vx.i3 = find(z0 <= vx.Z3 & x0 >= x0(TC_ind), 1, 'first'); % Pivot pt
vx.i4 = find(z0 <= vx.Z4 & x0 >= x0(TC_ind), 1, 'first'); % Bottom of storm bar

% BEACHFACE ANGLE
Run1 = x0(vx.i2) - x0(vx.i1)   ;
Rise1 =  z0(vx.i1) - z0(vx.i2) ;
Slope1 = Rise1 / Run1;
Beta1 = atan2d(Rise1, Run1);
% Beta = atan2d( RISE[positive] / RUN(positive) )
% Beta_check = atan2d(1,10)

%% z1 = TRANSLATE BEACHFACE ONSHORE (distance dX)
z1 = z0;

% SHIFT BEACHFACE (vx.i1 : vx.2 -> onshore by dX)
z1(vx.i1 + dX : vx.i2 + dX) = interp1(x0([vx.i1,vx.i2]),[vx.Z1,vx.Z2], x0(vx.i1:vx.i2)  );

% JOIN BOTTOM OF BEACHFACE TO PIVOT POINT
z1(vx.i2 + dX : vx.i3) = interp1(x0([vx.i2 + dX, vx.i3]), [vx.Z2, vx.Z3], x0(vx.i2 + dX : vx.i3));

%% z2 = SLUMP DUNE
OPT = ST_OPT_defaults;
OPT.duneSlope = duneSlope;
[z2, ~, ~] = ST_SLUMP(x0, z1, z1, OPT);

%% z_final = FIND VOLUME -> add it back OFFSHORE (SINE CURVE)
V_initial = trapz(x0(1:vx.i3), z0(1:vx.i3) - vx.Z3); 
    % offset z to vx.Z3 = base of storm cut
V_final = trapz(x0(1:vx.i3), z2(1:vx.i3) - vx.Z3);
dV = V_final - V_initial;

% ERODE(negative deepOn_V) / ACCRETE(positive) -> THE LOWER SHOREFACE
z_temp = z2;
lower_len = x0(vx.i3) - x0(vx.i4);
%     dz_lowerEro = - deepOn_V / lower_len;
%     z_temp(DoC_ind+1 : DoC2_ind) = z_temp(DoC_ind+1 : DoC2_ind) + dz_lowerEro;

bot_ind   = [vx.i3 : vx.i4];
x_bot     = x0(bot_ind);
x_bot0    = x_bot - x_bot(1); % x-values for the bottom half of the profile, offset to start at 0
a = (pi * -dV) / (2*lower_len); % amplitude
dz_lowerEro =  a .* sin(x_bot0 .* (pi / (lower_len)));
z_temp(bot_ind) = z_temp(bot_ind) + dz_lowerEro;
z_final = z_temp;


%% PUT ALL IN "vx" STRUCTURE
vx.z0 = z0;
vx.z1 = z1;
vx.z2 = z2;
vx.z_final = z_final;
vx.BeachSlope = Slope1;
vx.Beta = Beta1;
vx.dV = dV; 






%%









%%












%%









%%