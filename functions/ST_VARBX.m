%% ST_VARBX --> Calculate storm demand impact on profile
% 5/4/2021
% Follows the method in Kinsela et al. (2017)
    % Kinsela, M. A., Morris, B. D., Linklater, M., & Hanslow, D. J. (2017). 
    % Second-pass assessment of potential exposure to shoreline change in New South Wales,
    % Australia, using a sediment compartments framework. Journal of Marine Science and Engineering, 5(4), 61.
    % https://www.mdpi.com/2077-1312/5/4/61
%
% INPUTS
% x0 = input x-values for profile
% x0 = input z-values for profile
% vx.dV_target  % STORM DEMAND VOLUME (m3/m)
% vx.Z1 = 2;    % TOP of beachface
% vx.Z2 = 0;    % BOTTOM of BEACHFACE (take angle from TOP to MID and shift onshore)
% vx.Z3 = -1;   % PIVOT pt (join translated MID to original PIVOT pt)
% vx.Z4 = -6;   % BOTTOM of STORM BAR pt (lost volume shifted as a sine curve offshore to 
% OPT
% OPT.duneSlope = slump slope for dune
% OPT.toeCrest_level

function [vx] = ST_VARBX(x0, z0, vx, OPT)

dV = 9999;
dV_target = vx.dV_target;
vx.dX = 0;

%% IF NO WALL PRESENT
if OPT.wallSwitch == 0
    while dV > dV_target

        vx.dX = vx.dX - 1;
        vx = ST_VARBX_2(x0, z0, vx, OPT);
        dV = vx.dV;
        disp(['VARBX: For dX = ' num2str(vx.dX) ' m,  dV = ' ...
            num2str(vx.dV, '%0.1f') ' m3/m']);
    end
    vx.dV = dV;
end


%% IF WALL PRESENT
if OPT.wallSwitch == 1
    while dV > dV_target
        vx.dX = vx.dX - 1;
%       vx.dX = -21

        vx = ST_VARBX_2(x0, z0, vx, OPT);
        dV_noWall = vx.dV;
        
        vx.z_noWall = vx.z_final;
        z_temp = vx.z_final; 
        z_temp(x0<x0(OPT.wall_ind) ) = z0(x0 < x0(OPT.wall_ind) );
        
        [vx.z_final, ~] = ST_WALL_VOL...
            (x0, z0, z_temp, vx.z_noWall, OPT, vx.dX);
        
        V_initial = trapz(x0(1:vx.i3), z0(1:vx.i3) - vx.Z3)  ;
        % offset z to vx.Z3 = base of storm cut
        V_final = trapz(x0(1:vx.i3), vx.z_final(1:vx.i3) - vx.Z3);
        dV = V_final - V_initial  ;
        
        disp(['VARBX: For dX = ' num2str(vx.dX) ' m,  dV = ' ...
            num2str(dV, '%0.1f') ' m3/m']);
    end
     vx.dV = dV;
end


%%







%%









%%