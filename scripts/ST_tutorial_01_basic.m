%% ST_tutorial_01_1
% May 2021
% Some examples of the basic function of SHORETRANS, 
% ... using "ENCROACHMENT" type translation (Type 4 in McCarroll et al., 2021)

close all, clear all, clc
ST_dir = 'D:\Dropbox\7_MODELS\013_ShrTrns'; % set local directory where ST is located
cd(fullfile(ST_dir, 'data'));
load('tutorial_01_data_x0z0.mat','x0','z0');

%% -------- TUT01-00: Input settings  ------------ %%
OPT = ST_OPT_defaults; % generate OPT (default settings)
   % ST_OPT_defaults contains a summary of  

% Change some settings
OPT.dS = 0.8 % dS = amount of sea level change (m)
OPT.DoC = -12; % (UPPER) Depth of closure (default = 10 m)
OPT.toeCrest_level = 2.5; % for dunes -> toe elevation (for low barriers -> crest elevation)

%% -------- TUT01-01: Run ST, 1 time step ------------ %%
[outProf1,~, OPT1] = ST_MAIN(x0, z0, 0, OPT);
z1 = outProf1.z_final; % translated profile

% Plot output
close all, figure, figpos, hold on;
title(['Profile translation, SLR = ' num2str(OPT.dS) ' m']);
plot(x0,z0, 'k');
plot(x0,z1, 'r:');
legend('Initial profile','Translated profile');
xlabel('Cross-shore distance (m)');
ylabel('Elevation (m)');

%% -------- TUT01-02: 1 time step, sediment budget deficit ------------ %%
dV_input = -200; % sediment budget (deficit) = -100 m3/m
[outProf2,~, OPT2] = ST_MAIN(x0, z0, dV_input, OPT);
z2 = outProf2.z_final; % translated profile

% Plot output
close all, figure, figpos, hold on;
title(['Profile translation']);
plot(x0,z0, 'k');
plot(x0,z1, 'r:');
plot(x0,z2, 'b-.');
legend('Initial profile','SLR = 0.8 m','... & budget = -100 m^3/m');
xlabel('Cross-shore distance (m)');
ylabel('Elevation (m)');


%% -------- TUT01-03: Add SEAWALL ------------ %%
clc
OPT = ST_OPT_defaults;
OPT.toeCrest_level = 2.5; % for dunes -> toe elevation (for low barriers -> crest elevation)
OPT.dS = 0.8;
OPT.duneSlope = 35;
OPT.wallSwitch = 1;
OPT.wall_x = 195;
[outProf3,~, OPT3] = ST_MAIN(x0, z0, 0, OPT);
z3 = outProf3.z_final;

% Plot output
close all, figure, figpos, hold on;
title(['Profile translation: Add a seawall']);
plot(x0,z0, 'k');
plot(x0,z1, 'r:');
plot(x0,z3, 'b-.');
plot([OPT.wall_x OPT.wall_x], [-2 8], 'k--','linewidth',2)
legend('Initial profile','SLR','SLR + WALL', 'Wall location');
xlabel('Cross-shore distance (m)');
ylabel('Elevation (m)');
xlim([50 600]);
ylim([-6 10]);

%% -------- TUT01-04: STORM DEMAND ------------ %%
vx.dV_target = -100;  % STORM DEMAND VOLUME (m3/m)
vx.Z1 = 2;           % TOP of beachface
vx.Z2 = 0;           % BOTTOM of BEACHFACE (take angle from TOP to MID and shift onshore)
vx.Z3 = -1;          % PIVOT pt (join translated MID to original PIVOT pt)
vx.Z4 = -6;          % BOTTOM of STORM BAR pt (lost volume shifted as a sine curve offshore to 

% Use output from TUT01-01 (SLR = 0.8 m) as input to storm demand function
[vx1] = ST_VARBX(x0, z1, vx, OPT1);
z4 = vx1.z_final;

% Plot output
close all, figure, figpos, hold on;
title(['Profile translation: SLR + STORM DEMAND']);
plot(x0,z0, 'k');
plot(x0,z1, 'b:');
plot(x0,z4, 'r-.');
legend('Initial profile','SLR','SLR + STORM');
xlabel('Cross-shore distance (m)');
ylabel('Elevation (m)');
xlim([50 600]);
ylim([-6 10]);




%%












