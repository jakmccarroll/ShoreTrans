%% ST_tutorial_02_1
% May 2021
% SHORETRANS translation types:
    % TYPE 1 - ROLLOVER (onshore transport)
    % TYPE 2 - CREST ROLLOVER, CREST HEIGHT KEEP UP WITH SLR
    % TYPE 3 - CREST ROLLOVER, MAINTAIN INITIAL CREST HEIGHT
    % TYPE 4 - ENCROACHMENT (offshore transport)
% ... using "ENCROACHMENT" type translation (Type 4 in McCarroll et al., 2021)

close all, clear all, clc
ST_dir = 'D:\Dropbox\7_MODELS\013_ShrTrns'; % set local directory where ST is located
cd(fullfile(ST_dir, 'data'));
load('idealised_profiles.mat');

% SELECT PROFILE (BARRIER)
x0 = PROF3(12).xi;
z0 = PROF3(12).zi;


%% SETTINGS
% duneToe = 3;
duneToe = max(z0);
DoC     = -10;
dS      = 1; % SLR
lagoonToe_z = -1.5
lagoonToe_ind = find(z0 > lagoonToe_z, 1, 'first')


%% ------------ RUN SHORETRANS -------------- %%
% TYPE 1 - FULL ROLLOVER (3) (use lagoon toe as toeCrest_ind)
[out(1).outProf, out(1).OPT]= ST_MAIN(x0, z0, 0, ...
    'toeCrest_ind', lagoonToe_ind , 'DoC',DoC, 'dS', dS,...
    'rollover', 1, 'roll_backSlope',4);
z1 = out(1).outProf.z_final;

% TYPE 2 - (crest keep up with SLR)
[out(2,1).outProf, out(2,1).OPT]= ST_MAIN(x0, z0, 0,...
    'toeCrest_level', duneToe , 'DoC',DoC, 'dS', dS,...
    'rollover', 1, 'roll_backSlope',4);
z2 = out(2).outProf.z_final;

% TYPE 3 - (restict to initial height, don't keep up with SLR)
[out(3,1).outProf, out(3,1).OPT]= ST_MAIN(x0, z0, 0,...
    'toeCrest_level', duneToe , 'DoC',DoC, 'dS', dS,...
    'rollover', 2, 'roll_backSlope',4);
z3 = out(3).outProf.z_final;

% TYPE 4 - ENCROACHMENT / NO ROLLOVER (encroach/erode barrier)
[out(4,1).outProf, out(4,1).OPT]= ST_MAIN(x0, z0, 0,...
    'toeCrest_level', duneToe , 'DoC',DoC, 'dS', dS,...
    'rollover', 0);
z4 = out(4).outProf.z_final;

%% PLOT OUTPUTS
close all, figure, figpos, hold on;
plot(x0,z0, 'k')
plot(x0,z1, ':b');
plot(x0,z2, ':', 'color', [0 .5 0]);
plot(x0,z3, ':', 'color', [1 .5 0]);
plot(x0,z4, ':r');
legend('Initial','Type 1: Rollover', 'Type 2: Crest keep up',...
    'Type 3: Crest restrict to initial height', 'Type 4: Encroach');














