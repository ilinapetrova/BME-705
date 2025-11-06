%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BME705: Rehabilitation Engineering
% Lab 3: Center of Pressure (COP) and Sway Analysis 
%
% Created by:  Devon Santillo and Dharmendra Gurve, 2020
% Updated by: Dan Genkin, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name:
% Student ID:
%
%
%%
clear; clc; close all;

% Load the data and define the sampling frequency

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyzeFP function: "Given" (do not edit): convert raw forceplate data
% INPUT = 16 channels of raw forceplate data;
% OUTPUT: xCOPL, yCOPL, xCOPR, yCOPR, Fz_R, Fz_L
% Note: Rv_R and Rv_L are vertical forces of right and left forceplate 

% You can analyze the eyes open (EO) and eyes closed (EC) data seperately in a repeated process or in a logic loop
% but remember you cant have matrices as the iterative variable in a loop
%% 1.1) Calcualte and Plot the COPnet of the eyes open data
% Apply the provided function to obtain necessary outputs for further analysis
[xCOPL, yCOPL, xCOPR, yCOPR, Rv_R, Fz_L] = analyzeFP();

% Shift all X-coordinates to the middle of the split force plate by +/- 12.5825cm

% Filter all data
% Use of a Low-pass filter: a Butterworth is recommended
% define corner frequency and normalized cut-off frequency
fc=10;
Wn = 
[b, a] = 

% Calculate COPnet for overall system using Equation A from the lab manual (Winter et al., 2003 [1]): 

% Plot net COP (the stabilograms)

%% 1.2) Calcualte and plot the COPnet of the eyes closed data

%% 1.3) Find total excursions and mean velocity using equations B-I
% Using the method described in Prieto el al., 1994 paper [2]
% Select 20 seconds of the data in the time frame from 40-60 seconds
T=20;
start_point = fs*40;
end_point = fs*60;

% Downsample A/P and M/L axes to 100 Hz
fs_new = 100;
downsample();


% Zero-mean AP and ML : Equations B & C (Prieto el al., 1994)

% Total Excursion: Equations D-F (Prieto el al., 1994)
% Calculate TOTEX ML, AP, Net for Open and Closed

TOTEXap = 0;
TOTEXml = 0;
TOTEX = 0;

for i = 1:length(x)-1
    TOTEXap = TOTEXap + ...
    TOTEXml = TOTEXml + ...
    TOTEX = TOTEX + ...
end

% Mean velocity: Equations G-I (Prieto el al., 1994)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2.1) Processing EMG data: Consider the timing of activation of the Sol. and TA muscles with respect to the COP sway in the AP direction.

% Define time vectors 

% Rectification 

% Normalization

% Signal Envelope Filtering


%% 2.2) Processing sway data (ML, AP) - Repeat force plate analysis
% This section should follow the first part of the analysis
% call the analyzeFP function again (do not edit it)

% Apply the provided function to obtain necessary outputs for further analysis

% Shift X-coordinates to the middle of the split force plate by +/- 12.5825cm

% Filter all the slow & fast, L & R COP data

% Calculate COPnet for overall system of both fast and slow sway using Equation A from the lab manual

% Plot Open and Closed COPnet Stabilograms 
figure;
subplot(2,1,1);

subplot(2,1,2);

%% 2.3) Plot the slow and fast EMG and the Sway data together

% slow
figure;
subplot(2,1,1);

subplot(2,1,2);

% fast
figure;
subplot(2,1,1);

subplot(2,1,2);
