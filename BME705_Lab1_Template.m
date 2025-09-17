%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BME705: Rehabilitation Engineering
% Lab 1: Applications of FES and EMG in Rehabilitation
%
% Created by: Devon Santillo, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Names: Tanvir Hassan, Ilina Petrova
 
% Student IDs: 501104056
%

close all
clear all

%% Initializations
% Load all relevent files and, if necessary, save the components of them
% into variables

%Fes-10
%Ta-1

TA = Data = ("//path") %TA data
load (fatigue data);
load(stimulation data)

%define the three sampling frequencies for each data set



%% Part 1: Introduction to EMG analysis
% separate EMG and force data using the dot operator
EMG = (Data.increase_ta_emg);
Force = (Data.increase_ta_force);

% Define time as a vector


% Processing: 

%1)original graphs


% 2) normalizations


%3) rectification; 


%4) Filtered EMG signal
fc = % corner frequency
%establishing the transfer function of a 4th order butterworth filter


%filtering emg and graphing output



%5) calculating iEMG



%6) Dividing signal into individual contractions
% using the reshape function to evenly split each matrix into even columns
% each column should be a contraction event
EMGr = reshape();
Forcer = reshape();

%%% This next code is optional to compute the envelopes automatically %%%
%establishing blank arrays to input RMS and median frequency data

%Arrays 1 and 2 are established as sets of zeros equivalent to the # of
%envelopes

%using a for loop to go through every column of the reshaped EMG matrix

for i = %i variable should go from 1 to your the total number of envelopes
    
    %Calculating i-th RMS variable corresponding to i-th column 
    %in reshaped EMG matrix
    % Array1(i) = rms of EMGr(i)
   
    %plotting the EMG signal with RMS for the i-th column
    %recommended to use 'num2str()' function for a dynamic plot title
    
    
    %Calculating the corresponding average force for the i-th contraction
    % Array2(i) = rms of EMGr(i)
end

%plot emg vs. force magnitude relationship

%% Part 2: Effects of stimulation frequency on muscle fatigue

% Plot the raw data
%30Hz raw data

%60Hz raw data

%
 


 
% Time synchronizing - 
% This step has been preformed in advance for you by your TAs

 
% Synchronize averaging
 

 
% New means and Standard Deviations of synchronized signals


% Define each signal time, as they are sampled over different durations

 
% Plot: Time-synchronized forces and Average force for: 30Hz and 60Hz and
% Use holds to plot both signals on the same chart

 

%% Part 3: Stimulation input - motor threshold investigation

% Define time of signal

% Plot original force signal

 
% Zero-mA Force: Subtract the dc-offset at 0mA stimulation


% Plot the force and stimulation curves on the same chart. 

