%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BME705: Rehabilitation Engineering
% Lab 2: Balance Board Signal Analysis
%
% Created by: Dharmendra Gurve, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Names:
%
% Student IDs: SOLUTION
%
%%

clear all
close all

%%%%%%%%%%%%%%% Part 1 : %%%%%%%%%%%%%%%
%% Load the reference accelerometer data for part 1

ref_data = xlsread( "Reference_Part1.xlsx" );

% Define sampling frequency 
fs = 30;


% Separate the sample_index from the excel file (column 1) to create time axis

% Define time
dlen = length(ref_data)
t= [0:dlen-1]/fs;

% Separate the X-axis and Y- axis from the excel file (column 2 and column 3)
% colum 1 is sample number
ref_X = ref_data(:, 2);
ref_Y = ref_data(:, 3);


%Plot the X and Y axis data against time as separate lines on the same grid
%plotting the x and y samples against t
plot(t, ref_X, 'b', t, ref_Y, 'r');
xlabel('Time (s)');
ylabel('Acceleration (g)');
title('Reference Accelerometer Data');
legend('X-axis', 'Y-axis');
grid on;


%% Start the signal from the actual activity 
% 35 sec * 30 Hz = 1050 (+ start point) e.g. start point for referance file= 75
r1_start_point = 75; % Starting point for Subject 1
r1_X_aligned = ref_X(r1_start_point:end); % Aligning X-axis data
r1_Y_aligned = ref_Y(r1_start_point:end); % Aligning Y-axis data


% Define time
dlen = length(r1_X_aligned);
time_aligned = [0:dlen-1]/fs;


% plot the newly aligned reference X and Y axis data against time aligned
% Plot the aligned reference X and Y axis data against time aligned
figure;
plot(time_aligned, r1_X_aligned, 'b', time_aligned, r1_Y_aligned, 'r');
xlabel('Time (s)');
ylabel('Acceleration (g)');
title('Aligned Reference Accelerometer Data');
legend('Aligned X-axis', 'Aligned Y-axis');
grid on;




%%%%%%%%%% Measuring Signal Similarities using Cross-correlation between signals ref_X to sub_X, ref_Y to sub_Y 

%% Cross-correlation for Sub1 

sub1_part1 = xlsread( "Subject1_Part1.xlsx" );
sub2_part1 = xlsread("Subject2_Part1.xlsx");
sub3_part1 = xlsread("Subject3_Part1.xlsx");

%Subject 1 data 
start_time_seconds = .1;% Change this according to the start time of the sample
s1_start_point = (30 * start_time_seconds);
s1_X_aligned = sub1_part1(s1_start_point:end, 2); % Aligning X-axis data for Subject 1
s1_Y_aligned = sub1_part1(s1_start_point:end, 3); % Aligning Y-axis data for Subject 1
dlen = length(s1_X_aligned);
time_s11 = [0:dlen-1]/fs;

%Subject 2 data
start_time_seconds2 = 0.1;
s2_start_point = 30 * start_time_seconds2;
s2_X_aligned = sub2_part1(s2_start_point:end, 2); % Aligning X-axis data for Subject 2
s2_Y_aligned = sub2_part1(s2_start_point:end, 3); % Aligning Y-axis data for Subject 2
dlen = length(s2_X_aligned);
time_s21 = [0:dlen-1]/fs;

%subject 3 data
start_time_seconds3 = 0.1;
s3_start_point = 30 * start_time_seconds3; % Aligning X-axis data for Subject 3
s3_X_aligned = sub3_part1(s3_start_point:end, 2); % Aligning X-axis data for Subject 3
s3_Y_aligned = sub3_part1(s3_start_point:end, 3); % Aligning Y-axis data for Subject 3
dlen = length(s3_X_aligned);
time_s31 = [0:dlen-1]/fs;

% Plotting Subject data in a loop
subjects = {s1_X_aligned, s1_Y_aligned; s2_X_aligned, s2_Y_aligned; s3_X_aligned, s3_Y_aligned};
time_data = {time_s11, time_s21, time_s31};
titles = {'Subject 1 Accelerometer Data', 'Subject 2 Accelerometer Data', 'Subject 3 Accelerometer Data'};

figure;
for i = 1:3
    subplot(3, 1, i);
    plot(time_data{i}, subjects{i, 1}, 'b', time_data{i}, subjects{i, 2}, 'r');
    xlabel('Time (s)');
    ylabel('Acceleration (g)');
    title(titles{i});
    legend('X-axis', 'Y-axis');
    grid on;
end


% Find the cross-correlation between Subject 1 and reference signal


% Calculate the Correlation coefficients


% Plot the correlation overlap of x-axis    


% Plot the correlation overlap of y-axis 



% repeat the above process for all the subjects, either using repeat coding
% or a logic loop (for, while, etc.)

%%%%%%%%%%%%%%% Part 2 : %%%%%%%%%%%%%%%

%% load the data from all subjects

ref_data=xlsread('Reference_Part2.xlsx');

sub1_part2=xlsread(%assigned dataset 1);

sub2_part2=xlsread(%assigned dataset 2);

sub3_part2=xlsread(%assigned dataset 3);


% Separate the  X-axis and Y- axis from the excel file (column 2 and column 3)





%% Calculate root mean square error between Ref file and subject files (X-X, Y-Y)

 