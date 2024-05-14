% Load data from file
data1 = load('Output3_CL200.txt');
time1 = data1(:, 1);
voltage1 = data1(:, 2);

data2 = load('Output3_CL1000.txt');
time2 = data2(:, 1);
voltage2 = data2(:, 2);

data3 = load('Output3_CL2000.txt');
time3 = data3(:, 1);
voltage3 = data3(:, 2);

% Plot the action potential waveform for data1
figure(1)
plot(time1, voltage1); 
xlabel('Time (ms)');
ylabel('Vm (mV)'); 
xlim([1.995E5 2E5]);
title('Action Potential CL200');
set(gca, 'fontsize', 14);
set(gcf, 'color', 'w');
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesFontSize', 20);

% APD90 Calculation for data1 (200CL)
% Enable zoom
zoom on;

% Wait for the user to finish zooming (can add a pause here if desired)
% disp('Zoom in on the plot and press Enter when done...');
% pause;

% Turn off zoom
zoom off;

% Get two points from the user
disp('Select two points on the plot:');
[x, y] = ginput(2);

% Calculate the time difference between the two selected points
time_diff = abs(x(2) - x(1));

%Calculate 90% of time difference 
depolarization90 = 0.9*time_diff; 

% Display the result
fprintf('Time difference between the selected points: %f seconds\n', depolarization90);


% Plot the action potential waveform for data2
figure(2)
plot(time2, voltage2); 
xlabel('Time (ms)');
ylabel('Vm (mV)'); 
title('Action Potential CL1000');
set(gca, 'fontsize', 14);
set(gcf, 'color', 'w');
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesFontSize', 20);

% Find Peaks for data2
maxpeakheight2 = 30; % Adjust as needed
[peaks2, peak_locs2] = findpeaks(voltage2, 'MinPeakHeight', maxpeakheight2);

% APD90 Calculation for data2
tolerance2 = 0.5; % Adjusted!
apd90_values2 = zeros(length(peaks2), 1); % Initialize array to store APD90 values
for i = 1:length(peaks2)  % Loop through each peak index
  peak_voltage = voltage2(peak_locs2(i)); % Corrected indexing

  % Calculate 90% repolarization level relative to the peak
  repolarization_level = peak_voltage - 0.9 * (peak_voltage - min(voltage2(peak_locs2(i):end)));

  % Define tolerance - should we do this as a percentage of repolarization
  % level?
  lower_bound = repolarization_level - tolerance2;
  upper_bound = repolarization_level + tolerance2;

  % Find indices where voltage falls within the specified range (starting from peak)
  start_indices = find(voltage2(peak_locs2(i):end) >= lower_bound & voltage2(peak_locs2(i):end) <= upper_bound);
  start_indices = start_indices + peak_locs2(i) - 1;  % Adjust indices based on starting point

  % Calculate APD90 for this peak
  apd90_values2(i) = time2(start_indices(1)) - time2(peak_locs2(i)); % Relative to peak

  % Display APD90 for this peak
  disp(['APD90 for peak ' num2str(i) ': ' num2str(apd90_values2(i)) ' ms']);
end

% Plot the action potential waveform for data3
figure(3)
plot(time3, voltage3); 
xlabel('Time (ms)');
ylabel('Vm (mV)'); 
title('Action Potential CL2000');
set(gca, 'fontsize', 14);
set(gcf, 'color', 'w');
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesFontSize', 20);

% Find Peaks for data3
maxpeakheight3 = 35; % Adjust as needed
[peaks3, peak_locs3] = findpeaks(voltage3, 'MinPeakHeight', maxpeakheight3);

% APD90 Calculation for data3
tolerance3 = 0.5; % Adjusted!
apd90_values3 = zeros(length(peaks3), 1); % Initialize array to store APD90 values
for i = 1:length(peaks3)  % Loop through each peak index
  peak_voltage = voltage3(peak_locs3(i)); % Corrected indexing

  % Calculate 90% repolarization level relative to the peak
  repolarization_level = peak_voltage - 0.9 * (peak_voltage - min(voltage3(peak_locs3(i):end)));

  % Define tolerance 
  lower_bound = repolarization_level - tolerance3;
  upper_bound = repolarization_level + tolerance3;

  % Find indices where voltage falls within the specified range (starting from peak)
  start_indices = find(voltage3(peak_locs3(i):end) >= lower_bound & voltage3(peak_locs3(i):end) <= upper_bound);
  start_indices = start_indices + peak_locs3(i) - 1;  % Adjust indices based on starting point

  % Calculate APD90 for this peak
  apd90_values3(i) = time3(start_indices(1)) - time3(peak_locs3(i)); % Relative to peak

  % Display APD90 for this peak
  disp(['APD90 for peak ' num2str(i) ': ' num2str(apd90_values3(i)) ' ms']);
end

% Extract the last APD90 value for each dataset
last_apd90_data1 = depolarization90;
last_apd90_data2 = apd90_values2(end);
last_apd90_data3 = apd90_values3(end);

% Define cycle lengths
cycle_lengths = [200, 1000, 2000];

% Create a figure
figure(4);
plot(cycle_lengths, [last_apd90_data1, last_apd90_data2, last_apd90_data3], 'o', 'MarkerSize', 10);
xlabel('Cycle Length (ms)');
ylabel('APD90 (ms)'); 
set(gca, 'fontsize', 14);
set(gcf, 'color', 'w');
set(groot, 'defaultLineLineWidth', 2);
xlim([0, 2000]);
ylim([0, 700]);

% Fit a curved line (polynomial) to the data
degree = 2; % Degree of the polynomial (quadratic in this case)
p = polyfit(cycle_lengths, [last_apd90_data1, last_apd90_data2, last_apd90_data3], degree);

% Generate points for the fitted curve
x_fit = linspace(0, 2000, 100); % Adjust range and number of points according to your preference
y_fit = polyval(p, x_fit);

% Plot the fitted curve
hold on;
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
legend('Data', 'Fitted Curve', 'Location', 'best');
hold off;
