% Load data from file
data1 = load('Output3_CL200.txt');
time1 = data1(:, 1);
voltage1 = data1(:, 2);

figure()
plot(time1, voltage1);
xlim([1.995E5 2E5]);

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