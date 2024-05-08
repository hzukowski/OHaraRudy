% Load data from file
data = load('results_project3.txt');
time = data(:, 1);
voltage = data(:, 2);

%Find Peaks
[peaks, peak_locs] = findpeaks(voltage, 'MinPeakHeight', 30); % Using MATLAB's built-in function, looking for above a certain y value

% Plot the action potential waveform 
figure
plot(time, voltage); 
xlabel('Time (ms)');
ylabel('Vm (mV)'); 
title('Action Potential, Ko=5.4 mV');
set(gca, 'fontsize', 14);
set(gcf, 'color', 'w');
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesFontSize', 20);

% APD90 Calculation
for i = 1:length(peaks)  % Loop through each peak index
  peak_voltage = voltage(peak_locs(i));

  % Calculate 90% repolarization level relative to the peak
  repolarization_level = peak_voltage - 0.9 * (peak_voltage - min(voltage(peak_locs(i):end)));

  % Define tolerance 
  tolerance = 0.05;
  lower_bound = repolarization_level - tolerance;
  upper_bound = repolarization_level + tolerance;

  % Find indices where voltage falls within the specified range (starting from peak)
  start_indices = find(voltage(peak_locs(i):end) >= lower_bound & voltage(peak_locs(i):end) <= upper_bound);
  start_indices = start_indices + peak_locs(i) - 1;  % Adjust indices based on starting point

  % Calculate APD90 for this peak
  apd90_values(i) = time(start_indices(1)) - time(peak_locs(i)); % Relative to peak

  % Display APD90 for this peak
  disp(['APD90 for peak ' num2str(i) ': ' num2str(apd90_values(i)) ' ms']);
end

cyclelength = [200, 1000, 2000];
apd90 = [200, 265, 270];

plot(cyclelength, apd90); 
