clear all;
close all;
clc;

%% ===================================================== LOAD SIGNALS =====================================================
% LOAD THE SIGNALS
% EMG sample values are in mV
% First column is time; second column is force; third column is EMG signal in mV
emg = load('EMGforce1.txt');
time = emg(:, 1); % Extract the time from the loaded data
force = emg(:, 2); % Extract the force signal from the loaded data 
emg_sig = emg(:, 3); % Extract the EMG signal from the loaded data

%% ===================================================== FORCE NORMALIZATION =====================================================
% NORMALIZE THE FORCE SIGNAL such that the minimum value is 0 and the maximum value is 100
normalized_force = normalize(force, 'range', [0, 100]);

fs = 2000; % Sampling rate
N = length(force); % Number of samples
f = 1:N;
f = fs/N*f; % Frequency range in Hz, corresponding to real time frequency

% Plot the EMG signal and the normalized force signal
figure(1);
subplot(3, 1, 1);
plot(time, emg_sig, 'b'); % EMG signal in blue
xlabel('Time (s)'); ylabel('EMG (mV)'); title('Original EMG Signal');
subplot(3, 1, 2);
plot(time, normalized_force, 'r'); % Normalized force signal in red
xlabel('Time (s)'); ylabel('Normalized Force (%MVC)'); title('Normalized Force Signal');

%% ===================================================== EMG SIGNAL FILTERING =====================================================
% FOURIER TRANSFORM OF THE EMG SIGNAL
emg_sig_fft = fft(emg_sig);

% FILTER THE EMG SIGNAL using a stacked high-pass and low-pass filter
fc_lp = 55; % Cutoff frequency of Signal 1
% fc_lp = 60; % Cutoff frequency of Signal 2
% fc_lp = 25; % Cutoff frequency of Signal 3
% fc_lp = 40; % Cutoff frequency of Signal 4
Wn_lp = 2*pi*fc_lp; % Cutoff frequency in rad/s

fc_hp = 10; % Cutoff frequency of Signal 1
% fc_hp = 4; % Cutoff frequency of Signal 2
% fc_hp = 40; % Cutoff frequency of Signal 3
% fc_hp = 20; % Cutoff frequency of Signal 4
Wn_hp = 2*pi*fc_hp; % Cutoff frequency in rad/s

% LOW-PASS FILTER DESIGN
[b, a] = butter(4, Wn_lp, 'low', 's'); % Butterworth filter returns in continuous time
lowpass_ct = tf(b, a);
[mag_lp, phase_lp] = bode(lowpass_ct, 2*pi*f);
mag_lp = squeeze(mag_lp);
phase_lp = squeeze(phase_lp);
H_freq_lp = mag_lp.*exp(1i*phase_lp/180*pi); % Convert the magnitude and phase to complex number
Y_lp = H_freq_lp.*emg_sig_fft; % Calculate the expected spectrum of the output signal
H_lp_dt = c2d(lowpass_ct,1/fs);

% HIGH-PASS FILTER DESIGN
[b, a] = butter(4, Wn_hp, 'high', 's'); % Butterworth filter returns in continuous time
highpass_ct = tf(b, a);
[mag_hp, phase_hp] = bode(highpass_ct, 2*pi*f);
mag_hp = squeeze(mag_hp);
phase_hp = squeeze(phase_hp);
H_freq_hp = mag_hp.*exp(1i*phase_hp/180*pi); % convert the magnitude and phase to complex number
Y_hp = H_freq_hp.*emg_sig_fft; % calculate the expected spectrum of the output signal
H_hp_dt = c2d(highpass_ct,1/fs);

% Apply the stacked high-pass and low-pass filter to the signal (in the freq domain)
% Apply the high-pass filter 
emg_fft_filtered_hp = emg_sig_fft.*squeeze(freqresp(H_hp_dt,2*pi*f));
% Apply the low-pass filter
emg_fft_filtered_hplp = emg_fft_filtered_hp.*squeeze(freqresp(H_lp_dt,2*pi*f));
% Take the inverse Fourier transform to get the filtered signal
emg_filtered_hplp = ifft(emg_fft_filtered_hplp);

% Calculate the combined magnitude of the high-pass and low-pass filters
combined_mag = 20*log10(mag_hp) + 20*log10(mag_lp);

% Plot the frequency spectrum of the EMG signals
figure(2);
subplot(3,1,1)
semilogx(f(1:N/2), 20*log10(abs(emg_sig_fft(1:N/2))));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Spectrum of (Original) EMG signal');
subplot(3,1,2)
semilogx(f(1:N/2), combined_mag(1:N/2));
title('Combined Frequency Response of Highpass and Lowpass Filters');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
subplot(3,1,3)
semilogx(f(1:N/2),20*log10(abs(emg_fft_filtered_hplp(1:N/2))))
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Frequency Spectrum of HP then LP Filtered EMG Signal')

% Plot the filtered EMG signal  
figure(1)
subplot(3, 1, 3);
plot(time, emg_filtered_hplp, 'b')
xlabel('Time (s)'); ylabel('EMG (mV)'); title('Filtered EMG Signal');

%% ===================================================== SIGNAL SEGMENTATION =====================================================
% AUTOMATIC SEGMENTATION WITH A PEAK DETECTION APPROACH
% Find local maxima of the normalized force signal
[pks, locs_peaks] = findpeaks(normalized_force); % Stores amplitude and location of peaks
% Find local minima of the normalized force signal
[trough, locs_troughs] = findpeaks(-normalized_force); % Stores amplitude and location of troughs
trough = -trough; % Change the sign of troughs for easier manipulation on the coming steps

% Set a threshold to filter out the unnecessary local maxima and minima of the signal
% Create a new arrays to store the peaks and their corresponding locations that are bigger than the threshold
new_peaks = zeros(length(locs_peaks), 1);
new_locs_peaks = zeros(length(locs_peaks), 1);
count_peaks = 0;

% If peaks are bigger then 10, it is considered a local maxima. 
for i = 1:length(locs_peaks)
    if pks(i) > 10
        count_peaks = count_peaks + 1;
        new_peaks(count_peaks) = pks(i);
        new_locs_peaks(count_peaks) = locs_peaks(i);
    end
end

new_peaks = new_peaks(1:count_peaks);
new_locs_peaks = new_locs_peaks(1:count_peaks);

% Create a new array to store the troughs and their corresponding locations that are bigger than the threshold
new_trough = zeros(length(locs_troughs), 1);
new_locs_trough = zeros(length(locs_troughs), 1);
count_trough = 0;

% If troughs are less then 10, it is considered a local minima. 
for i = 1:length(locs_troughs)
    if trough(i) < 10 
        count_trough = count_trough + 1;
        new_trough(count_trough) = trough(i);
        new_locs_trough(count_trough) = locs_troughs(i);
    end
end

new_trough = new_trough(1:count_trough);
new_locs_trough = new_locs_trough(1:count_trough);

% Plot the new peaks and troughs
figure(3);
plot(time, normalized_force, 'b', 'LineWidth', 1.5); 
xlabel('Time (s)'); ylabel('Normalized Force (%MVC)'); title('Normalized Force Signal'); 
hold on;
plot(time(new_locs_peaks), new_peaks, 'r*', 'LineWidth', 1.5);
hold on;
plot(time(new_locs_trough), new_trough, 'g*', 'LineWidth', 1.5);
hold off;

% Merge the peaks and troughs locations (new_locs)
new_locs = [new_locs_peaks; new_locs_trough];
new_locs = sort(new_locs);
new_locs = new_locs(1:count_peaks+count_trough);
new_locs = sort(new_locs);

% Assign 1 to peaks and 0 to troughs
new_vals = zeros(length(new_locs), 1);
for i = 1:length(new_locs)
    if ismember(new_locs(i), new_locs_peaks)
        new_vals(i) = 1;
    end
end

% Generate the segments
segments_loc = zeros(length(new_locs), 1);
count = 1;
extend = 0.25;
for i = 1:length(new_locs) - 1
    if new_vals(i) ~= new_vals(i+1)
        if new_vals(i) == 0
            segments_loc(count) = new_locs(i+1)-(extend*(new_locs(i+1)-new_locs(i)));
        else
            segments_loc(count) = new_locs(i)+(extend*(new_locs(i+1)-new_locs(i)));
        end
        count = count + 1;
    end
end
segments_loc = segments_loc(1:count-1);

% Calculate the mean force for each segment
segments_mean_force = zeros(length(segments_loc)/2, 1);
for i = 1:2:length(segments_loc)
    segments_mean_force((i+1)/2) = mean(normalized_force(segments_loc(i):segments_loc(i+1)));
end

% Adjust the starting point to be based on mean force
for i = 1:2:length(segments_loc)
    segments_loc(i) = segments_loc(i) - extend*(segments_mean_force((i+1)/2));
    segments_loc(i+1) = segments_loc(i+1) + extend*(segments_mean_force((i+1)/2));
end
    
% Plot the segments on the force signal
figure(4);
subplot(2, 1, 1);
plot(time, normalized_force, 'b'); % EMG signal in blue
xlabel('Time (s)'); ylabel('Normalized Force Signal (%MVC)'); title('Normalized Force Signal');
hold on;
% Highlight the segments (in red) on the force signal
for i = 1:2:length(segments_loc)
    plot(time(segments_loc(i):segments_loc(i+1)), normalized_force(segments_loc(i):segments_loc(i+1)), 'r', 'LineWidth', 1.5); 
end
% Plot mean force for each segment on the force signal
for i = 1:2:length(segments_loc)
    plot(time(segments_loc(i):segments_loc(i+1)), segments_mean_force((i+1)/2)*ones(length(emg(segments_loc(i):segments_loc(i+1), 1)), 1), 'k', 'LineWidth', 1.5);
end
hold off;

% Plot the segments on the filtered EMG signal
subplot(2, 1, 2);
plot(time, emg_filtered_hplp, 'b'); % EMG signal in blue
xlabel('Time (s)'); ylabel('EMG (mV)'); title('HPLP Filtered EMG Signal');
hold on;
for i = 1:2:length(segments_loc)
    plot(emg(segments_loc(i):segments_loc(i+1), 1), emg_filtered_hplp(segments_loc(i):segments_loc(i+1)), 'r'); % EMG signal in red
end
hold off;

%% ===================================================== PARAMETER CALCULATION =====================================================
% CALCULATION OF SEVERAL PARAMETERS BASED ON THE EMG SIGNAL SEGMENTS BEFORE AND AFTER FILTERING
% Create several empty arrays to store the value of the calculated parameters
segments_unfiltered_rms = zeros(length(segments_loc)/2, 1);
segments_rms = zeros(length(segments_loc)/2, 1);
segments_unfiltered_mean_freq = zeros(length(segments_loc)/2, 1);
segments_mean_freq = zeros(length(segments_loc)/2, 1);
segments_unfiltered_median_freq = zeros(length(segments_loc)/2, 1);
segments_median_freq = zeros(length(segments_loc)/2, 1);
segments_unfiltered_fractal_dim = zeros(length(segments_loc)/2, 1);
segments_fractal_dim = zeros(length(segments_loc)/2, 1);
segments_avgforce = zeros(length(segments_loc)/2, 1);

% Calculate the parameters for each segment before and after filtering
for i = 1:2:length(segments_loc)
    % 1. RMS value (A measure of the overall magnitude of the signal within that segment)
    segments_unfiltered_rms((i+1)/2) = rms(emg_sig(segments_loc(i):segments_loc(i+1)));
    segments_rms((i+1)/2) = rms(emg_filtered_hplp(segments_loc(i):segments_loc(i+1)));

    % 2. Mean frequency (A measure of the average frequency content of the signal within that segment)
    segments_unfiltered_mean_freq((i+1)/2) = meanfreq(emg_sig(segments_loc(i):segments_loc(i+1)), fs);
    segments_mean_freq((i+1)/2) = meanfreq(emg_filtered_hplp(segments_loc(i):segments_loc(i+1)), fs);

    % 3. Median frequency (A measure of the frequency content of the signal within that segment)
    segments_unfiltered_median_freq((i+1)/2) = medfreq(emg_sig(segments_loc(i):segments_loc(i+1)), fs);
    segments_median_freq((i+1)/2) = medfreq(emg_filtered_hplp(segments_loc(i):segments_loc(i+1)), fs);

    % 4. Higuchi Fractal Dimension (Describes the degree of self-similarity along a 1D time series, according to a sequencing integer, k)
    segments_unfiltered_fractal_dim((i+1)/2) = HigFracDim(emg_sig(segments_loc(i):segments_loc(i+1)), 10);
    segments_fractal_dim((i+1)/2) = HigFracDim(emg_filtered_hplp(segments_loc(i):segments_loc(i+1)), 10);

    % 5. Average force in %MVC of each segment (Average amplitude of each segment)
    segments_avgforce((i+1)/2) = mean(normalized_force(segments_loc(i):segments_loc(i+1)));
end

%% ============================================ LINEAR FITTING AND FEATURES CALCULATION =====================================================
% Use the polyfit function to obtain a straight-line (linear) fit to represent the variation
% of each EMG parameters versus force for both unfiltered and filtered EMG Signal
% Calculate its correlation coefficient as well
% 1. RMS value
p_unfiltered_rms = polyfit(segments_avgforce, segments_unfiltered_rms, 1);
p_rms = polyfit(segments_avgforce, segments_rms, 1);
r_unfiltered_rms = corrcoef(segments_avgforce, segments_unfiltered_rms);
r_rms = corrcoef(segments_avgforce, segments_rms);
% 2. Mean frequency
p_unfiltered_mean_freq = polyfit(segments_avgforce, segments_unfiltered_mean_freq, 1);
p_mean_freq = polyfit(segments_avgforce, segments_mean_freq, 1);
r_unfiltered_mean_freq = corrcoef(segments_avgforce, segments_unfiltered_mean_freq);
r_mean_freq = corrcoef(segments_avgforce, segments_mean_freq);
% 3. Median frequency
p_unfiltered_median_freq = polyfit(segments_avgforce, segments_unfiltered_median_freq, 1);
p_median_freq = polyfit(segments_avgforce, segments_median_freq, 1);
r_unfiltered_median_freq = corrcoef(segments_avgforce, segments_unfiltered_median_freq);
r_median_freq = corrcoef(segments_avgforce, segments_median_freq);
% 4. Fractal dimension
p_unfiltered_fractal_dim = polyfit(segments_avgforce, segments_unfiltered_fractal_dim, 1);
p_fractal_dim = polyfit(segments_avgforce, segments_fractal_dim, 1);
r_unfiltered_fractal_dim = corrcoef(segments_avgforce, segments_unfiltered_fractal_dim);
r_fractal_dim = corrcoef(segments_avgforce, segments_fractal_dim);
% Slope and intercept 
slope_unfiltered = [p_unfiltered_rms(1), p_unfiltered_mean_freq(1), p_unfiltered_median_freq(1), p_unfiltered_fractal_dim(1)]';
slope_filtered = [p_rms(1), p_mean_freq(1), p_median_freq(1), p_fractal_dim(1)]';

% Plot the values of the parameters computed above versus force (in %MVC)
% Superimpose the linear models (straight-line fits) obtained on the plots of the parameters in the preceding step. 
% Unfiltered EMG Signal
figure(5);
sgtitle('UNFILTERED EMG SIGNAL')
% 1. RMS value;
subplot(2, 2, 1);
plot(segments_avgforce, segments_unfiltered_rms, 'o');
hold on;
plot(segments_avgforce, p_unfiltered_rms(1)*segments_avgforce+p_unfiltered_rms(2), 'r-'); 
hold off;
xlabel('Average Force (%MVC)'); ylabel('RMS Value'); title('Force (%MVC) vs RMS Value');
% 2. Mean frequency
subplot(2, 2, 2);
plot(segments_avgforce, segments_unfiltered_mean_freq, 'o');
hold on;
plot(segments_avgforce, p_unfiltered_mean_freq(1)*segments_avgforce+p_unfiltered_mean_freq(2), 'r-'); hold off;
hold off;
xlabel('Average Force (%MVC)'); ylabel('Mean Frequency'); title('Force (%MVC) vs Mean Frequency');
% 3. Median frequency
subplot(2, 2, 3);
plot(segments_avgforce, segments_unfiltered_median_freq, 'o');
hold on;
plot(segments_avgforce, p_unfiltered_median_freq(1)*segments_avgforce+p_unfiltered_median_freq(2), 'r-'); hold off;
hold off;
xlabel('Average Force (%MVC)'); ylabel('Median Frequency'); title('Force (%MVC) vs Median Frequency');
% 4. Fractal dimension
subplot(2, 2, 4);
plot(segments_avgforce, segments_unfiltered_fractal_dim, 'o');
hold on;
plot(segments_avgforce, p_unfiltered_fractal_dim(1)*segments_avgforce+p_unfiltered_fractal_dim(2), 'r-'); hold off;
hold off;
xlabel('Average Force (%MVC)'); ylabel('Fractal Dimension'); title('Force (%MVC) vs Fractal Dimension');

% Filtered EMG Signal
figure(6);
sgtitle('FILTERED EMG SIGNAL')
% 1. RMS value;
subplot(2, 2, 1);
plot(segments_avgforce, segments_rms, 'o');
hold on;
plot(segments_avgforce, p_rms(1)*segments_avgforce+p_rms(2), 'r-'); 
hold off;
xlabel('Average Force (%MVC)'); ylabel('RMS Value'); title('Force (%MVC) vs RMS Value');
% 2. Mean frequency
subplot(2, 2, 2);
plot(segments_avgforce, segments_mean_freq, 'o');
hold on;
plot(segments_avgforce, p_mean_freq(1)*segments_avgforce+p_mean_freq(2), 'r-'); hold off;
hold off;
xlabel('Average Force (%MVC)'); ylabel('Mean Frequency'); title('Force (%MVC) vs Mean Frequency');
% 3. Median frequency
subplot(2, 2, 3);
plot(segments_avgforce, segments_median_freq, 'o');
hold on;
plot(segments_avgforce, p_median_freq(1)*segments_avgforce+p_median_freq(2), 'r-'); hold off;
hold off;
xlabel('Average Force (%MVC)'); ylabel('Median Frequency'); title('Force (%MVC) vs Median Frequency');
% 4. Fractal dimension
subplot(2, 2, 4);
plot(segments_avgforce, segments_fractal_dim, 'o');
hold on;
plot(segments_avgforce, p_fractal_dim(1)*segments_avgforce+p_fractal_dim(2), 'r-'); hold off;
hold off;
xlabel('Average Force (%MVC)'); ylabel('Fractal Dimension'); title('Force (%MVC) vs Fractal Dimension');

%% Tabulation of the slope, intercept, and correlation coefficient of the parameters 
fprintf('UNFILTERED EMG SIGNAL \n')
T1 = table({'RMS'; 'Mean frequency'; 'Median frequency'; 'Fractal dimension'}, ...
    [p_unfiltered_rms(1); p_unfiltered_mean_freq(1); p_unfiltered_median_freq(1); p_unfiltered_fractal_dim(1)], ...
    [p_unfiltered_rms(2); p_unfiltered_mean_freq(2); p_unfiltered_median_freq(2); p_unfiltered_fractal_dim(2)], ...
    [r_unfiltered_rms(1, 2); r_unfiltered_mean_freq(1, 2); r_unfiltered_median_freq(1, 2); r_unfiltered_fractal_dim(1, 2)], ...
    'VariableNames', {'EMG_Parameters', 'Slope', 'Interception', 'Correlation_Coefficient'});
disp(T1);

fprintf('FILTERED EMG SIGNAL \n')
T2 = table({'RMS'; 'Mean frequency'; 'Median frequency'; 'Fractal dimension'}, ...
    [p_rms(1); p_mean_freq(1); p_median_freq(1); p_fractal_dim(1)], ...
    [p_rms(2); p_mean_freq(2); p_median_freq(2); p_fractal_dim(2)], ...
    [r_rms(1, 2); r_mean_freq(1, 2); r_median_freq(1, 2); r_fractal_dim(1, 2)], ...
    'VariableNames', {'EMG_Parameters', 'Slope', 'Interception', 'Correlation_Coefficient'});
disp(T2);