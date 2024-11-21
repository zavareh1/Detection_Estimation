% Simulation Parameters
K = 50;  % Number of secondary data vectors
N = 10;  % Dimension of each data vector
P_fa = 1e-6;  % Probability of false alarm
P_d = 0.9;  % Target probability of detection
SNR_dB_range = -10:2:20;  % SNR range in dB
C = 6;  % Loss factor from the paper (depends on specific cases)
num_trials = 1000;  % Number of trials for averaging

% Pre-allocate arrays for storing results
P_d_lrt = zeros(size(SNR_dB_range));  % Probability of detection for LRT
P_d_mf = zeros(size(SNR_dB_range));   % Probability of detection for matched filter

% Threshold for the likelihood ratio test (based on P_fa)
threshold_lrt = (1 / P_fa)^(1 / (K + 1 - N));

% Loop over SNR range and compute P_d for both tests
for idx = 1:length(SNR_dB_range)
    SNR_dB = SNR_dB_range(idx);
    snr_linear = 10^(SNR_dB / 10);  % Convert SNR from dB to linear scale
    
    % Simulate LRT (Likelihood Ratio Test)
    P_d_lrt(idx) = simulate_lrt(snr_linear, num_trials, K, N, threshold_lrt, C);
    
    % Simulate matched filter
    P_d_mf(idx) = simulate_matched_filter(snr_linear, num_trials, K, N, P_fa);
end

% Plotting the results
figure;
plot(SNR_dB_range, P_d_lrt, '-o', 'LineWidth', 2, 'DisplayName', 'LRT');
hold on;
plot(SNR_dB_range, P_d_mf, '--s', 'LineWidth', 2, 'DisplayName', 'Matched Filter');
xlabel('SNR (dB)');
ylabel('Probability of Detection (P_d)');
title('Probability of Detection vs SNR for LRT and Matched Filter');
grid on;
legend('Location', 'best');
hold off;

% Functions are declared below

% Function to simulate likelihood ratio test (LRT)
function P_d = simulate_lrt(SNR_linear, num_trials, K, N, threshold, C)
    detection_count = 0;
    for trial = 1:num_trials
        % Generate secondary data (noise only)
        noise_cov = eye(N);
        secondary_data = generate_complex_gaussian(zeros(1, N), noise_cov, K);
        
        % Compute S as the sum of the outer products of the secondary data vectors
        S = zeros(N, N);
        for k = 1:K
            S = S + (secondary_data(k, :)' * secondary_data(k, :));
        end
        
        % Generate primary data (signal + noise)
        signal_amplitude = sqrt(SNR_linear);
        signal = signal_amplitude * ones(N, 1);
        noise = generate_complex_gaussian(zeros(1, N), noise_cov, 1);
        primary_data = signal + noise.';

        % Compute the detection statistic (likelihood ratio test)
        S_inv = inv(S);
        z_s_inv_z = real(primary_data' * S_inv * primary_data);
        s_s_inv_s = real(signal' * S_inv * signal);
        eta = z_s_inv_z / (1 + z_s_inv_z - abs(signal' * S_inv * primary_data)^2 / s_s_inv_s);

        % Apply decision rule
        if eta > threshold
            detection_count = detection_count + 1;
        end
    end
    P_d = detection_count / num_trials;
end

% Function to simulate matched filter
function P_d = simulate_matched_filter(SNR_linear, num_trials, K, N, P_fa)
    detection_count = 0;
    threshold = compute_threshold_matched_filter(P_fa, N);  % Compute threshold based on P_fa
    
    for trial = 1:num_trials
        % Generate secondary data (noise only)
        noise_cov = eye(N);
        secondary_data = generate_complex_gaussian(zeros(1, N), noise_cov, K);
        
        % Generate primary data (signal + noise)
        signal_amplitude = sqrt(SNR_linear);
        signal = signal_amplitude * ones(N, 1);
        noise = generate_complex_gaussian(zeros(1, N), noise_cov, 1);
        primary_data = signal + noise.';

        % Compute matched filter statistic
        mf_output = abs(primary_data' * signal)^2;
        
        % Apply decision rule
        if mf_output > threshold
            detection_count = detection_count + 1;
        end
    end
    P_d = detection_count / num_trials;
end

% Function to compute the threshold for the matched filter
function threshold = compute_threshold_matched_filter(P_fa, N)
    threshold = chi2inv(1 - P_fa, 2 * N);  % Threshold from chi-squared distribution
end

% Helper function to generate complex Gaussian vectors
function data = generate_complex_gaussian(mean_vec, cov_matrix, num_samples)
    real_part = mvnrnd(real(mean_vec), real(cov_matrix), num_samples);
    imag_part = mvnrnd(imag(mean_vec), imag(cov_matrix), num_samples);
    data = real_part + 1i * imag_part;
end
