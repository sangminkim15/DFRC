function [] = main ()

% Figure 2, 3

% Parameter Settings
% Communication Settings
p.K = 4;                            % # of Users
p.N = 16;                           % # of Antennas per Each Users (ULA)
p.L = 20;                           % # of Communication Frame
p.Pt = 1;                           % Total Power Constraint
p.N0dB = 2 : -2 : -12;              % Noise Settings
p.N0 = 10.^(p.N0dB ./ 10);  
p.SNR = p.Pt ./ p.N0;
p.SNRdB = 10 * log(p.SNR) / log(10);    % SNR Settings

% Radar Settings
p.theta = -pi/2 : pi/180 : pi/2;        % Radar ULA Angle Settings
p.theta_target = [-pi*10/180, -pi*5/180, 0, pi*5/180, pi*10/180];
p.target_DoA = [-pi/3,0,pi/3];

p.beam_width= 9;
p.l=ceil((p.target_DoA + pi/2 * ones(1, length(p.target_DoA)))/(pi/180) + ones(1, length(p.target_DoA)));
p.Pd_theta = zeros(length(p.theta), 1);

for idx = 1:length(p.target_DoA)
    p.Pd_theta(p.l(idx)-(p.beam_width-1)/2 : p.l(idx)+(p.beam_width-1)/2, 1) = ones(p.beam_width, 1);
end

p.c = 3e8;
p.fc = 3.2e9;
p.lambda = p.c / p.fc;
p.spacing = p.lambda / 2;

% Tradeoff Settings
p.rho = 0.2;                            % Weighting Factor

% Omni-Directional Beampattern
OmniRd = (p.Pt / p.N) * eye(p.N, p.N);

% Simulation Settings
p.montecarlo = 1000;

OmniStrictCapacityArray = zeros(p.montecarlo, length(p.SNRdB));
OmniTradeoffCapacityArray = zeros(p.montecarlo, length(p.SNRdB));
OmniStrictBPArray = zeros(p.montecarlo, length(p.theta));
OmniTradeoffBPArray = zeros(p.montecarlo, length(p.theta));

for idx = 1 : p.montecarlo
    % Channel Realization
    H = (1/sqrt(2)) * (randn([p.K, p.N]) + 1i * randn([p.K, p.N]));

    % Desired Signal Matrix - 4QAM Modulation
    S = (1/sqrt(2)) * ((2 * randi([0 1], p.K, p.L) - ones(p.K, p.L)) + 1i * ((2 * randi([0 1], p.K, p.L) - ones(p.K, p.L))));

    % Optimal Waveform Design (Omni-Directional)
    F = chol(OmniRd);                   % Cholesky Factorization
    [U, ~, V] = svd(F * H' * S);        % SVD (singular value decomposition)
    
    OmniXStrict = sqrt(p.L) * F' * U * eye(p.N, p.L) * V';

    % Trade-off Between Radar and Communication Performances
    Q = p.rho * (H' * H) + (1 - p.rho) * eye(p.N, p.N);
    G = p.rho * H' * S + (1 - p.rho) * OmniXStrict;

    [P, LAMBDA] = eig(Q);               % Eigenvalue, Eigenvector of Matrix Q

    lambda_low = - min(diag(LAMBDA));
    lambda_high = - min(diag(LAMBDA)) + sqrt(p.N / p.Pt) * max(max(abs(P' * G)));

    OmniXTradeoff = BisectionSearch(Q, G, lambda_low, lambda_high, p);
    
    % Communication Capacity
    for jdx = 1 : length(p.N0dB)
        OmniEStrict = H * OmniXStrict - S;
        OmniETradeoff = H * OmniXTradeoff - S;
        OmnigammaStrict = 1 ./ (mean(abs(OmniEStrict).^2, 2) + p.N0(jdx));
        OmnigammaTradeoff =  1 ./ (mean(abs(OmniETradeoff).^2, 2) + p.N0(jdx)); 
        
        for kdx = 1 : p.K
            OmniStrictCapacityArray(idx, jdx) = OmniStrictCapacityArray(idx, jdx) + log(1 + OmnigammaStrict(kdx)) / log(2);
            OmniTradeoffCapacityArray(idx, jdx) = OmniTradeoffCapacityArray(idx, jdx) + log(1 + OmnigammaTradeoff(kdx)) / log(2);
        end
    end
    
    % Radar Beampattern
    for jdx = 1 : length(p.theta)
        a = zeros(p.N, 1);
        
        for kdx = 1 : p.N
            a(kdx, 1) = exp(1i * pi * (kdx - ceil(p.N / 2)) * sin(p.theta(jdx)));
        end
        
        OmniStrictBPArray(idx, jdx) = a' * (OmniXStrict * OmniXStrict') * a / real(trace(OmniXStrict * OmniXStrict'));
        OmniTradeoffBPArray(idx, jdx) = a' * (OmniXTradeoff * OmniXTradeoff') * a / real(trace(OmniXTradeoff * OmniXTradeoff'));
    end
end

% Communication Plot
AWGNCapacity = p.K * log(1 + p.SNR) / log(2);
OmniStrictCapacity = mean(OmniStrictCapacityArray);
OmniTradeoffCapacity = mean(OmniTradeoffCapacityArray);

figure
plot(p.SNRdB, AWGNCapacity, 'r--v', p.SNRdB, OmniStrictCapacity, 'k-x', p.SNRdB, OmniTradeoffCapacity, 'k--o', 'LineWidth', 1.5);
xlabel('Transmit SNR (dB)');
ylabel('Average Achievable Sum Rate (bps/Hz)');
legend('AWGN Capacity', 'Omni-Strict', 'Omni-Tradeoff (\rho = 0.2)', 'Location', 'northwest');
grid on

% Radar Plot
OmniStrictBP = real(OmniStrictBPArray(1,:));
OmniStrictBP = 10 .* log10(OmniStrictBP);
OmniTradeoffBP = real(OmniTradeoffBPArray(1,:));
OmniTradeoffBP = 10 .* log10(OmniTradeoffBP);

p.theta_deg = p.theta * (180/pi);

figure
plot(p.theta_deg, OmniStrictBP, 'r--', p.theta_deg, OmniTradeoffBP, 'k-', 'LineWidth', 1.5);
xlabel('\theta (deg)');
ylabel('Beampattern');
xlim([-90 90]);
legend('Omni-Strict', 'Omni-Tradeoff (\rho = 0.2)');
grid on

end