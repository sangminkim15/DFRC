function [] = main ()

% Parameter Settings
p.K = 30;                           % # of Users
p.N = 64;                           % # of Antennas per Each Users (ULA)
p.L = 70;                           % # of Communication Frame
p.Pt = 1;                           % Total Power Constraint
p.N0dB = 5 : -1 : -5;               % Noise Settings
p.N0 = 10.^(p.N0dB ./ 10);  
p.SNR = p.Pt ./ p.N0;
p.SNRdB = 10 * log(p.SNR ./ 10);    % SNR Settings

p.theta = -pi/4 : pi/180 : pi/4;    % Radar ULA Angle Settings

p.rho = 0.1;                        % Weighting Factor

% Simulation Settings
p.iterations = 1000;

OmniStrictCapacityArray = zeros(p.iterations, length(p.SNRdB));
OmniTradeoffCapacityArray = zeros(p.iterations, length(p.SNRdB));
OmniStrictBPArray = zeros(p.iterations, length(p.theta));
OmniTradeoffBPArray = zeros(p.iterations, length(p.theta));

for idx = 1 : p.iterations
    % Channel Realization
    H = (1/sqrt(2)) * (randn([p.K, p.N]) + 1i * randn([p.K, p.N]));

    % Desired Signal Matrix - 4QAM Modulation
    S = (1/sqrt(2)) * ((2 * randi([0 1], p.K, p.L) - ones(p.K, p.L)) + 1i * ((2 * randi([0 1], p.K, p.L) - ones(p.K, p.L))));

    % Omni-Directional Beampattern Design
    OmniRd = (p.Pt / p.N) * eye(p.N, p.N);

    % Optimal Waveform Design
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
    
    
    
    % Radar Beampattern
    for jdx = 1 : length(p.theta)
        a = zeros(p.N, 1);
        
        for kdx = 1 : p.N
            a(kdx, 1) = exp(1i * pi * (kdx - 1) * sin(p.theta(jdx)));
        end
        
        OmniStrictBPArray(idx, jdx) = a' * (1/p.L) * (OmniXStrict * OmniXStrict') * a;
        OmniTradeoffBPArray(idx, jdx) = a' * (1/p.L) * (OmniXTradeoff * OmniXTradeoff') * a;
        
    end
    
end

% Communication Plot
AWGNCapacity = log(1 + p.SNRdB);
OmniStrictCapacity = mean(OmniStrictCapacityArray);
OmniTradeoffCapacity = mean(OmniTradeoffCapacityArray);

% figure
% plot(p.SNRdB, AWGNCapacity, p.SNRdB, OmniStrictCapacity, p.SNRdB, OmniTradeoffCapacity, 'LineWidth', 1.5);
% xlabel('Transmit SNR (dB)');
% ylabel('Average Achievable Rate (bps/Hz/user)');
% legend('AWGN Capacity', 'Omni-Strict', 'Omni-Tradeoff (\rho = 0.1)');

% Radar Plot
OmniStrictBP = mean(abs(OmniStrictBPArray));
OmniTradeoffBP = mean(abs(OmniTradeoffBPArray));

figure
plot(p.theta, OmniStrictBP, p.theta, OmniTradeoffBP, 'LineWidth', 1.5);
xlabel('Angle (rad)');
ylabel('Beampattern');
legend('Omni-Strict', 'Omni-Tradeoff (\rho = 0.1)');

end