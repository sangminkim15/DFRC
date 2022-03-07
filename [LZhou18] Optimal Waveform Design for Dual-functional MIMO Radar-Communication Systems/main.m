function [] = main ()

% Figure 2, 3

% Parameter Settings
p.K = 30;                           % # of Users
p.N = 64;                           % # of Antennas per Each Users (ULA)
p.L = 70;                           % # of Communication Frame
p.Pt = 1;                           % Total Power Constraint
p.N0dB = 5 : -1 : -5;               % Noise Settings
p.N0 = 10.^(p.N0dB ./ 10);  
p.SNR = p.Pt ./ p.N0;
p.SNRdB = 10 * log(p.SNR) / log(10);    % SNR Settings

p.theta = -pi/2 : pi/180 : pi/2;        % Radar ULA Angle Settings

p.rho = 0.1;                            % Weighting Factor

% Omni-Directional Beampattern
OmniRd = (p.Pt / p.N) * eye(p.N, p.N);

% Directional Beampattern (needs to bge fixed)
DirectRd = zeros(p.N, p.N);

for idx = 1 : p.N
    for jdx = 1 : p.N
        DirectRd(idx, jdx) = 0.9999.^abs(idx - jdx);
    end
end

a1 = zeros(p.N, 1);
a2 = zeros(p.N, 1);

for idx = 1 : p.N
    a1(idx, 1) = exp(1i * pi * (idx - 1) * sin(+pi/3));
end

for idx = 1 : p.N
    a2(idx, 1) = exp(1i * pi * (idx - 1) * sin(-pi/3));
end

DirectRd1 = DirectRd * (a1 * a1');
DirectRd2 = DirectRd * (a2 * a2');
DirectRd = DirectRd + DirectRd1 + DirectRd2;
DirectRd = DirectRd * (sqrt(p.Pt / p.N) / norm(DirectRd, 'fro'));

% Simulation Settings
p.iterations = 1000;

OmniStrictCapacityArray = zeros(p.iterations, length(p.SNRdB));
OmniTradeoffCapacityArray = zeros(p.iterations, length(p.SNRdB));
DirectStrictCapacityArray = zeros(p.iterations, length(p.SNRdB));
DirectTradeoffCapacityArray = zeros(p.iterations, length(p.SNRdB));

OmniStrictBPArray = zeros(p.iterations, length(p.theta));
OmniTradeoffBPArray = zeros(p.iterations, length(p.theta));
DirectStrictBPArray = zeros(p.iterations, length(p.theta));
DirectTradeoffBPArray = zeros(p.iterations, length(p.theta));

for idx = 1 : p.iterations
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
    
    % Optimal Waveform Design (Directional)
    [L, D] = ldl(DirectRd);
    F = sqrt(D) * L';                                   % Cholesky Factorization
    [U, ~, V] = svd(F * H' * S);                        % SVD (singular value decomposition)
    
    DirectXStrict = sqrt(p.L) * F' * U * eye(p.N, p.L) * V';

    % Trade-off Between Radar and Communication Performances
    Q = p.rho * (H' * H) + (1 - p.rho) * eye(p.N, p.N);
    G = p.rho * H' * S + (1 - p.rho) * DirectXStrict;

    [P, LAMBDA] = eig(Q);               % Eigenvalue, Eigenvector of Matrix Q

    lambda_low = - min(diag(LAMBDA));
    lambda_high = - min(diag(LAMBDA)) + sqrt(p.N / p.Pt) * max(max(abs(P' * G)));

    DirectXTradeoff = BisectionSearch(Q, G, lambda_low, lambda_high, p);
    
    % Communication Capacity
    for jdx = 1 : length(p.N0dB)
        OmniEStrict = H * OmniXStrict - S;
        OmniETradeoff = H * OmniXTradeoff - S;
        OmnigammaStrict = 1 ./ (mean(abs(OmniEStrict).^2, 2) + p.N0(jdx));
        OmnigammaTradeoff =  1 ./ (mean(abs(OmniETradeoff).^2, 2) + p.N0(jdx));
        
        DirectEStrict = H * DirectXStrict - S;
        DirectETradeoff = H * DirectXTradeoff - S;
        DirectgammaStrict = 1 ./ (mean(abs(DirectEStrict).^2, 2) + p.N0(jdx));
        DirectgammaTradeoff =  1 ./ (mean(abs(DirectETradeoff).^2, 2) + p.N0(jdx));        
        
        for kdx = 1 : p.K
            OmniStrictCapacityArray(idx, jdx) = OmniStrictCapacityArray(idx, jdx) + (1/p.K) * log(1 + OmnigammaStrict(kdx)) / log(2);
            OmniTradeoffCapacityArray(idx, jdx) = OmniTradeoffCapacityArray(idx, jdx) + (1/p.K) * log(1 + OmnigammaTradeoff(kdx)) / log(2);
            
            DirectStrictCapacityArray(idx, jdx) = OmniStrictCapacityArray(idx, jdx) + (1/p.K) * log(1 + DirectgammaStrict(kdx)) / log(2);
            DirectTradeoffCapacityArray(idx, jdx) = OmniTradeoffCapacityArray(idx, jdx) + (1/p.K) * log(1 + DirectgammaTradeoff(kdx)) / log(2);
        end
    end
    
    % Radar Beampattern
    for jdx = 1 : length(p.theta)
        a = zeros(p.N, 1);
        
        for kdx = 1 : p.N
            a(kdx, 1) = exp(1i * pi * (kdx - 1) * sin(p.theta(jdx)));
        end
        
        OmniStrictBPArray(idx, jdx) = 10 * log(a' * (1/p.L) * (OmniXStrict * OmniXStrict') * a) / log(10);
        OmniTradeoffBPArray(idx, jdx) = 10 * log(a' * (1/p.L) * (OmniXTradeoff * OmniXTradeoff') * a) / log(10);
        DirectStrictBPArray(idx, jdx) = 10 * log(a' * (1/p.L) * (DirectXStrict * DirectXStrict') * a) / log(10);
        DirectTradeoffBPArray(idx, jdx) = 10 * log(a' * (1/p.L) * (DirectXTradeoff * DirectXTradeoff') * a) / log(10);        
    end
end

% Communication Plot
AWGNCapacity = log(1 + p.SNR) / log(2);
OmniStrictCapacity = mean(OmniStrictCapacityArray);
OmniTradeoffCapacity = mean(OmniTradeoffCapacityArray);

DirectStrictCapacity = mean(DirectStrictCapacityArray);
DirectTradeoffCapacity = mean(DirectTradeoffCapacityArray);

figure
plot(p.SNRdB, AWGNCapacity, 'r--v', p.SNRdB, OmniStrictCapacity, 'b-x', p.SNRdB, OmniTradeoffCapacity, 'b--o', p.SNRdB, DirectStrictCapacity, 'k-x', p.SNRdB, DirectTradeoffCapacity, 'k--o', 'LineWidth', 1.5);
xlabel('Transmit SNR (dB)');
ylabel('Average Achievable Rate (bps/Hz/user)');
legend('AWGN Capacity', 'Omni-Strict', 'Omni-Tradeoff (\rho = 0.1)', 'Directional-Strict', 'Directional-Tradeoff (\rho = 0.1)', 'Location', 'northwest');
grid on

% Radar Plot
OmniStrictBP = mean(real(OmniStrictBPArray));
OmniTradeoffBP = mean(real(OmniTradeoffBPArray));
DirectStrictBP = mean(real(DirectStrictBPArray));
DirectTradeoffBP = mean(real(DirectTradeoffBPArray));

p.theta_deg = p.theta * (180/pi);

figure
plot(p.theta_deg, OmniStrictBP, 'r-', p.theta_deg, OmniTradeoffBP, 'r--', p.theta_deg, DirectStrictBP, 'b-', p.theta_deg, DirectTradeoffBP, 'b--', 'LineWidth', 1.5);
xlabel('\theta (deg)');
ylabel('Beampattern');
xlim([-90 90]);
legend('Omni-Strict', 'Omni-Tradeoff (\rho = 0.1)', 'Direct-Strict', 'Direct-Tradeoff (\rho = 0.1)');
grid on

end