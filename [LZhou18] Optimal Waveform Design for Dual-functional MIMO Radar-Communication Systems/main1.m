function [] = main1 ()

% Figure 4

% Parameter Settings
p.K = 25 : 5 : 35;                      % # of Users
p.N = 64;                               % # of Antennas per Each Users (ULA)
p.L = 70;                               % # of Communication Frame
p.Pt = 1;                               % Total Power Constraint
p.N0dB = -5;                            % Noise Settings
p.N0 = 10.^(p.N0dB ./ 10);  
p.SNR = p.Pt ./ p.N0;
p.SNRdB = 10 * log(p.SNR) / log(10);    % SNR Settings

p.theta = pi/5;
p.alphadB = 19;
p.alpha = 10.^(p.alphadB / 10);

p.rho = 0 : 0.1 : 1;                    % Weighting Factor

% Simulation Settings
p.iterations = 1000;

OmniRateArray = zeros(p.iterations, length(p.rho), length(p.K));
OmniProbabilityArray = zeros(p.iterations, length(p.rho), length(p.K));

for idx = 1 : length(iterations)
    for jdx = 1 : length(p.rho)
        for kdx = 1 : length(p.K)
            % Channel Realization
            H = (1/sqrt(2)) * (randn([p.K(kdx), p.N]) + 1i * randn([p.K(kdx), p.N]));

            % Desired Signal Matrix - 4QAM Modulation
            S = (1/sqrt(2)) * ((2 * randi([0 1], p.K(kdx), p.L) - ones(p.K(kdx), p.L)) + 1i * ((2 * randi([0 1], p.K(kdx), p.L) - ones(p.K(kdx), p.L))));

            % Omni-Directional Beampattern Design
            OmniRd = (p.Pt / p.N) * eye(p.N, p.N);

            % Optimal Waveform Design
            F = chol(OmniRd);                   % Cholesky Factorization
            [U, ~, V] = svd(F * H' * S);        % SVD (singular value decomposition)
    
            OmniXStrict = sqrt(p.L) * F' * U * eye(p.N, p.L) * V';

            % Trade-off Between Radar and Communication Performances
            Q = p.rho(jdx) * (H' * H) + (1 - p.rho(jdx)) * eye(p.N, p.N);
            G = p.rho(jdx) * H' * S + (1 - p.rho(jdx)) * OmniXStrict;

            [P, LAMBDA] = eig(Q);               % Eigenvalue, Eigenvector of Matrix Q

            lambda_low = - min(diag(LAMBDA));
            lambda_high = - min(diag(LAMBDA)) + sqrt(p.N / p.Pt) * max(max(abs(P' * G)));

            OmniXTradeoff = BisectionSearch(Q, G, lambda_low, lambda_high, p);
            
            % Communication Rate
            OmniETradeoff = H * OmniXTradeoff - S;
            OmnigammaTradeoff =  1 / (mean(abs(OmniETradeoff).^2, 2) + p.N0);
            
            for ldx = 1 : p.K
                OmniRateArray(idx, jdx, kdx) = OmniRateArray(idx, jdx, kdx) + log(1 + OmnigammaTradeoff(kdx)) / log(2);
            end
            
            % Detection Probability
            
            
        end
    end
end
    
    
    
    
end


end