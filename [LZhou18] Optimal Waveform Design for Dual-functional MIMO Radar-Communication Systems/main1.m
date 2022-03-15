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
p.alphadB = -19;
p.alpha = 10.^(p.alphadB / 10);
p.falsealarm = 1e-7;

p.rhodB = [-30, -25, -20, -15, -10 : 2 : -2, -1, -0.5];
p.rho = 10.^(p.rhodB / 10);                     % Weighting Factor

% Simulation Settings
p.iterations = 5000;

OmniRateArray = zeros(p.iterations, length(p.rho), length(p.K));
OmniProbabilityArray = zeros(p.iterations, length(p.rho), length(p.K));

for idx = 1 : p.iterations
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
            OmnigammaTradeoff =  1 ./ (mean(abs(OmniETradeoff).^2, 2) + p.N0);
            
            temp = p.K(kdx);
            for ldx = 1 : temp
                OmniRateArray(idx, jdx, kdx) = OmniRateArray(idx, jdx, kdx) + log(1 + OmnigammaTradeoff(ldx)) ./ log(2);
            end
            
            OmniRateArray(idx, jdx, kdx) = OmniRateArray(idx, jdx, kdx) / p.K(kdx);
            
            % Detection Probability
            a = zeros(p.N, 1);
            for ldx = 1 : p.N
                a(ldx, 1) = exp(1i * pi * (ldx - 1) * sin(p.theta));
            end
            
            % Noncentrality needs to be fixed (line 75)
            p.noncentrality = p.alpha * p.L^2 * abs(a' * (1/p.L) * (OmniXTradeoff * OmniXTradeoff') * a);
            
            delta = chi2inv(1 - p.falsealarm, 2);
            OmniProbabilityArray(idx, jdx, kdx) = 1 - ncx2cdf(delta, 2, p.noncentrality);
        end
    end
end

OmniRate = mean(real(OmniRateArray));
OmniProbability = mean(real(OmniProbabilityArray));

OmniRate1 = OmniRate(:,:,1);
OmniRate2 = OmniRate(:,:,2);
OmniRate3 = OmniRate(:,:,3);

OmniProbability1 = OmniProbability(:,:,1);
OmniProbability2 = OmniProbability(:,:,2);
OmniProbability3 = OmniProbability(:,:,3);
    
figure
plot(OmniRate1, OmniProbability1, 'b-o', OmniRate2, OmniProbability2, 'k-x', OmniRate3, OmniProbability3, 'r-s', 'LineWidth', 1.5);
xlabel('Average Achievable Rate');
ylabel('Detection Probability');
legend('K=25', 'K=30', 'K=35', 'Location', 'southwest');
grid on
    
end