function [] = main1 ()

% Figure 4

% Parameter Settings
p.K = 4 : 2 : 8;                        % # of Users
p.N = 16;                               % # of Antennas per Each Users (ULA)
p.L = 20;                               % # of Communication Frame
p.Pt = 1;                               % Total Power Constraint
p.N0dB = -10;                           % Noise Settings
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

p.alphadB = -6;
p.alpha = 10.^(p.alphadB / 10);
p.falsealarm = 1e-7;

a = zeros(p.N, 1);
for idx = 1 : p.N
    a(idx, 1) = exp(1i * pi * (idx - ceil(p.N/2)) * sin(p.theta(127))); % p.theta = pi/5;
end

p.rho = zeros(1, 19);                   % Weighting Factor
for idx = 1 : length(p.rho)
    p.rho(idx) = idx / 20;
end

% Simulation Settings
p.montecarlo = 1000;

OmniRateArray = zeros(p.montecarlo, length(p.rho), length(p.K));
OmniProbabilityArray = zeros(p.montecarlo, length(p.rho), length(p.K));

for idx = 1 : p.montecarlo
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

            OmniXTradeoff = sqrt(p.N) * BisectionSearch(Q, G, lambda_low, lambda_high, p);
            
            % Communication Rate
            OmniETradeoff = H * (OmniXTradeoff / sqrt(p.N)) - S;
            OmnigammaTradeoff =  1 ./ (mean(abs(OmniETradeoff).^2, 2) + p.N0);
            
            temp = p.K(kdx);
            for ldx = 1 : temp
                OmniRateArray(idx, jdx, kdx) = OmniRateArray(idx, jdx, kdx) + log(1 + OmnigammaTradeoff(ldx)) ./ log(2);
            end
            OmniRateArray(idx, jdx, kdx) = OmniRateArray(idx, jdx, kdx) / p.K(kdx);
            
            Rs = OmniXTradeoff * OmniXTradeoff' / p.L;
            p.noncentrality = p.alpha * abs(a' * Rs.' * a)^2;
            
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
plot(OmniRate1, OmniProbability1, 'b', OmniRate2, OmniProbability2, 'k', OmniRate3, OmniProbability3, 'r', 'LineWidth', 1.5);
xlabel('Average Achievable Rate (bps/Hz/user)');
ylabel('Detection Probability');
legend('K=4', 'K=6', 'K=8', 'Location', 'southwest');
grid on
    
end