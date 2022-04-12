function [] = main2 ()

% Figure 5

% Communication Settings
p.K = 4 : 2 : 8;                        % # of Users
p.N = 16;                               % # of Antennas per Each Users (ULA)
p.L = 20;                               % # of Communication Frame
p.Pt = 1;                               % Total Power Constraint
p.N0dB = -10;                           % Noise Settings
p.N0 = 10.^(p.N0dB ./ 10);  
p.SNR = p.Pt ./ p.N0;
p.SNRdB = 10 * log(p.SNR) / log(10);    % SNR Settings

% Radar Settings
p.theta = -pi/2 : pi/180 : pi/2;
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

p.rhodB = [-30 -25 -20 -15 -10 -8 -6 -4 -2 -1];                                 % Weighting Factor
p.rho = 10.^(p.rhodB ./ 10);

% Directional Beampattern
DirectRd = directbeampattern(p);

% Simulation Settings
p.montecarlo = 1000;

DirectBPArray = zeros(p.montecarlo, length(p.rho), length(p.K));
DirectStrictBPArray = zeros(p.montecarlo, length(p.rho), length(p.K));
DirectTradeoffBPArray = zeros(p.montecarlo, length(p.rho), length(p.K));

for idx = 1 : p.montecarlo
    for jdx = 1 : length(p.rho)
        for kdx = 1 : length(p.K)
            % Channel Realization
            H = (1/sqrt(2)) * (randn([p.K(kdx), p.N]) + 1i * randn([p.K(kdx), p.N]));

            % Desired Signal Matrix - 4QAM Modulation
            S = (1/sqrt(2)) * ((2 * randi([0 1], p.K(kdx), p.L) - ones(p.K(kdx), p.L)) + 1i * ((2 * randi([0 1], p.K(kdx), p.L) - ones(p.K(kdx), p.L))));

            % Optimal Waveform Design
            F = chol(DirectRd);                 % Cholesky Factorization
            [U, ~, V] = svd(F * H' * S);        % SVD (singular value decomposition)
    
            DirectXStrict = sqrt(p.L) * F' * U * eye(p.N, p.L) * V';

            % Trade-off Between Radar and Communication Performances
            Q = p.rho(jdx) * (H' * H) + (1 - p.rho(jdx)) * eye(p.N, p.N);
            G = p.rho(jdx) * H' * S + (1 - p.rho(jdx)) * DirectXStrict;
            
            [P, LAMBDA] = eig(Q);               % Eigenvalue, Eigenvector of Matrix Q

            lambda_low = - min(diag(LAMBDA));
            lambda_high = - min(diag(LAMBDA)) + sqrt(p.N / p.Pt) * max(max(abs(P' * G)));

            DirectXTradeoff = sqrt(p.N) * BisectionSearch(Q, G, lambda_low, lambda_high, p);
            
            % Radar Beampattern
            for ldx = 1 : length(p.theta)
                a = zeros(p.N, 1);
        
                for lldx = 1 : p.N
                    a(lldx, 1) = exp(1i * pi * (kdx - ceil(p.N / 2)) * sin(p.theta(ldx)));
                end
                
                DirectStrictBPArray(idx, jdx) = a' * (DirectXStrict * DirectXStrict') * a / real(trace(DirectXStrict * DirectXStrict'));
                DirectTradeoffBPArray(idx, jdx) = a' * (DirectXTradeoff * DirectXTradeoff') * a / real(trace(DirectXTradeoff * DirectXTradeoff'));
            end
            
            
        end
    end
end




end