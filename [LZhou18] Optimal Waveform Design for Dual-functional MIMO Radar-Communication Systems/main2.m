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

a = zeros(p.N, 1);
for idx = 1 : p.N
    a(idx, 1) = exp(1i * pi * (idx - ceil(p.N/2)) * sin(p.theta(127)));     % p.theta = pi/5;
end

end