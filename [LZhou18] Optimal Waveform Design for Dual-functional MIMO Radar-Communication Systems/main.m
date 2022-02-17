function [] = main ()

% Parameter Settings
p.K = 30;       % # of Users
p.N = 64;       % # of Antennas per Each Users (ULA)
p.T = 70;       % # of Communication Frame

p.Pt = 1;                   % Total Power Constraint
p.N0dB = -5 : 1 : 5;        % Noise Settings
p.N0 = 10.^(p.N0dB / 10);   % SNR Settings


% Channel Realization
H = (1/sqrt(2)) * (randn([p.K, p.N]) + 1i * randn([p.K, p.N]));

% Omni-Directional Beampattern Design
Rd = () * eye(p.N);


end