function [] = main1 ()

% Figure 4

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

% Simulation Settings
p.iterations = 1000;

OmniStrictCapacityArray = zeros(p.iterations, length(p.SNRdB));
OmniTradeoffCapacityArray = zeros(p.iterations, length(p.SNRdB));

OmniStrictBPArray = zeros(p.iterations, length(p.theta));
OmniTradeoffBPArray = zeros(p.iterations, length(p.theta));

end