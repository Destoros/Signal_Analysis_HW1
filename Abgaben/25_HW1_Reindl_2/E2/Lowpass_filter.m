function Hd = Lowpass_filter
%LOWPASS_FILTER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.3 and Signal Processing Toolbox 7.5.
% Generated on: 29-Nov-2019 16:25:39

% Equiripple Lowpass filter designed using the FIRPM function.

% All frequency values are normalized to 1.

N     = 100;   % Order
Fpass = 0.48;  % Passband Frequency
Fstop = 0.52;  % Stopband Frequency
Wpass = 1;     % Passband Weight
Wstop = 1;     % Stopband Weight
dens  = 20;    % Density Factor

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, [0 Fpass Fstop 1], [1 1 0 0], [Wpass Wstop], {dens});
Hd = dfilt.dffir(b);

% [EOF]
