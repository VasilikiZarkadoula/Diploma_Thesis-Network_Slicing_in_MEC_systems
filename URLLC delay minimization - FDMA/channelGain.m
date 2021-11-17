function g = channelGain(radius,a,K)

% distance between user - BS
d = radius*rand(K,1);
% exponential channel gain corresponding to Rayleigh fading
h = exprnd(2,K,1);
% channel gain
g = (d.^(-a)).*(abs(h).^2);