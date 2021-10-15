function [gm,gu] = channelGain(radius_m,radius_u,a,Bn,No)

rng shuffle;

dm = radius_m*rand;
du = radius_m + (radius_u - radius_m)*rand;
h = exprnd(2);
gm = ((dm.^(-a)).*(abs(h).^2))./(Bn*No);
gu = ((du.^(-a)).*(abs(h).^2))./(Bn*No);