function g_k = channelGain(radius,a,K)

d_k = radius*rand(K,1);
h_k = exprnd(2,K,1);
g_k = (d_k.^(-a)).*(abs(h_k).^2);
