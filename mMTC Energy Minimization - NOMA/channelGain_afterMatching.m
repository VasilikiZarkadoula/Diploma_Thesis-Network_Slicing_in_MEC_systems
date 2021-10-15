function g_matched = channelGain_afterMatching(SC_matched,channels,channel_gain)

A = zeros(1,length(SC_matched));
for i = 1:length(SC_matched)
    A(i) = (i-1)*channels + SC_matched(i);
end
g_matched = channel_gain(A);