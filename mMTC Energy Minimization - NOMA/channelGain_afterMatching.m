function g_matched = channelGain_afterMatching(SC_matched,channels,channel_gain)
% This function returns the channel gains of users assigned in subchannels,
% after the optimal assignment has been made
% In detail, the first element of vector g_matched conntains the channel
% gain of the user assigned to subchannel 1, the second element of g_matched 
% conntains the channel gain of the user assigned to subchannel 2 etc.

A = zeros(1,length(SC_matched));
for i = 1:length(SC_matched)
    A(i) = (i-1)*channels + SC_matched(i);
end
g_matched = channel_gain(A);
