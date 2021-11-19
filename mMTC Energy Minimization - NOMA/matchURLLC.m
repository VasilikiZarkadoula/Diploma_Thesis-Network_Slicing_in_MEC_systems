function g_URLLC_matched = matchURLLC(users,d_URLLC,a,Bs,No)
% This functions matches URLLC users to subchannels based on their channel gain

rng shuffle;

subchannels = users/2;
h_URLLC = zeros(users/2,subchannels); 
g_URLLC = zeros(users/2,subchannels); 
for k=1:users/2
    for s=1:subchannels
        % exponential channel gain corresponding to Rayleigh fading
        h_URLLC(k,s) = exprnd(2);
    end
    % channel gain
    g_URLLC = ((d_URLLC(k).^(-a)).*(abs(h_URLLC).^2))./(Bs*No);
end

% matching
G = g_URLLC;
SC_URLLC_unmatched = zeros(1,subchannels);
while any(SC_URLLC_unmatched == 0)
    maxG = max(G(:));
    [channel,user] = find(G == maxG);
    SC_URLLC_unmatched(channel) = user;
    G(:, user) = 0;
    G(channel, :) = 0;
end
SC_URLLC_matched = SC_URLLC_unmatched;

g_URLLC_matched = channelGain_afterMatching(SC_URLLC_matched,subchannels,g_URLLC);

