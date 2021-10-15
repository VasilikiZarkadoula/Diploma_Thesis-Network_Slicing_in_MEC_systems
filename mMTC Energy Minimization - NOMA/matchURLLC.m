function g_weak_matched = matchURLLC(users,channels,d_weak,a,Bs,No)

rng shuffle;
h_weak = zeros(users/2,channels); 
g_weak = zeros(users/2,channels); 
for k=1:users/2
    for s=1:channels
        h_weak(k,s) = exprnd(2);
    end
    g_weak = ((d_weak(k).^(-a)).*(abs(h_weak).^2))./(Bs*No);
end

%%% match weak users to sub-channels according to their channel gains %%%
G = g_weak;
SC_weak_unmatched = zeros(1,channels);
while any(SC_weak_unmatched == 0)
    maxG = max(G(:));
    [channel,user] = find(G == maxG);
    SC_weak_unmatched(channel) = user;
    G(:, user) = 0;
    G(channel, :) = 0;
end
SC_weak_matched = SC_weak_unmatched;

g_weak_matched = channelGain_afterMatching(SC_weak_matched,channels,g_weak);
