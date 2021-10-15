function [g_strong_Matched_hung]  = matchMMTC(users,channels,g_URLLC_matched,d_mMTC,p_mMTC,p_URLLC,L,a,Bs,No)

rng default;
%%% Construct cost matrix %%%

h_mMTC = zeros(users/2,channels);
g_mMTC = zeros(users/2,channels);
Energy_mMTC = zeros(users/2,channels);
for k=1:users/2
    for s=1:channels
        h_mMTC(k,s) = exprnd(2);
    end
    g_mMTC(k,:) = ((d_mMTC(k).^(-a)).*(abs(h_mMTC(k,:)).^2))./(Bs*No);
    Energy_mMTC(k,:) = (p_mMTC*L)./(Bs*log2(1+(p_mMTC*g_mMTC(k,:))./(p_URLLC*g_URLLC_matched(k) + 1)));
end


%%% Μatching mMTC users with subchannels/URLLC users %%%

% Hungarian Methodd
hung_matching = matchpairs(Energy_mMTC,1);
hung_matching = hung_matching(:,1);
g_strong_Matched_hung = channelGain_afterMatching(hung_matching,channels,g_mMTC);

% % Random Method
% random_matching = randperm(channels);
% g_strong_Matched_rand = channelGain_afterMatching(random_matching,channels,g_mMTC);

% % Exhaustive Method
% exhaustive_matching = exh_match(Energy_mMTC);
% g_strong_Matched_ex = channelGain_afterMatching(exhaustive_matching,channels,g_mMTC);




