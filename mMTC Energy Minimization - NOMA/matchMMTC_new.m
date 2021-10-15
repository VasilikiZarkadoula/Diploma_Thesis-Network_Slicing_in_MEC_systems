function [hung_matching,g_mMTC_Matched_hung]  = matchMMTC_new(users,channels,h_mMTC,g_URLLC_matched,d_mMTC,p_mMTC,p_URLLC,L,a,Bs,No)

rng default;
%%% Construct cost matrix %%%

g_mMTC = zeros(users/2,channels);
Energy_mMTC = zeros(users/2,channels);
for k=1:users/2
    g_mMTC(k,:) = ((d_mMTC(k).^(-a)).*(abs(h_mMTC(k,:)).^2))./(Bs*No);
    Energy_mMTC(k,:) = (p_mMTC(k)*L)./(Bs*log2(1+(p_mMTC(k)*g_mMTC(k,:))./(p_URLLC(k)*g_URLLC_matched(k) + 1)));
end


%%% Μatching mMTC users with subchannels/URLLC users %%%

% Hungarian Methodd
hung_matching = matchpairs(Energy_mMTC,1);
hung_matching = hung_matching(:,1);
g_mMTC_Matched_hung = channelGain_afterMatching(hung_matching,channels,g_mMTC);
