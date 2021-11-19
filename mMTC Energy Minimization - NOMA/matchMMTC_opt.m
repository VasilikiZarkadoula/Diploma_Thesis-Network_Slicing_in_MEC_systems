function [matching,g_mMTC_Matched]  = matchMMTC_opt(users,h_mMTC,g_URLLC_matched,d_mMTC,p_mMTC,p_URLLC,L,a,Bs,No,method)
% This function returns the optimal assignment of mMTCs to subchannels and
% their channel gains after matching

channels = users/2;

%%% Construct cost matrix %%%
g_mMTC = zeros(users/2,channels);
Energy_mMTC = zeros(users/2,channels);
for k=1:users/2
    g_mMTC(k,:) = ((d_mMTC(k).^(-a)).*(abs(h_mMTC(k,:)).^2))./(Bs*No);
    Energy_mMTC(k,:) = (p_mMTC(k)*L)./(Bs*log2(1+(p_mMTC(k)*g_mMTC(k,:))./(p_URLLC(k)*g_URLLC_matched(k) + 1)));
end

%%% Μatching mMTC users with subchannels/URLLC users %%%

if method  == 0 % Hungarian Algorithm
    
    hung_matching = matchpairs(Energy_mMTC,1);
    matching = hung_matching(:,1);
    g_mMTC_Matched = channelGain_afterMatching(matching,channels,g_mMTC);
    
elseif method  == 1 % Exhaustive search
    
    matching = exh_match(Energy_mMTC);
    g_mMTC_Matched = channelGain_afterMatching(matching,channels,g_mMTC);
    
else
    
    disp('\nPlease type a valid matching method\n')
    
end
