function [random_matching,g_mMTC_Matched_rand]  = matchMMTC_rand(users,h_mMTC,d_mMTC,a,Bs,No)
% This function returns the random assignment of mMTCs to subchannels and
% their channel gains after matching
rng default;

channels = users/2;
% Construct cost matrix 
g_mMTC = zeros(users/2,channels);
for k=1:users/2
    g_mMTC(k,:) = ((d_mMTC(k).^(-a)).*(abs(h_mMTC(k,:)).^2))./(Bs*No);
end

% Μatching mMTC users with subchannels/URLLC users - Random Method
random_matching = randperm(channels);
g_mMTC_Matched_rand = channelGain_afterMatching(random_matching,channels,g_mMTC);
