function exhaustive_matching = exh_match(cost_matrix)
% Exhaustive Search Algorithm

rng shuffle;

% Create a list of possibile indicies in the n x n matrix
index_range = 1:length(cost_matrix);

% Create a permutation of array indicies
permutation = perms(index_range);

lowest = -1;
temp = zeros(1,length(cost_matrix));          % Temporary tuple storage

% Loop through all permutations possible
for i = 1:length(permutation)
    
    index_assignment = permutation(i,:);
    cost = 0;
    
    % Loop through all index assignments
    first_perm = permutation(1,:);
    for j = 1:length(first_perm)
        user = first_perm(j);
        cost = cost + cost_matrix(user,index_assignment(user));
        
    end
    if cost < lowest || lowest == -1
        lowest = cost;
        temp = index_assignment;
    end
end
us = index_range;
arr = [us' temp'];
arr = flip(arr');
out = sortrows(arr.',1).';
exhaustive_matching = out(2,:);


