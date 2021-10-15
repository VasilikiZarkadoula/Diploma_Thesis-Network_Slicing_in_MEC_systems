function [f,g] = obj_scaling(x)
% objective f
f = x(1);

if nargout > 1 % gradient 
    g = [1; zeros(length(x)-1,1)];
 
end

