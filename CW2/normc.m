function [B] = normc(B)
%Function that normalize the column of each matrix
B = sqrt(B.^2 ./ sum(B.^2)) .* sign(B);
end

