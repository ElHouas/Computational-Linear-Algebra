function b = Backward_substitution(U,b,n)

%LECTURE 3

%This program compute solution "y" with the forward substitution applying
% Ux=b with U as a upper triangular matrix.


%We compute the last element first, which is the n element to avoid 
%singularity
b(n) = b(n)/U(n,n);

%Loop backwards to go from the last row 
for i=n-1:-1:1

%We compute xi as a function of the values of xj, j>i, as
%they have been calculated in the previous steps. It is not 
%the real value of xi yet.
    
b(i) = (b(i)-inner_product(U(i,i+1:n),b(i+1:n),n-i))/U(i,i);
    
end

end