function C = matrix_subtraction(A,B,n)
% This function computes the subtraction between two matrices
% A - The first input matrix
% B - The second input matrix
% n - The size of the matrices A and B
for i=1:n
for j=1:n
C(i,j) = A(i,j)-B(i,j); %a simple matrix subtraction procedure
end
end
