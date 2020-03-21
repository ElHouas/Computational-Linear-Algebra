%%Matrix1D function to compute sparse matrix B  with a defined number of
%%points n
function [B] = Matrix1D(n)
%The function receives the n points of the domain. As it is a square
%matrix, we have the same number of points in x and y diretion

%We create a column vector of ones of length n-1 in order to be able to
%define the diagonals of the matrix
Ones = ones(n-2, 1);

%Define tridiagonal matrix of size (n-1)X(n-1) with the main diagonal of 2 and the upper and
%lower diagonal of -1

%We use the Matlab function diag, which allows us to create a matrix with
%a diagonal that where we can specify the length and the position. 
%Then, combining the three matrices we obtain the tridiagonal matrix.

B = diag(2 * Ones, 0) - diag(Ones(1:n-3), -1) - diag(Ones(1:n-3), 1);

%Using Matlab sparse function to store pnly the non-zero values of the tridiagonal
%matrix
B = sparse(B);
end