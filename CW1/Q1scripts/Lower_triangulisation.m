function L = Lower_triangulisation(A,n)

% The function computes the lower triangularisation of a given
%square matrix A using the Gauss elimination method (GEM).
% A - the coefficient matrix
% n - the size of the matrix A
%Initialisation of the identity (I) and upper-triangular (U) matrices

for i=1:n
    for j=1:n
        if (i==j)
            I(i,j)=1.0;
            U(i,j) = A(i,j);
        else
            I(i,j) = 0.0;
            U(i,j) = A(i,j);
        end
    end
end

%The Lower triangularisation matrix needs to be computed by adding some 
%more steps to the uppertriangularisation code.

L=I;  %Initialize L as the identity matrix

for k=1:n-1 % the main loop of the GEM starts here

% First step in GEM is to construct Gauss vector
gv = gauss_vector(U,k,n);

%Secondly, I construct e vector for the corresponding stage of the GEM
evec = e_vector(k,n);

%Here, I take outer product between the Gauss and e vectors for
%the given stage
LE = outer_product(gv,evec,n);

%Here, I construct the Gauss transformation matrix by subtracting
%the previously computed outer product from the identity matrix
J = matrix_subtraction(I,LE,n);

%Here, I compute the inverse of the matrix J in order to calculate L
J_inv=inv(J);

%We compute the updated upper-triangular matrix through matrix
%multiplication with the Gauss transformation matrix.
%It will be used in the loop to calculate J and therefore: L
U = matrix_matrix_mult(J,U,n);

%According to lecture 4, L=inv(J1)*...*inv(k), with k?n-1
L=matrix_matrix_mult(L,J_inv,n);

end % the main loop of the GEM ends here
end