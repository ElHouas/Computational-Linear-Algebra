function [P,D,X] = RITZ(A,B,X0,tol)

 %Get size of Matrix A
    [m,n] = size(A);
        if m~=n
        	disp('A is not a square matrix') ;
        	return;
        end

%As the matrices are tridiagonal, they can be stored in sparse format to
%take adavantage of theiR shape and only compute with the nonzero entries.
%Storing the matrices in this format will speed up the computations.

A=sparse(A);
B=sparse(B);

%The matrices from the eigenvalue problem
%K_*P=M_*P*Delta are defined

[B_inv] =  Inverse_Tridiag(B);

X = B_inv*A*X0;
A_bar = X'*A*X;
B_bar = X'*B*X;
              
% The product matrix M^{-1}*K is computed

L=A_bar\B_bar;

%P and D are computed using QR iteration which contains 
%the eigenvectors and eigenvalues respectively
[P,D,~] = QRITER(L,tol);


end

