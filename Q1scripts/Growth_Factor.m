clear all;

%Define size of matrix A
n=32;

%We create A as a function of 'n'

A = eye(n);
A = A+tril(-ones(n,n),-1);

for i=1:n
    A(i,n) = 1; 
end

[~,n]=size(A);

Pivf=eye(n);        %Initialize the permutation matrix as identity matrix
L=zeros(n);         %Initialize the lower triangular matrix
U=A;                %Initialize A matrix as U matrix before modifying it 
I=eye(n);           %Initialize the identity matrix

for k=1:n
    
P=Partial_pivoting(U,k,n); %Perform partial pivoting to obtain 
                           %the permutation matrix
U=P*U;                     %Compute the U matrix with pivoting matrix
L=P*L;                     %Compute the L matrix with pivoting matrix
Pivf = P*Pivf;
gv = gauss_vector(U,k,n); % the first step in GEM is to construct Gauss 
                          %vector
evec = e_vector(k,n); %construct e vector for the corresponding stage 
                      %of the GEM
LE = outer_product(gv,evec,n); %Here, I take outer product between the 
                               %Gauss and e vectors for the given stage
J = matrix_subtraction(I,LE,n); %Here, I construct my Gauss transformation 
                           %matrix by subtracting the previously 
                           %computed outer product from the identity matrix

U = matrix_matrix_mult(J,U,n); %Here, I compute the updated 
                    %upper-triangular matrix through matrix multiplication  
                    %with the Gauss transformation matrix.
L=L+LE;
end

%L is the sum of l and e vectors
L=L+eye(n);

%We calculate the maximum absolute value of U
U_max = max(max(abs(U)));
%We calculate the maximum absolute value of A
A_max = max(max(abs(A)));
%Thus, we get rho according to the formula
rho = U_max/A_max;