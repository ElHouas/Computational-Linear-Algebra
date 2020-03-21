clear all
%%CLA_Q1

%Declare the matrix of the system Ax=b given in the coursework

A=[1 0 0 0 1; -1 1 0 0 1; -1 -1 1 0 1; -1 -1 -1 1 1; -1 -1 -1 -1 1];
b=[2;1;0;-1;-3];

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

%We now calculate LY=Pb
y = Forward_substitution(L,P*b,n);

%Then, Ux=Y
x = Backward_substitution(U,y,n);

%And we finally get the unknown vector X