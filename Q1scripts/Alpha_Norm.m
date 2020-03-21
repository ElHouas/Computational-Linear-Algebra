clear all;
%Define size of matrix A
n=8;

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

L_abs=abs(L);%We define |L|

U_abs=abs(U);%We define |U|

R=L_abs*U_abs; %We compute the product |L|.|U|

%We create the vector which will contain the value of the sum of each row
S = zeros(1,n);

%We sum the coefficents row by row

for i=1:n
    v=0;
        for j=1:n
            v=v+R(i,j);
        end
        
        %The value of the sum of the ith row is saved as the ith value of c
S(i)=v;
end

%We get the highest value at S_max
S_max=max(S);


A_abs=abs(A);%We compute its absolute value

%We create the Sum vector which will contain the value of the sum of 
%each row

T = zeros(1,n);
%We sum the coefficents row by row

for i=1:n
    c=0;
        for j=1:n
            c=c+A_abs(i,j);
        end
        
        %The value of the sum of the ith row is saved as the ith value of c
T(i)=c;
end

Norm_max = max(T); %n is the result of the matrix norm for A
alpha = S_max/Norm_max; %alpha is the formula set by the coursework statement