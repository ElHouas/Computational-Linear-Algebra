function [dL_K,dL_M,S_K,S_M] = SENS(A,dK,dM)

[m,n] = size(A);
        if m~=n
        	disp('A is not a square matrix') ;
        	return;
        end

%Compute the eigenvectors and eigenvalues of our matrix to use them later
%in the sensitivity analysis

[V,D]=eig(A);

%Initialize the sensitivity matrix that or K and M that has the same size 
%as A

S_K=zeros(m,n);
S_M=zeros(m,n);

%Initialize the loop to go through the whole points and construct the
%sensitivity matrix.

%The S matrix is build by S=[dL_1/da_1.....dL_n/da_n], every row is formed
%by the multiplication of X'*dK/da*X. So in the first row, we will have
%the elements resutl of multiplying the first eigenvector by the
%derivatives of every K or M matrix.
                                      
for i=1:n
    for j=1:n
        S_K(i,j) = V(:,i)'*dK(:,:,j)*V(:,i);
        S_M(i,j) = V(:,i)'*(-D(i,i)*dM(:,:,j))*V(:,i);
    end
end

%Define an error in the Young Modulus E and the rho value, assigning to 
%DeltaE and DeltaRho to check the effect in the increment of eigenvalues
%DeltaL

%The problem to be solved is DL=S*DE

%Initialize vector of input errors that will multiply the
%sensitivity matrices
D_E = zeros(n,1);
D_Rho = zeros(n,1);

%Define the magintude of the error that we want to study   
e_E = 10^-2;
e_Rho = 10^-2;


%Define the vector of errors for E and rho
D_E(:,1) = e_E;
D_Rho(:,1) = e_Rho;


%Compute the value of the deviation experienced in the eigenvalues
%for the stifness and mass matrix   
dL_K= S_K*D_E;
dL_M= S_M*D_Rho;


end

