clear all;
close all;

%% Computational Linear Algebra CW2
%% Variables initialization to use the algebraic methods
% Given matrices by Structural Dynamics course for N=4

M=[3.2525 0.6486 0 0;
   0.6486 2.0868 0.4151 0
   0 0.4151 1.3320 0.2642
   0 0 0.2642 0.4709];

K=[9.4955 -4.0024 0 0;
   -4.0024 6.9048 -2.9024 0
   0 -2.9024 4.9956 -2.0932
   0 0 -2.0932 2.0932];

%As the matrices are tridiagonal, they can be stored in sparse format to
%take adavantage of theiR shape and only compute with the nonzero entries.
%Storing the matrices in this format will speed up the computations.

M=sparse(M);
K=sparse(K);

%Own function to inverse a tridiagonal matrix and solve  the system
%Ax=B
[M_inv] =  Inverse_Tridiag(M);

%Using the inverse computed function to solve the system W=M^(⁻1)*K
W = M_inv*K;

%Get size of matrix M
[m,n]=size(W);

% Error tolerance to check convergence
tol=10^-4;

%Trial matrix used in SSI method
X0=rand(m); 
%% Rayleigh Quotient Iteration and Householder deflation method

% Use the Rayleigh Equation to get the eigenvalues and eiegenvectors 
[x,lambda,RQI_iterations,RQI_timeElapsed] = RQI(W,tol);

% Compute the vector perpendicular to the hyperplane used in the 
% Householder method to compute the first Hermitanian matrix 
[v] =  vector_perpendicular(x,m,1);

% Apply the Householder method and obtain the eigenvalues and
% eigenvectors which are stored in H_eig_values,H_eig_vector respectively

[H_eig_values,H_eig_vector,H_timeElapsed] =  HHD(W,v,lambda,tol);  

% Apply Qriteration with the W matrix
[QR_eig_vector,QR_eig_values,QR_timeElapsed] =  QRITER(W,tol);

% Apply the SSI method using M, K and a trial random matrix X0
[SSI_eig_vector,SSI_eig_values,SSI_timeElapsed] =  SSI(M,K,X0,tol);

%% Sensitivity analysis for each iterative method

% Introduce some noise in the matrix W=M^(⁻1)*K, to check how the 
% eigenvalues and eigenvectors vary for each method. Instead of introducing
% the noise directly to the result matrix W, it is more realistic and
% meaningful to introduce noise in M and K where data from measures are
% stored and that error could be introduced in the process of comopute them

noise = 10^-1;
M_noise = M+noise;
K_noise = K+noise;

[M_inv_noise] =  Inverse_Tridiag(M_noise);

%Using the inverse computed function to solve the system W=M^(⁻1)*K
W_noise = M_inv_noise*K_noise;

%Repeat each iterative method for the error matrix

% Use the Rayleigh Equation to get the eigenvalues and eiegenvectors 
[x,lambda,RQI_iterations,~] = RQI(W_noise,tol);

% Compute the vector perpendicular to the hyperplane used in the 
% Householder method to compute the first Hermitanian matrix 
[v] =  vector_perpendicular(x,m,1);

% Apply the Householder method and obtain the eigenvalues and
% eigenvectors which are stored in H_eig_values,H_eig_vector respectively

[H_eig_values_noise,H_eig_vector_noise,~] =  HHD(W_noise,v,lambda,tol);  

[QR_eig_vector_noise,QR_eig_values_noise,~] =  QRITER(W_noise,tol);

[SSI_eig_vector_noise,SSI_eig_values_noise,~] =  SSI(M_noise,K_noise,X0,tol);

Delta_HHD_eig_values=round(norm(H_eig_values-H_eig_values_noise),4);
Delta_HHD_eig_vector=round(norm(H_eig_vector-H_eig_vector_noise),4);

Delta_QRITER_eig_values=round(norm(QR_eig_values-QR_eig_values_noise),4);
Delta_QRITER_eig_vector=round(norm(QR_eig_vector-QR_eig_vector_noise),4);

Delta_SSI_eig_values=round(norm(SSI_eig_values-SSI_eig_values_noise),4);
Delta_SSI_eig_vector=round(norm(SSI_eig_vector-SSI_eig_vector_noise),4);

%% Computational cost and accuracy for each iterative method

% In this section the computational cost and time for each method is
% compared with the inbuild matlab function eig(A) which computes the
% eigenvalues and eigenvectors for a given matrix

% Compute eigenvectors and eigenvalues from Matlab function
[M_eig_vectors, M_eig_values] = eig(W);

Error_HHD_eig_values=abs(M_eig_values-H_eig_values_noise);
Error_HHD_eig_vector=abs(M_eig_vectors-H_eig_vector_noise);

Error_QRITER_eig_values=abs(M_eig_values-QR_eig_values_noise);
Error_QRITER_eig_vector=abs(M_eig_vectors-QR_eig_vector_noise);

Error_SSI_eig_values=abs(M_eig_values-SSI_eig_values_noise);
Error_SSI_eig_vector=abs(M_eig_vectors-SSI_eig_vector_noise);

%% Sensitivity matrix analysis

%Use the small matrices given to assemble a vector 3D that contains the
%derivatives of the global matri for each small matrix

[dK,dM] =  Assembly;

%Compute the sensitivity matrices and the increment in the eigenvalues
%due to an arbitrary error introduced in the Young Modulus or Density
[dL_K,dL_M,S_K,S_M] = SENS(A,dK,dM);