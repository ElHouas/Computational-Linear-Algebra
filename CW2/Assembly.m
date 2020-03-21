function [dK,dM] =  Assembly

%Initialize the variables of Young Modulus and Density
Eo = 4;
rho = 3;

%Local stifness matrices to build K
k1=[5.4931 -5.4931;
    -5.4931 5.4931];
k1=k1/Eo;

k2=[4.0024 -4.0024;
    -4.0024 4.0024];
k2=k2/Eo;

k3=[2.9024 -2.9024;
    -2.9024 2.9024];
k3=k3/Eo;

k4=[2.0932 -2.0932;
    -2.0932 2.0932];
k4=k4/Eo;

%Local mass matrices to build M 

m1=[2.2484 1.0087;
    1.0087 1.8052];
m1=m1/rho;

m2=[1.4473 0.6486;
    0.6486 1.1594];
m2=m2/rho;

m3=[0.9273 0.4151;
    0.4151 0.7410];
m3=m3/rho;

m4=[0.5910 0.2642;
     0.2642 0.4709];
m4=m4/rho;

%Initilalize empty matrices for K and M, and use them as tool to assembly
%the derivatives for each local matrix

K_global = zeros(4, 4);
M_global = zeros(4, 4);

%Get the size of the K global matirx
[m,n] = size(K_global);

%Create each derivative matrix of the local stifness matrices, assigning
%the right indices to the position specified in the report
dk1 = zeros(4, 4);
dk1(1,1)=K_global(1,1)+k1(2,2);

dk2 = zeros(4, 4);
dk2(1:2,1:2)=K_global(1:2, 1:2)+k2;

dk3 = zeros(4, 4);
dk3(2:3,2:3)=K_global(2:3,2:3)+k3;

dk4 = zeros(4, 4);
dk4(3:4,3:4)=K_global(3:4,3:4)+k4;

%Initialize 3D array to store each derivative matrix
dK = zeros(n,n,n);

%Assign each derivative to a position of the 3D array
dK(:,:,1)=dk1;
dK(:,:,2)=dk2;
dK(:,:,3)=dk3;
dK(:,:,4)=dk4;

%Create each derivative matrix of the local mass matrices, assigning
%the right indices to the position specified in the report
dm1 = zeros(4, 4);
dm1(1,1)=M_global(1,1)+m1(2,2);

dm2 = zeros(4, 4);
dm2(1:2,1:2)=M_global(1:2, 1:2)+m2;

dm3 = zeros(4, 4);
dm3(2:3,2:3)=M_global(2:3,2:3)+m3;

dm4 = zeros(4, 4);
dm4(3:4,3:4)=M_global(3:4,3:4)+m4;

%Initialize 3D array to store each derivative matrix
dM = zeros(n,n,n);

%Assign each derivative to a position of the 3D array
dM(:,:,1)=dm1;
dM(:,:,2)=dm2;
dM(:,:,3)=dm3;
dM(:,:,4)=dm4;

end
