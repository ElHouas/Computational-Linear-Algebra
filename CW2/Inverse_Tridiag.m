function[A_inv] =  Inverse_Tridiag(A)

[m,n] = size(A);
        if m~=n
        	disp('A is not a square matrix') ;
        	return;
        end

%This function computes the inverse of the matrix A using the principal
%minors

%Initialize an empty inverse matrix of the same size
A_inv=zeros(m,n);

a = diag(A);
b = diag(A,1);
c = diag(A,-1);

%Define vector of prinipal minors
theta = zeros(1,n+1);
%Define parameter phi to use in the recurrence formula
phi = zeros(1,n);

%Initialize initial conditions or theta and phi
theta(1) = 1;
theta(2) = a(1);
phi(n+1) = 1;
phi(n) = a(n);

%Loop to compute principal minors, where we need the first two initial
%values to compute the third one and so on so for
for i=3:n+1
    theta(i) = a(i-1)*theta(i-1)-b(i-2)*c(i-2)*theta(i-2);
end

%Backwards loop to compute parameter phi
for i=n-1:-1:1
    phi(i) = a(i)*phi(i+1)-b(i)*c(i)*phi(i+2);
end

%Main loop to construct the inverse matrix of A.
%During the loop the algoithm states to be aware of three possible
%solutions and compute them seprately.
for i=1:n
    for j=1:n
        if (i<j)
            A_inv(i,j) = (-1)^(i+j)*prod(b(i:j-1))*theta(i)*(phi(j+1)/theta(n+1));
        elseif (i==j)
            A_inv(i,j) = theta(i)*(phi(j+1)/theta(n+1));
        else
            A_inv(i,j) = (-1)^(i+j)*prod(c(j:i-1))*theta(j)*(phi(i+1)/theta(n+1));
        end 
    end
end

end