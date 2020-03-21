function[x, lambda, RQI_iterations,timeElapsed] =  RQI(A,tol)
tic;
%Get size of Matrix A
[m,n] = size(A);
        if m~=n
        	disp('A is not a square matrix') ;
        	return;
        end

% Initialize vector where the eigenvectors are gonna be stored 
x = ones(m,1);

%Initialize eigenvalues
lambda = 1;
k=1; %Counter

    while norm( (A-lambda * eye(m,n))*x )  > tol
        lambda = (x'*A*x)/(x'*x); %Compute the Rayleigh quotient
        y = (A-lambda * eye(m,n))\x; %Shifting
        x = y/norm(y); %Normalize the eigenvectors
        k=k+1;
    end
RQI_iterations = k;
timeElapsed = toc;
end