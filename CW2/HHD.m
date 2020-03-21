function[H_eig_values,H_eig_vector,H_timeElapsed] =  HHD(A,v,lambda,tol)
tic ;   
%Get size of Matrix A
    [m,n] = size(A);
        if m~=n
        	disp('A is not a square matrix') ;
        	return;
        end
        
%Apply the Householder method in order to obtain the first 
%Hermitanian Matrix, before to start the loop 

        A_prev=A; %Store A value to use it later in the inverse iteration
        H = eye(n)-2*(v*v');
        A = H*A*H';
        B = A(2:m,2:n);
                
        HH = zeros(n,n,n); %3D array which stores the Hermitanian matrices
        HH(:,:,1) = H; %First H matrix computed and store it
        H_eig_values=zeros(n,1);
        H_eig_values(1)=lambda;
        A=B;

%Apply the Householder method from the second point until the size of
%the matrix 

for i=2:m
    x0 = ones(m+1-i,1);
    
    %Apply the Rayleigh Quotient Iteration method to obtain eigenvalues 
    %and eigenvectors
    [x,lambda,iterations] = RQI(A,tol);
      
    %Compute the perpendicular vector to the hyperplane with the new 
    %eigenvectors x
    [v] =  vector_perpendicular(x,m,i);
    
    %Compute another time the Hermitanian matrix with the new perpendicular
    %vector
    H = eye(size(A))-2*(v*v');
    A = H*A*H';
    B = A(2:m+1-i,2:n+1-i);  
    
    %Compute each H matrix and store them in the 3D array HH. 
    %At each iteration, a smaller matrix is being store it, so a 
    % mathematical ralation rule has been used in the indeces to be able 
    % to store them properly
    
    HH(1:m+1-i,1:n+1-i,i) = H;
    
    %Store the eigenvectos and eigenvalues obtained
    %H_eig_vector(1:m,i)=x;
    H_eig_values(i)=lambda;
    
    %Reasign the D matrix to A to continue the loop computations
    A=B;
end

%Instead of computing the eigenvectors inside the Householder loop, an
%inverse iteration is performed afterwards once the eigenvalues are
%obtained in order to obtain the eigenvectors faster

%Initialize variables of vectors to compute later the eigenvectors 
V=zeros(n,n);
H_eig_vector=ones(n,n);

for i = 1:n
    error=1;
    %While loop until convergence is achieved
    while error>10^-4
        %Compute the eigenvectors with shifting using the eigenvalues
        %previously computed
        V(:,i) = (A_prev-H_eig_values(i)* eye(n,n))\H_eig_vector(:,i);
        %Normalize the obtained eigenvectors
        V(:,i) = V(:,i)/norm(V(:,i));
        %Compute the error between the actual eigenvectors and the previous
        %ones to check if convergenced has been achieved
        error = max(abs(H_eig_vector(:,i))-abs(V(:,i)));
        H_eig_vector(:,i) = V(:,i);
    end   
end

H_timeElapsed = toc;
end