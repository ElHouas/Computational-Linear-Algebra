function[Q,R] =  QR_Decomposition(A)

%Get size of Matrix A
    [m,n] = size(A);
        if m~=n
        	disp('A is not a square matrix') ;
        	return;
        end
        
%Initialize matrices
     U = zeros(m,n);
     e = zeros(m,n);
     Q = zeros(m,n);

%Compute the first element for matrices U, e and Q
     U(:,1) = A(:,1);
     e(:,1) = U(:,1)/norm(U(:,1)); 
     Q(:,1) = e(:,1);
     
for k=2:n
    %Assign A values to matrix U
    U(:,k) = A(:,k);
    s=1; 
    
    %In order to perform the right concatenate substraction
    %a while loop is need to substract 
    %u_n = a_n-...-proj_e_(n-1)*a_n
	while (s~=k)
        t = k-s;
        %Perform concatenate substraction using own made projection
        %function
        U(:,k) = U(:,k)-projection(A(:,k),e(:,t)); 
        s=s+1;
    end
    
    %Normalize the e vectors that are each column of U
    e(:,k) = U(:,k)/norm(U(:,k));
    
    %Assign each e vector to a column of Q
    Q(:,k) = e(:,k);
end

%Compute R with the Q obtained and the initial matrix A
R = Q'*A;

end
