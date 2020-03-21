function[QR_eig_vector,QR_eig_values,QR_timeElapsed] =  QRITER(A,tol)
tic;
%Get size of Matrix A
    [m,n] = size(A);
        if m~=n
        	disp('A is not a square matrix') ;
        	return;
        end
    
    %Initialize matrix Qv
    Qv = eye(m,n);
    a=0;%Break loop condition

    %Convergence criteria that checks if the norm of the subdiagonal of  
    %the matrix A is less than the tolerance entered by the user.
    %As the elements of the subdiagonal are the largest ones, check  
    %if they are under the tolerance is enough to know if the method has
    %converged
    
    %Perform QR decomposition before entering the loop
    [Q,R]=QR_Decomposition(A);
    
    while  (a==0)
        %Compute A with Q and R to perform another time QR_Decomposition
        A_new=R*Q;
        
        %Multiply Q by the matrix Qv
        Qv=Qv*Q;
        
         if (norm((diag(A_new,-1))) < tol)
            a=1;
         end
        
        %Perform QR_Decomposition
        [Q,R]=QR_Decomposition(A_new);
        
    end      

    %In the case that A is not symmetric, the eigenvectos are not yet 
    %found. So, a backwards substitution is made in order to find the right
    %eigenvectors
    
    if (issymmetric(A)==0)
        Rv=zeros(n);
        l = diag(R);
        for j=n:-1:1
            for i = n:-1:1
                if (i==j)
                    Rv(i,j) = 1;
                 elseif (j<i)
                     Rv(i,j) = 0;
                else
                    Rv(i,j) = dot(R(i,i:n),Rv(i:n,j))/(l(j)-R(i,i));
                end
            end
            Rv(:,j)=Rv(:,j)/norm(Rv(:,j));
        end
        Qv=Qv*Rv;
    end

    %The eigenvectors are the columns of matrix 
    QR_eig_vector = Qv;
    
    %The eigenvalues are on the diagonal of R
    QR_eig_values = diag(diag(R));
    QR_timeElapsed = toc;
end