function [SSI_eig_vector,SSI_eig_values,SSI_timeElapsed] = SSI(A,B,X0,tol)
tic;   
    %Get size of Matrix A
    [m,n] = size(A);
        if m~=n
        	disp('A is not a square matrix') ;
        	return;
        end

    a=0;%Variable used to break the loop in case of convergence
    Dprev=0;
    k=1; %Counter
    
    while (a==0)
        %Solve the Ritz eigenvalue problem and obtain
        %the matrix P of eigenvectors the D of eigenvalues and Xk+1
        [P,D,X] = RITZ(A,B,X0,tol);
         q = X*P;
         X0 = q;
         q = normc(q);
         if(norm(diag(D)-diag(Dprev)) < tol)
            a=1;
         end
         Dprev = D;
         k=k+1;
    end
       
    %The eigenvectors are stored in q and assigned to SSI_eig_vector 
    SSI_eig_vector = q;
    
    %Normalize the eigenvectors obtained 
    [SSI_eig_vector] = normc(SSI_eig_vector);
    
    %The eigenvalues are stored in D and assigned to SSI_eig_values 
    SSI_eig_values = D;
    SSI_timeElapsed = toc;
end
