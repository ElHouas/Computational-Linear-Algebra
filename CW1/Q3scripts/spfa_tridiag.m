function A = spfa_tridiag(A)
%Returns the matrix U computed via the Cholesky factorisation
%of a symmetric definite positive matrix A, writing the values of U
%directly in the upper triangular part of A
%based on Dr Kevin Gouder script for AEM-ADV08 2016-2017

[m,n]=size(A);
if (m~=n || max(max(abs((A+A')/2-A)))~=0 )
        error('error: A is not symmetric');
end

%Use  function to compute only the non-zero values
nZ=find(A(1,1:n)~=0);

b=nZ(end)-1;
for j = 1:n
    s = 0.0;   
        for k = max(1,j-b):j-1
            t = A(k,j);
            t = t/A(k,k);
            A(k,j) = t;
            s = s+t*t;
        end
        
        q = min(j+b,n);
        s = A(j,j)-s;
        if s<0
              error('error: A is not definite positive');
        end
        A(j,j) = sqrt(s);
end

end
