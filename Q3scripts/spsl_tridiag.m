function b = spsl_tridiag(U,b)
%Solve Ax=b for symmetric definite positive matrices after 
%the the upper triangular matrix U was computed by spfa.m
%based on Dr Kevin Gouder script for AEM-ADV08 2016-2017
[~,n]=size(U);

%Solve U'y=b the values of y are written directly in b to save memory;
b(1)=b(1)/(U(1,1));
for k = 2:n
    if(k-1 ~= 0)
        t = dot(U(k-1,k),b(k-1));
        b(k) = (b(k)-t)/U(k,k);
    end
end

%Solve Ux=y, 
%to save memory the values of x are written inside the array b.
for k = n:-1:2
    if(k-1 ~= 0)
        b(k) = b(k)/U(k,k);
        b(k-1) = b(k-1)-b(k)*U(k-1,k);
    end
end
b(1)=b(1)/(U(1,1));


end 

