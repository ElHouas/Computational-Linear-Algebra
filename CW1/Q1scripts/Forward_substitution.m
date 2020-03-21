function b = Forward_substitution(L,b,n)

%LECTURE 3

%This program compute solution "y" with the forward substitution 
% applying Ly=b with L as a lower triangular matrix.

b(1)=b(1)/L(1,1);

for i=2:n

%We start with the first row which contains only one value so there is 
%only one division.    
        
%We compute xi as a function of the values of xj, j<i, as they
% have been calculated in the previous steps. It is not the real
% value of xi yet.
                
b(i) = (b(i)-inner_product(L(i,1:i-1),b(1:i-1),i-1))/L(i,i);
   
%Forward substitution is an O(n^{2}) process

end

end