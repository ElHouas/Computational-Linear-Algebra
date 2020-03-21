function [I] = Partial_pivoting(A,k,n)

%This script compute the matrix I (in general P, the permutation matrix)
%which swaps rows.
%'a' and 'b' are the positions of the rows that we want to swap.

%First, we define I as the identity matrix.
I=eye(n);

%We search the pivot for the column j
a=k;
max=0;

for j=k:n
    if abs(A(j,k))>max
        max=abs(A(j,k));
        index=j;
    end
end

b=index;

%The justification of this process is explained with more details
%in the report.

I(a,a)=0;%the coefficient 'Iaa' and 'Ibb' are now 0.
I(a,b)=1;%'Iab' is now 1.
I(b,b)=0;
I(b,a)=1;%And 'Iba'=1.

%We multiply the matrix for which we want to swap rows with I
%and get A with rows 'a' and 'b' switched.
end