function[T,timeElapsed,iterations] = SolveMatrixA2D(A,Nx)

%A time counter will be used to determine the needed time for the iterative
%solver to converge. It is used the tic-toc Matlab feature.

tic; %The time counter is activitated

%Get matrix size, by obtaining the number of points in x and y
[ny,nx]=size(A);


q = 1;                      %Source term
b = -q*ones(ny,1);          %Value resulted for multiplying q and DeltaX
T = zeros(ny,1);            %Initialize the Temperature columb vector
                            %with zeros

%The problem is going to use an iterative solver which checks that
%the convergence of the solution is under a certain tolerance, 
%defined by the user

tol = 10^-6;                  %Error tolerance used for check convergence

%http://employees.oneonta.edu/GoutziCJ/fall_1999/math323/matlab/lesson_05.pdf

%The itertive solver will use the LDU decomposition, which allows to 
%solve the problem Ax=b, using Gauss-Seidel forward substitution,
%by decomposing A=L+D+U. Then, we have the following eq:
%(L+D)x^(k+1) = b-Ux^k)

%First, the lower and upper triangular matrices and a main diagonal are 
%obtained matrix uisng defined Matlab functions

L = tril(A,-1);     %Lower triangular matrix using Matlab function tril()
D = diag(diag(A));  %Main diagonal matrix using Matlab function diag()
U = triu(A,1);      %Upper triangular matrixusing Matlab function triu()

Tn = zeros(ny,1);            %Preallocating Tn
a=0; %Condition to break loop
m=1; %Counter of the matrices computed inside the loop

while (a~=1)
    Tn=T;
    T = (L+D)\(b-U*Tn);
    %Error betwen p and pn
        dif = abs(T - Tn);
        eps= dif;
        
        for i=2:ny-1
            eps(i,1) = abs(dif(i,1)/T(i,1));
        end
    
        e_max(m) = max(eps(:));
        
        if(e_max(m) < tol)
            a=1;
        end
        m=m+1;
        
end
       
     timeElapsed = toc;
     iterations = m;
    
end