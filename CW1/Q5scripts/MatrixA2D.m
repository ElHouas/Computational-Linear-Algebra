%%MatrixA2D function to compute sparse matrix B  with a defined number of
%%points n
function [A,N,pos] = MatrixA2D(Nx)

%Define number of points and length and we substract two points due to 
%the boundaries which are zero
N=(2*(Nx)+1)-2;

%Define the longitude of the domain without the boundaries and define the 
%width space in x.

Lx = 2*Nx+1;
DeltaX = 2/Lx;

%Create vector column of ones of length N
Ones = ones(N,1);

%Define tridiagonal matrix of size Nx^2 X Ny^2 with the main diagonal of  
%4 and the upper and lower diagonal of -1

%We use the Matlab function diag, which allows us to create a matrix with
%a diagonal that where we can specify the length and the position. 
%Then, combining the 5 matrices we obtain the desired matrix.

A= diag((4)/(DeltaX^2) * Ones, 0) + diag((-1)/(DeltaX^2)*Ones(1:N-1), -1) + diag((-1)/(DeltaX^2)*Ones(1:N-1), 1);

A = kron(eye(N), A);

%Shift the other non zero diagonal Ny positions and the length of these  
%diagonals is Nx-1 X Ny

off_diag = ((-1)/(DeltaX^2))*ones((N-1)*(N),1);

A = A + diag(off_diag, N)+ diag(off_diag, -N) ;
 
%As the problem is defined over a L shaped domain, the points of the 
%square matrix that are out this domain have to be deleted from the matrix
%Delete these points now will make our computations faster 

%The zone to elimiate it's defined as a small in the top right corner 
%square starting at (n-1) to (2*n-1)

for j = (2*Nx-1):-1:(Nx) %It loops backwards from the top-right point to  
                         %the half of the domain in the x-direction
    
    for i = (2*Nx-1):-1:(Nx)%It loops backwards from the top-right point   
                           %to the half of the domain in the x-direction
        
         p = i+(j-1)*(2*Nx-1);%This equation gives us the position of the 
                             %points we want to eliminate from our matrix
        
            A(p,:) = [];%The rows and columns that have a dependence with
            A(:,p) = [];%the target points are eliminated
        
            pos(i,j)=p;
    end
end

%Now we are storing a smaller sparse matrix than if we apply the BC later
A = sparse(A);
end