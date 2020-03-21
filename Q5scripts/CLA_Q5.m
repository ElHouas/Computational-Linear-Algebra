clear all;
%CLA Q5

%Number of grid points, as the matrix is square we don't need to define Ny
Nx=5;

%Create sparse matrix of A with size (2n-2)
[A,N,pos] = MatrixA2D(Nx); 

%Visualize the sparsity pattern of the matrix with the function spy()

figure(1)
spy(A);

%We compute the temperature vector using the function SolveMatrixA2D that
%receives the matrix A and the number of points N.
%The function returns the number the iterations needed and the timeelapsed
%This will be useful when to compare the result for different number of
%grid points

[T,timeElapsed,iterations] = SolveMatrixA2D(A,N);

%Create a Matrix M where to allocate the temperature values obtained from
%the function SolveMatrixA2D and plot it afterwards
M = zeros(N,N);

%This is a counter used to allocate the values from the column vecctor T
%in the right place of the matrix M
k = 1;

%Loop over the Matrix M to allocate the points properly
for j = 1:N
    for i=1:N
        
        %The matrix "pos", obtained from the function MatrixA2D, has 0 
        %values where we don't want to eliminate values and the cells
        %with values contains the index of the the Temperature value to
        %delete
        
        if(pos(i,j) == 0)    %If the value at the cell is 0, we want  
            M(i,j) = T(k,1); %to save this point in our matrix M to plot 
                             %it afterwards.
            k=k+1;           %We add 1 to the counter to go through the 
                             %values of T
        else
            M(i,j) = 0;      %If the value at the cell is not 0, we want
                             %to make this point 0 in our matrix M to 
                             %discard the points that are outside our
                             %domain and then get the L-shape.
        end    
    end
end


%We add the boundaries conditions as we know that the initial and final
%point are 0

Mi=0;
Mf=0;

%We add a row and column of zeros at the beginning and at the end of the
%matrix
M = [zeros(N,1) M zeros(N,1)];
M = [zeros(1,N+2);M; zeros(1,N+2)];


%Code to plot the surface and iscontours of the temperature distribution
figure(2)
surf(M)
title({'2-D Poisson''s equation surface';['{\itNumber of Grid Points} = ',num2str(Nx)];['{\itNumber of iterations} = ',num2str(iterations)]})
xlabel('Spatial co-ordinate (x) \rightarrow')
ylabel('{\leftarrow} Spatial co-ordinate (y)')
zlabel('Solution profile (T) \rightarrow')

figure(3)
contour(M)
title({'2-D Poisson''s equation isocontours';['{\itNumber of Grid Points} = ',num2str(Nx)];['{\itNumber of iterations} = ',num2str(iterations)]})
xlabel('Spatial co-ordinate (x) \rightarrow')
ylabel('{\leftarrow} Spatial co-ordinate (y)')
zlabel('Solution profile (T) \rightarrow')


%Code to plot the runt time versus the number of points

%time = [0.004, 0.011, 0.035];
iter =[64, 228, 795, 1184];
N=[5, 10, 20,25];
loglog(iter,N,'- .','Markersize',20,'MarkerEdgeColor','black');
title('Run time versus \Delta x in log scale')
xlabel('log(T) [s]')
ylabel('log(\epsilon_{max})')
