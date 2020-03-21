clear all
%% Question 3
% Heat equation using Crank-Nicholson discretization 

%% Compute the Concentration matrix for the analytical solution
%Define the constant values N,D,L,Co,t,n

N = 40;            %Number of points in the domain
D = 1;             %Diffusivity constant
L = 1;             %Length of the domain
Co = 1;            %Initial Concentration
t = 0.1;           %Discretization time
n = 100;           %Number of time steps
DeltaT = t/n;      %Increment of time
DeltaX = L/(N-1);  %Increment of space
S = 0;             %Initialazation of a variable to Summation term

%Loop over the Fourier series 
for j = 1:n %We loop over each time step
    T = DeltaT*j;
    for i = 1:N-2 
        x = DeltaX*i;
        for m = 1:n  %Compute sumation terms from m=1 truncating at m=100 
            %We assign each term of the summation factor to a variable, in
            %order to be have a clear code
            Y = (D*((2*m)-1)^2*(pi^2)*T/(L^2)); 
            P = (((2*m)-1)*pi*x)/(L);
            Q = 1/(((2*m)-1)*pi);
            S = S + (Q*exp(-Y)*sin(P));
            Y = 0;
        end
        C_a(i,j) = 4*Co*S;
        S = 0;
    end
end

%Add boundary conditions in the first and last point, 
%which are specified as zero

Ci=0;
Cf=0;
C_a = [Ci*ones(1,n);C_a;Cf*ones(1,n)];

X = linspace(0,1,N);  %Space vector to plot the specified number of points

figure (1)
hold on

for j=10:16:100 %Plot concentration distribution for different time steps
 T = DeltaT*j;
 time = ['t = ', num2str(T)];
 plot(X, C_a(:,j), 'DisplayName',time)
 title('Concentration distribution for analytical solution N=40');
 ylabel('Concentration');
 xlabel('Distance');
end
hold off
legend show

%% Compute the Concentration matrix for the numerical solution

B = Matrix1D(N);                    % Sparse Matrix using  Matrix1D.m

sigma = D * (DeltaT/((DeltaX)^2));  %Sigma value otbtained from DT and DX
I = eye(size(B,1));                 %Identity matrix of the size of B
Z1 = (I +(sigma/2) * B);            %Matrix multipliying C^(n+1)
Z2 = (I -(sigma/2) * B);            %Matrix multipliying C^(n)
C_n=Co*ones(N-2,n);                 %Initialize concentration values at 0

U = spfa_tridiag(Z1);               %Obtain U matrix to solve Ax=b problem
                                    %It is computed using spfa_tridiag.m

                           
for j = 1:n-1                       %Loop until number of time steps=100
    b = Z2*C_n(:,j);                %Compute column vecctor solution b
    C_n(:,j+1) = spsl_tridiag(U, b);%Compute the new concentration solving   
end                                 %Ax=b using spsl_triridiag.m
                            
%Add boundary conditions in the first and last point, 
%which are specified as zero

Ci = 0;
Cf = 0;
C_n = [Ci*ones(1,n);C_n;Cf*ones(1,n)];

figure (2)
hold on

for j=10:15:100 %Plot concentration distribution for different time steps
 T = DeltaT*j;
 time = ['t = ', num2str(T)];
 plot(X, C_n(:,j), 'DisplayName',time)
 title('Concentration distribution for numerical solution N=40');
 ylabel('Concentration');
 xlabel('Distance');
end
hold off
legend show

%% Plots

figure (3)

subplot(3,2,1); plot(X, C_a(:,10),'b',X,C_n(:,10),'r'), title('t=0.01s'),axis([0 1 0 1]),xticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]), yticks([0 0.2 0.4 0.6 0.8 1])
ylabel('Concentration');
xlabel('Distance');

subplot(3,2,2); plot(X, C_a(:,30),'b',X,C_n(:,30),'r'), title('t=0.03s'),axis([0 1 0 1]),xticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]), yticks([0 0.2 0.4 0.6 0.8 1])
ylabel('Concentration');
xlabel('Distance');

subplot(3,2,3); plot(X, C_a(:,50),'b',X,C_n(:,50),'r'), title('t=0.05s'),axis([0 1 0 1]),xticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]), yticks([0 0.2 0.4 0.6 0.8 1])
ylabel('Concentration');
xlabel('Distance');

subplot(3,2,4); plot(X, C_a(:,70),'b',X,C_n(:,70),'r'), title('t=0.07s'),axis([0 1 0 1]),xticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]), yticks([0 0.2 0.4 0.6 0.8 1])
ylabel('Concentration');
xlabel('Distance');

subplot(3,2,5); plot(X, C_a(:,90),'b',X,C_n(:,90),'r'), title('t=0.09s'),axis([0 1 0 1]),xticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]), yticks([0 0.2 0.4 0.6 0.8 1])
ylabel('Concentration');
xlabel('Distance');

subplot(3,2,6); plot(X, C_a(:,100),'b',X,C_n(:,100),'r'), title('t=0.1s'),axis([0 1 0 1]),xticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]), yticks([0 0.2 0.4 0.6 0.8 1])
ylabel('Concentration');
xlabel('Distance');

sgtitle('Concentration for analytical and numerical solution with N=40')
