clear all
%% Question 4
% Demonstrate that the matrix is diagonal dominant 

a=0;   %Initialize the first value of alpha
f = 0; %Initialize condition that states when alpha makes the matrix non diagonal dominant

while (f~=1) %Condition to break the loop. If f=1, alpha will not satisfy
             %the condition for a non diagonal dominant matrix
             
B=[1 a 0; a 1 0; 0 a 1]; %The matrix defined in the coursework

[m,n]=size(B); %Get the size paramters of B

R=zeros(1,3); %Vector which will store the sum of the absolute values of 
              %each row
              
    for i=1:m     %We loop until the end of the row
        S=0;      %Initialize the sum of the non-diagonal elements 
        for j=1:n
            %Sum the absolute values of a row without the values
            %on the main diagonal.
            if i~=j %Check if the value is not placed in the main diagonal.
                S=S+abs(B(i,j)); %Sum the absolute value
            end
        end
        %We save each sum in the vector R.
        R(i)=S;
    end

    Delta = 0.0001; %Increment values for alpha

    for i=1:m
        %Compare the sum of the absolute values on the row with the absolute 
        %values on the main diagonal.
        if R(i)>B(i,i)
            f=1; %An alpha value which makes the matrix non diagonal dominant
                 %has been FOUND
        else
            a=a+Delta; %Delta is added to alpha until the condition of no diagonal 
                       %dominant is not satisfied
        end
    end
end