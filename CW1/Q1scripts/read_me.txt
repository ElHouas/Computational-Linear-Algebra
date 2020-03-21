1. upper_triangularisation.m is the main file.

2. In order to run this on your machine, copy all these files in a folder and point your MATLAB directory to it. Then do the following:
	
(i) Define your matrix.
	
(ii) Call your upper_triangulisation subroutine by passing your matrix and its size. 
	
For example:
	
>> clear
	
>> A= [1,4,7; 2,5,8; 3,6,10]

;	
>> U = upper_triangulisation(A,3)
	
	U =
	
	     1     4     7
	     0    -3    -6
	     0     0     1
	
	
>> 

3. The main function (i.e. upper_triangulisation.m) calls five subroutines to perform the assigned duties. These subroutines are in-line with the theory explained in class, and the details of these subroutines are clearly written in the respective files.

4. New scripts are expected to meet the standards maintained in the scripts uploaded here (i.e. a neat presentation and proper commenting of what you are doing in a particular line of the script). 
Best wishes!
