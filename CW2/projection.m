function u = projection(a,e)

%This function defines the projection of vector "a" over "e", 
%proj = (<a,e>/<e,e>) *e. It computes the inner product of <a,e> divided
%by the inner produt of <e,e> multiplied by the vector e.

u=((a'*e)/(e'*e))*e;
end