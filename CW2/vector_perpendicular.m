function[v] =  vector_perpendicular(x,m,i)
    %Get the sign of the first element of the eigenvector found by RQI.
    %This is important in terms of stability. Because the angle between the
    %hyperplane and the vector will be small which could let to
    %cancellation errors.
    p=sign(x(1));
    %Get the euclidian longitude of the vector, by using the parameter 2  
    %in norm(x,2) and multiply it by the sign of the first element
    alpha=norm(x,2)*p;
    %Construct the perpendicular vector shrinking the size at each
    %iteration
    v=x+alpha*eye(m+1-i,1);
    %Normalize the perpendicular vector just computed
    v_norm=norm(v);
    %Divide the vector by its norm
    v=v/v_norm;
end