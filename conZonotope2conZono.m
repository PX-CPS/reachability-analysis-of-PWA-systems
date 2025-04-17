function conZonoObj = conZonotope2conZono(conZonotopeObj)
    % Converts conZonotope object to conZono object
    %
    % Inputs:
    %    conZonotopeObj - conZonotope object to be converted
    %
    % Outputs:
    %    conZonoObj - conZono object
    
    % Extract the Z matrix, constraint matrix, and constraint vector from conZonotope
    Z = conZonotopeObj.Z;  % Z matrix containing center and generators
    A = conZonotopeObj.A;  % Constraint matrix
    b = conZonotopeObj.b;  % Constraint vector

    % Separate the Z matrix into center and generator matrix
    c = Z(:, 1);        % First column is the center
    G = Z(:, 2:end);    % Remaining columns are the generators

    % Construct the conZono object using G, c, A, and b
    conZonoObj = conZono(G, c, A, b);
end
