function zonotopeObj = zono2zonotope(zonoObj)
    % Extract center and generator from zono object
    c = zonoObj.c;
    G = zonoObj.G;
    
    % Combine center and generator matrix into one matrix as expected by zonotope class
    Z = [c, G];  % First column as center, remaining as generators
    
    % Create zonotope object
    zonotopeObj = zonotope(c, G);
end
