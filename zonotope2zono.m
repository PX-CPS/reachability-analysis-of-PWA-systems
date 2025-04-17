function zonoObj = zonotope2zono(zonotopeObj)
    % Extract the full matrix [c, G] from zonotope object
    Z = zonotopeObj.Z;
    
    % Separate center vector and generator matrix
    c = Z(:,1);      % Center vector
    G = Z(:,2:end);  % Generator matrix
    
    % Create zono object
    zonoObj = zono(G, c);
end
