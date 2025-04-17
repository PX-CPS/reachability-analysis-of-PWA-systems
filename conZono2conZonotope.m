function czZ = conZono2conZonotope(cz)
    % Validate the input
    if ~isa(cz, 'conZono')
        error('Input must be a conZono object.');
    end

    % Extract properties
    G = cz.G;
    c = cz.c;
    A = cz.A;
    b = cz.b;

    % Create the conZonotope object
    czZ = conZonotope(c, G, A, b);
end
