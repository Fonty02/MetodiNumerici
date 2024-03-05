function yy = lagrange(x, y, xx)
    % INPUT
    % x: vettore dei nodi
    % y: vettore delle ordinate
    % xx: vettore di ascisse in cui valutare pn(x)
    
    % OUTPUT
    % yy: vettore delle valutazioni di pn(x) nelle ascisse xx(i)

    n = length(x);  % numero di nodi
    yy = zeros(size(xx));
    
    for k = 1:n
        Lk = 1;
        for i = 1:n
            if i ~= k
                Lk = Lk .* ((xx - x(i)) / (x(k) - x(i)));
            end
        end
        yy = yy + (y(k) * Lk);
    end
end