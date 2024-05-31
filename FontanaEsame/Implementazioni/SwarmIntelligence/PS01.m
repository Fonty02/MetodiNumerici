function [x, v, pb, n_it, g_b, radius, lista] = min_swarm(n, x, pb, v, func, max_it, ub, lb, tol)
    % Inizializzazione della lista
    lista = {};
    n_it = 0;
    % Calcolo del miglior punto globale iniziale
    [~, best_index] = min(feval(func,pb));
    g_b = pb(best_index); 

    lista{1} = x;

    % Prima iterazione
    c1 = 0.5;
    c2 = 0.3;
    v = v + c1 * rand(1, n) .* (pb - x) + c2 * rand(1, n) .* (g_b - x);
    x = x + v;

    % Controllo se x va oltre i limiti dello spazio di ricerca
    if min(x) < lb
        ind_m = find(x < lb, 1);
        x(ind_m) = lb;
    end
    if max(x) > ub
        ind_M = find(x > ub, 1);
        x(ind_M) = ub;
    end

    lista{2} = x;

    bool_vec = feval(func,x) < feval(func,pb);
    pb(bool_vec) = x(bool_vec);
    [~, best_index] = min(feval(func,pb));
    g_b = pb(best_index);

    radius = (max(abs(x - g_b))) / (abs(ub - lb));

    while (radius > tol) && (n_it < max_it)
        temp_x = x;
        temp_v = v;
        temp_pb = pb;
        v = temp_v + c1 * rand(1, n) .* (temp_pb - temp_x) + c2 * rand(1, n) .* (g_b - temp_x);
        x = temp_x + v;
        if min(x) < lb
            ind_m = find(x < lb, 1);
            x(ind_m) = lb;
        end
        if max(x) > ub
            ind_M = find(x > ub, 1);
            x(ind_M) = ub;
        end

        lista{end + 1} = x;
        n_it = n_it + 1;
        bool_vec = func(x) <feval(func,pb);
        pb(bool_vec) = x(bool_vec);
        [~, best_index] = mmin(feval(func,pb));
        g_b = pb(best_index);
        radius = max(abs(x - g_b)) / (abs(ub - lb));
    end

end
