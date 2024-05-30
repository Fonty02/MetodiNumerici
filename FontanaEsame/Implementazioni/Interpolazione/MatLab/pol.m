function polynomial = pol(coef,x)
    % Controllo che la lunghezza di coef e x sia coerente
    assert(length(coef) == length(x), 'Il numero di coefficienti deve essere uguale al numero di elementi in x.');

    % Numero di coefficienti
    n = length(coef);

    % Inizializza la stringa della funzione polinomiale
    poly_string = sprintf('%f', coef(1));

    % Aggiungi i termini del polinomio
    for i = 2:n
        poly_string = strcat(poly_string, sprintf(' + %f', coef(i)));
        for j = 1:i-1
            poly_string = strcat(poly_string, sprintf(' .* (t - %f)', x(j)));
        end
    end

    % Valuta la stringa come una funzione anonima
    polynomial = str2func(['@(t)' poly_string]);
end
