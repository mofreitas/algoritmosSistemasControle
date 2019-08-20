function[] = observador(A, B, C, polos)
    
    V = zeros(length(polos), length(polos));
    for i=1:length(polos),
        V(i,:) = C*(A^(i-1));
    end
    
    printf("Matriz V/Wo ------------------------\n");
    disp(V);
    
    if(rank(V) == size(V, 'r')) then
        printf("\nSistema é observável\n");
    else
        halt("\nSistema NÃO é observável \n");
    end
    
    delta = poly(polos, 't');
    printf("\nDelta ------------------------\n");
    disp(delta);
    
    qo = zeros(length(polos), length(polos));
    for i=0:length(polos),
        qo = qo+(A^i)*coeff(delta, i);
    end
    
    printf("\nMatriz qo(A)/qo(G)---------------------\n");
    disp(qo);
    
    V_1 = inv(V);
    printf("\nMatriz V^(-1)/Wo^(-1)---------------------\n");
    disp(V_1);
    
    t = zeros(length(polos), 1);
    t($, 1) = 1;
    L = qo*V_1*t;
    printf("\nL---------------------\n");
    disp(L);
endfunction
