i=%i;

//Observador de estados
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

//Realimentador de estados
function[] = realimentador(A, B, C, polos)
    
    U = zeros(length(polos), length(polos));
    for i=1:length(polos),
        U(:,i) = (A^(i-1))*B;
    end
    
    printf("Matriz U/Wc ------------------------\n");
    disp(U);
    
    if(rank(U) == size(U, 'r')) then
        printf("\nSistema é controlável\n");
    else
        halt("\nSistema NÃO é controlável \n");
    end
    
    delta = poly(polos, 't');
    printf("\nDelta ------------------------\n");
    disp(delta);
    
    qc = zeros(length(polos), length(polos));
    for i=0:length(polos),
        qc = qc+(A^i)*coeff(delta, i);
    end
    
    printf("\nMatriz qc(A)/qc(G)---------------------\n");
    disp(qc);
    
    U_1 = inv(U);
    printf("\nMatriz U^(-1)/Wc^(-1)---------------------\n");
    disp(U_1);
    
    t = zeros(1, length(polos));
    t(1, $) = 1;
    K = -t*U_1*qc;
    printf("\nK---------------------\n");
    disp(K);
endfunction

//Obtem matrix de observabilidade
function[V] = obs(A, C)
    V = zeros(size(A, 'r'), size(A, 'c'));
    for i=1:size(A, 'r'),
        V(i, :) = C*A^(i-1);
    end
endfunction

//Obtem matrix de controlabilidade
function[U] = cont(A, B)
    U = zeros(size(A, 'r'), size(A, 'c'));
    for i=1:size(A, 'r'),
        U(:, i) = A^(i-1)*B;
    end
endfunction

//Seguidor de referência discreto
function[] = sr_d(G, H, C, polos)
    
    Ga = zeros(size(G, 'r')+1, size(G,'c')+size(H, 'c'));
    Ga(1:$-1, 1:size(G,'c')) = G;
    Ga(1:$-1, size(G, 'c')+1:$) = H;
    
    Ha = zeros(size(Ga, 'r'), 1);
    Ha($, 1) = 1;
    
    printf("Matriz Ga -------------------------- \n");
    disp(Ga);
    
    printf("\nMatriz Ha -------------------------- \n");
    disp(Ha);
    
    Wc = zeros(length(Ha), length(Ha));
    for i=1:length(Ha),
        Wc(:, i) = (Ga^(i-1))*Ha;
    end
    
    printf("\nMatrix Wc ------------------------ \n");
    disp(Wc);
    
    if(rank(Wc) == length(Ha)) then
        printf("\nRank cheio \n");
    else
        halt("\nNão tem Rank cheio \n");
    end
    
    printf("\nPolinomio delta -------------------- \n");
    delta = poly(polos, 'G');
    disp(delta);
    
    qc = zeros(size(Ga, 'r'), size(Ga, 'c'));
    coeficientes = coeff(delta);
    for i=1:length(coeficientes),
        qc = qc + (Ga^(i-1))*coeficientes(i);
    end
    printf("\nqc(Ga) ------------------------- \n");
    disp(qc);
    
    printf("\nWc inversa ---------------------- \n");
    disp(inv(Wc));
    
    m = zeros(1, size(qc, 'c'));
    m($) = 1;
    
    ka = -m*inv(Wc)*qc;
    printf("\nKa ------------------ \n");
    disp(ka);
    
    t1 = [G-eye(size(G, 'r'), size(G, 'r')) H; C*G C*H];
    t2 = zeros(1, length(ka));
    t2(1, $) = 1;
      
    printf("\nTa ------------------ \n");  
    disp(t1);
    
    k = (ka+t2)*inv(t1);
    printf("\nK ------------------ \n");
    disp(k);

endfunction

//Seguidor de referencia contínuo
function[] = sr_c(a, b, c, polos)
    
    Aa = zeros(size(c, 'r')+size(a, 'r'), size(c,'c')+1);
    Aa(1, 2:$) = c;
    Aa(2:$, 2:$) = a;
    
    Ba = [0; b]
    
    printf("Matriz Aa -------------------------- \n");
    disp(Aa);
    
    printf("\nMatriz Ba -------------------------- \n");
    disp(Ba);
    
    U = zeros(length(Ba), length(Ba));
    for i=1:length(Ba),
        U(:, i) = (Aa^(i-1))*Ba;
    end
    
    printf("\nMatrix U ------------------------ \n");
    disp(U);
    
    if(rank(U) == length(Ba)) then
        printf("\nRank cheio \n");
    else
        halt("\nNão tem Rank cheio \n");
    end
    
    printf("\nPolinomio delta -------------------- \n");
    delta = poly(polos, 'A');
    disp(delta);
    
    qc = zeros(size(Aa, 'r'), size(Aa, 'c'));
    for i=0:size(Aa, 'r'),
        qc = qc + (Aa^i)*coeff(delta, i);
    end
    printf("\nqc(Aa) ------------------------- \n");
    disp(qc);
    
    printf("\nU inversa ---------------------- \n");
    disp(inv(U));
    
    m = zeros(1, size(qc, 'c'));
    m($) = 1;
    
    ka = -m*inv(U)*qc;
    printf("\nKa ------------------ \n");
    disp(ka);

endfunction
