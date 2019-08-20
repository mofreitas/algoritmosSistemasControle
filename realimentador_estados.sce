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
    
//Se
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
    
    ka = m*inv(Wc)*qc;
    printf("\nKa ------------------ \n");
    disp(ka);
    
    t1 = [G-eye(size(G, 'r'), size(G, 'r')) H; C*G C*H];
    t2 = zeros(1, length(ka));
    t2(1, $) = 1;
    
    k = (ka+t2)*inv(t1);
    printf("\nK ------------------ \n");
    disp(k);

endfunction
