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

