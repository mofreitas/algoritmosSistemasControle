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
