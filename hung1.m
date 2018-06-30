function hung1
    clc();
    DEBUG = true;
    MAXIMIZATION = true;
    % Данные
    
    A = dlmread('given.txt');
    fprintf('Original:\n');
    disp(A);
    
    if (MAXIMIZATION)
        M = abs(A - max(A));
        %M=A;
        fprintf('Maximization:\n');
        %disp(M);
    else
        M = A;
        fprintf('Minimization:\n');
        %disp(M);
    end
    %fprintf('\n');
    disp(M);
    fprintf('Subtract by Column\n');
    
    M = M - min(M);
    %fprintf('\n');
    disp(M);
    M = M - repmat(min(M,[],2),1, size(M,1));
    fprintf('Subtract by Row\n');
    disp(M);
    
    
 end