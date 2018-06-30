function hungarian
    
    maxmin = false;
    debug = false;
    
    clc();
    
    % Данные
    A = dlmread('given.txt');

    % Проверка задачи максимизации или минимизации
    if (maxmin)
        M = abs(A - max(A));
    else
        M = A;
    end

    % ШАГ 1 (получение матрицы с нулями)
    M = M - min(M);
    M = M - repmat(min(M,[],2),1, size(M,1));
    zvezda = zeros(1,size(M,1));
    number = zeros(1,size(M,1));

    % ШАГ 2 (находим все 0*)
    for i = 1:size(M,1)
       nol = setdiff(find(M(i,:) == 0), zvezda);
       if (nol)
           zvezda(i) = nol(1);
       end
    end

    % ШАГ 3 (выделение покрытых строк и столбцов)
    rows = zeros(1,size(M,1));
    cols = zeros(1,size(M,1));
    for i = 1:length(zvezda)
        if (zvezda(i) != 0)
            cols(zvezda(i)) = zvezda(i);
        end
    end

    loop = 0;
    while (find(cols == 0))
        if (debug)
            fprintf('\nStart!\n\n');
            show(M, zvezda, number, rows, cols, ++loop);
        end

        l_row = 0;
        l_col = 0;
        while (l_row == 0)

            % ШАГ 4 (находим 0')
            while (!isempty(find(M(!rows,!cols) == 0, 1)))
                for i = 1:size(M,1)
                    if (rows(i) ~= 0)
                        continue;
                    end
                    zero = setdiff(find(M(i,:) == 0), cols);
                    if (numel(zero) > 0)
                        if (zvezda(i) == 0)
                            l_row = i;
                            l_col = zero(1);
                            number(i) = l_col;
                            break;
                        end

                        number(i) = zero(1);
                        rows(i) = i;
                        cols(zvezda(i)) = 0;

                        if (debug)
                            show(M, zvezda, number, rows, cols, ++loop);
                        end
                    end
                end

                if (l_row > 0)
                    break;
                end
            end

            % ШАГ 5 (если 0' не найдены)
            if (l_row == 0)
                MM = M(!rows,!cols);
                mn = min(MM(MM>0));
                M = M + (sign(cols).*mn);
                M = M - diag(!sign(rows).*mn)*ones(size(M));
            end
        end
        
        if (debug)
            show(M, zvezda, number, rows, cols, ++loop);
        end

        % ШАГ 6 (построение L-цепочки)
        while true
            if (isempty(find(zvezda(:) == l_col,1)))
                zvezda(l_row) = l_col;
                break;
            end

            n_row = find(zvezda(:) == l_col, 1);
            zvezda(l_row) = l_col;
            l_row = n_row;
            l_col = number(l_row);
        end

        % ШАГ 3 (выделение покрытых строк и столбцов)
        rows = zeros(1,size(M,1));
        cols = zeros(1,size(M,1));
        for i = 1:length(zvezda)
            if (zvezda(i) != 0)
                cols(zvezda(i)) = zvezda(i);
            end
        end
    end

    if (debug && loop!=0)
        fprintf('Finished at iteration no.%d\n\nFinish!\n\n',loop);
    elseif(debug && loop==0)
        fprintf('Start!\n\n');
        show(M, zvezda, number, rows, cols, ++loop);
        fprintf('Finished at iteration no.%d\n\nFinish!\n\n',loop);
    end

    % ШАГ 7 (построение матрицы назначений)
    R = zeros(size(M));
    for i = 1:size(M,1)
        R(i,zvezda(i)) = 1;
    end
      
    result(maxmin, A, R);
end

% Вывод итерации
function show(M, zvezda, number, rows, cols, loop)
  fprintf('Interation no.%d\n', loop);
  for i = 1:size(M,1)
      for j = 1:size(M,2)
          if(zvezda(i) == j)
             fprintf('%5d ',M(i,j));
          elseif(number(i) == j)
             fprintf('%5d ',M(i,j));
          elseif(rows(i) != 0 || cols(j) != 0)
             fprintf('%5d ',M(i,j));
          else
             fprintf('%5d ',M(i,j));
          end
      end
      fprintf('\n');
  end
  fprintf('\n');
end

% Вывод результата
function result(maxmin, A, R)
  if(maxmin)
     fprintf('\nMaximum assignment problem:\n');      
  else
     fprintf('\nMinimum assignment problem:\n');
  end
  fprintf('\nCost matrix:\n');
  disp(A);
  fprintf('\nAssignment matrix:\n');
  disp(R);
  fprintf('\nF(opt) = %d\n\n', sum(sum(A.*R)));
end
