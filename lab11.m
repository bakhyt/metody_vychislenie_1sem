function lab11()
        clc();
        DEBUG = true;
        MAXIMIZATION = 1;
                
        C_orig = getMatrix();
        dispC(C_orig, 'Original C');

        if MAXIMIZATION
            C = maximize(C_orig);
            dispC(C, 'Maximization task C')
             
        else
            C = C_orig;
        end
        C = subColumnMin(C);
        C = subRowMin(C);
        dispC(C, 'After subtraction', DEBUG);
         
        [C, k] = markStarZeros(C);
 
          while k~= C.n %%крутим программу пока "к" не стал равен размерности матрицы
              C = markColumnsWithStars(C);%%выделяем с "+" все столбцы с 0-звездочкой
              dispC(C, 'Marked *', DEBUG);
              
              while true %%крутим программу пока среди невыделенных элекментов не осталис просто нолики
                  [flag, posShtrih] = unmarkedHasZero(C);%% есть ли среди невыделенных элементов просто нолики
                  if flag %% есть ли среди невыделенных элементов просто нолики
                      C.shtrih(posShtrih.r, posShtrih.c) = true; %% если есть то отмечаем такой ноль с штрихом 
                      dispC(C, 'Marked ''', DEBUG);
                      [flag, posStar] = shtrihedRowHasStar(C, posShtrih.r); %% А есть в одной строке с ноль-штрих ноль-звезда
                      if flag %% А есть ли в одной строке с ноль-штрих ноль-звезда
                          C.markedColumns(posStar.c) = false; % если есть таково тогда снимаем "+" выделение с тех столбцов 
                          C.markedRows(posShtrih.r) = true; % если есть таково тогда отмечаем "+" эти строки 
                          dispC(C, 'Unmark the column and Mark the row', DEBUG);
                      else %% Если в одной строке с ноль-штрих ноль-звезда тогда строим нашу цепочку L
                          L = buildLChain(C, posShtrih); % строим цепочку L
                          dispL(L, DEBUG);
                          C = swapZeros(C, L); %%в пределах этой цепочки меняем все ноль-штрихи с ноль-звездами, а ноль-звезды но обычные нолики 
                          dispC(C, 'Swapped * and ''', DEBUG);
                          C = unmarkPlusAndShtrih(C);%% снимаем все выделении
                          dispC(C, 'Unmarked all', DEBUG);
                          k = k + 1; %% увеличиваем "к" на один ("к" это был количество ноли с звездами
                          break %% останавливаем второй цикл и идем на первый цикл где сравнивает количество "к" с размерности матрицы
                      end
                  else %% если среди невыделенных элементов нет нолики 
                      H = selectMinPositiveHFromUnmarked(C); %%тогда из невыделенных элементов находим минимальный элемент
                      C = subAndAddH(C, H);%% отнимает этот найденный минимальный элемент вычитаем от всех невыделенных строк и лобавляем ко всем выделенным столбцам.
                      dispC(C, 'Subtracted min', DEBUG);
                  end
              end
          end
           
          disp('Result:');
           
          Xopt = formXopt(C);%% отмечаем результать матрицы в виде нолики и однерки
          fopt = calcF(C_orig, Xopt);%%сложим знечение тех элементов где выше покажет однерки
           
          fprintf('X opt is:\n');
          disp(Xopt);
          fprintf('\nf_opt is:\n');
          disp(fopt);
end
 
%% сформировать структуру C из отдельных матриц
function C = formC(m,n,star,shtrih,markedRows,markedColumns)
         C = struct('m',m, 'n',n, 'star',star, 'shtrih',shtrih, 'markedRows',markedRows, 'markedColumns',markedColumns);
end

function P = formPoint(r,c)
         P = struct('r',r,'c',c);
end
 
%% отобразить матрицу стоимостей со всеми + и *
function dispC(C, msg, show)
  m = C.m;
  n = C.n;
  arg3 = (nargin == 3);%%если количество параметров функции меньше трех то arg3 = 0, но если равен трем то arg3 = 1   
  if (~arg3 || (arg3 && show == true))
        fprintf('%s:\n', msg);
        for r = 1:n
            for c = 1:n
                fprintf('  %4g', m(r,c));
                if C.star(r,c)
                    fprintf('*');
                elseif C.shtrih(r,c)
                    fprintf('''');
                else
                    fprintf(' ');
                end
            end
            if C.markedRows(r)
                fprintf('  +\n');
            else
                fprintf('   \n');
            end
        end
        for c = 1:n
            if C.markedColumns(c)
                p = '+';
            else
                p = ' ';
            end
            fprintf('  %4c ', p);
        end
        fprintf('\n\n');
    end
end
 
%% отобразить L-цепочку
function dispL(L, show)
    if show
        fprintf('L-chain:\n');
        disp(L);
    end
end
 
%% 1. получить матрицу стоимостей задачи о назначениях
function C = getMatrix()
        m = dlmread('given.txt');
         
        [rows, cols] = size(m);
        assert(rows == cols, 'Matrix should be square!');
         
        z = zeros(rows,cols);
        C = formC(m, rows, z, z, z(:,1), z(1,:));
end
 
%% 1.1 привести матрицу C к задаче о максимизации
function C = maximize(C)
          m = C.m;
          m = -m + max(max(m));
          C = formC(m, C.n, C.star, C.shtrih, C.markedRows, C.markedColumns);
end
 
%% 2. из каждого столбца матрицы вычесть минимальное значение
function res = subColumnMin(C)
      m = C.m;
      n = C.n;
       
      for i = 1:n
          col = m(:,i);
          m(:,i) = col - min(col);
      end
       
      res = C;
      res.m = m;
end
 
%% 3. из каждой строки матрицы вычесть минимальное
function res = subRowMin(C)
      m = C.m;
      n = C.n;
       
      for i = 1:n
          row = m(i,:);
          m(i,:) = row - min(row);
      end
       
      res = C;
      res.m = m;
end
 
%% 4. просмотреть полученную матрицу по столбцам в поисках нулей, в одной
% строке с которыми нет 0*; отметить такие нули
function [Cres, k] = markStarZeros(C)
    m = C.m;
    n = C.n;
    k = 0;
    star = C.star;
 
  function f = rowHasStar(r1)
    for c1 = 1:n
        if star(r1,c1)
            f = true;
            return;
        end
    end
    f = false;
  end
 
        for c = 1:n
            for r = 1:n
                if m(r,c) == 0 && ~rowHasStar(r)
                    star(r,c) = true;
                    k = k+1;
                    break;
                end
            end
        end
 
  Cres = formC(m, n, star, C.shtrih, C.markedRows, C.markedColumns);
end
 
%% 5. отмечаем столбцы, содержащие 0*
function C = markColumnsWithStars(C)
    m = C.m;
    n = C.n;
    markedColumns = C.markedColumns;
    for c = 1:n
        for r = 1:n
            if C.star(r,c)
                markedColumns(c) = true;
                break;
            end;
        end
    end
    C = formC(m, n, C.star, C.shtrih, C.markedRows, markedColumns);
end
 
%% 6. есть ли среди невыделенных столбцов просто 0
function [flag, posShtrih] = unmarkedHasZero(C)
        m = C.m;
        n = C.n;
        flag = false;
        posShtrih = formPoint(0,0);
         
        for c = 1:n
            if ~C.markedColumns(c)
                for r = 1:n
                    if ~C.markedRows(r)
                        if m(r,c) == 0 && ~C.star(r,c) && ~C.shtrih(r,c)
                            flag = true;
                            posShtrih = formPoint(r,c);
                            return
                        end
                    end
                end
            end
        end
end
 
%% 7. есть ли 0* в одном ряду с 0'
function [flag, posStar]  = shtrihedRowHasStar(C, r)
      n = C.n;
      flag = false;
      posStar = formPoint(0,0);
      for c = 1:n
          if C.star(r,c)
              flag = true;
              posStar = formPoint(r,c);
              return;
          end
      end
end
 
%% 8. строим L-цепочку от 0'
function L = buildLChain(C, posShtrih)
          n = C.n;
          star = C.star;
          shtrih = C.shtrih;
              function pos = findStarInCol(startR, c_)
                  pos = [0 0];
                  %вверх
                  for r_ = startR:-1:1
                      if star(r_,c_)
                          pos = [r_ c_];
                          return
                      end
                  end
                  %вниз
                  for r_ = startR:n
                      if star(r_,c_)
                          pos = [r_ c_];
                          return
                      end
                  end
              end
              function pos = findShtrihInRow(r_, startC)
                  pos = [0 0];
                  %влево
                  for c_ = startC:-1:1
                      if shtrih(r_,c_)
                          pos = [r_ c_];
                          return;
                      end
                  end
                  %вправо
                  for c_ = startC:n
                      if shtrih(r_,c_)
                          pos = [r_, c_];
                          return;
                      end
                  end
              end
           
          L = [];
          turn = 1; % turn%2=1 - искать *, turn%2=0 - искать '
          curPos = [posShtrih.r, posShtrih.c];
          while curPos
              L(turn,:) = curPos;
              if mod(turn,2)
                  curPos = findStarInCol(curPos(1), curPos(2));
              else
                  curPos = findShtrihInRow(curPos(1), curPos(2));
              end
              turn = turn+1;
          end
end
 
%% 9. заменяем 0' и 0* в L-цепочке
function Cres = swapZeros(C, L)
      [rows, ~] = size(L);
      star = C.star;
      shtrih = C.shtrih;
      for i = 1:rows
          pos = L(i,:);
          r = pos(1); c = pos(2);
          if star(r,c)
              star(r,c) = false;
              shtrih(r,c) = true;
          elseif C.shtrih(r,c)
              shtrih(r,c) = false;
              star(r,c) = true;
          end     
      end
      Cres = formC(C.m, C.n, star, shtrih, C.markedRows, C.markedColumns);
end
 
%% 10. снимаем выделения
function Cres = unmarkPlusAndShtrih(C)
      z = zeros(C.n, C.n);
      Cres = formC(C.m, C.n, C.star, z, z(:,1), z(1,:));
end
 
%% 11. выбираем из неотмеченных элементов наименьший положительный
function H = selectMinPositiveHFromUnmarked(C)
      m = C.m;
      n = C.n;
       
      H = max(max(m));
      for r = 1:n
          if ~C.markedRows(r)
              for c = 1:n
                  if ~C.markedColumns(c)
                      H = min(H, m(r,c));
                  end
              end
          end
      end
end
 
%% 12. вычитаем его из неотмеченных строк и добавляем к выделенным столбцам
function Cres = subAndAddH(C, H)
        m = C.m;
        n = C.n;
        for r = 1:n
            if ~C.markedRows(r)
                m(r,:) = m(r,:) - H;
            end
        end
        for c = 1:n
            if C.markedColumns(c)
                m(:,c) = m(:,c) + H;
            end
        end
        Cres = formC(m, n, C.star, C.shtrih, C.markedRows, C.markedColumns);
end
 
%% ##. сформировать матрицу нулей Xopt
function X = formXopt(C)
      n = C.n;
       
      X = zeros(n,n);
      for r = 1:n
          for c = 1:n
              if C.star(r,c)
                  X(r,c) = 1;
              end
          end
      end
end
 
%% ##. вычислить целевую функцию от оптимального решения
function f = calcF(C, X)
      n = C.n;
      f = 0;
       
      for r = 1:n
          for c = 1:n
              f = f + C.m(r,c) * X(r,c);
          end
      end
end
