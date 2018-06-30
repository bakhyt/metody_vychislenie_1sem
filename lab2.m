function PotentialMethod()
    %WeightMatrix = importdata(FileAdress);
    WeightMatrix = [390 80 60 170 80; 110 5 4 3 4; 190 3 2 5 5; 90 1 6 3 2];
    SelectionMatrix = CreateEmptyMatrix(WeightMatrix);
    outFigure = figure;
    Num = 0;
    PrintMatrix(WeightMatrix,'Исходная матрица',outFigure,Num);
    Num = Num + 1;
    SelectionMatrix = SelectStartSolution(SelectionMatrix);
    PrintMatrix(SelectionMatrix,'Начальное БДР',outFigure,Num);
    Num = Num + 1;
    ImprovedMatrix = ImproveSolution(SelectionMatrix,WeightMatrix,outFigure,Num);
    Num = ImprovedMatrix.Count;
    if (ImprovedMatrix.Result == 0 )
        PrintMatrix(ImprovedMatrix.Matrix,'Транспортная матрица',outFigure,Num);
    end
    Num = Num + 1;
    while (ImprovedMatrix.Result == 0)
        ImprovedMatrix = ImproveSolution(ImprovedMatrix.Matrix,WeightMatrix,outFigure,Num);
        Num = ImprovedMatrix.Count;
        if (ImprovedMatrix.Result == 0 )
            PrintMatrix(ImprovedMatrix.Matrix,'Транспортная матрица',outFigure,Num);
            Num = Num + 1;
        end        
    end
    F = CountF(ImprovedMatrix.Matrix,WeightMatrix);
    OutputText(strcat('F = ',num2str(F)),'Решение',outFigure,Num,size(WeightMatrix,1));
end
 
function PrintMatrix(Matrix, logText, outFigure,Num)
    [N, M] = size(Matrix);
    hT = zeros(N, M);
    K = 6;
    Y = 1.7 - (fix(Num / K)*0.6 + 0.1);
    X = (Num - K*fix(Num/K))*0.55 + 0.1;
    Height = 0.4;
    Width = 0.5;
    figure(outFigure);
    hA11 = axes('Position', [X Y Width Height], 'Units', 'centimeters',...
                'XLim', [0 M+1], 'YLim', [0 N+1], 'XTick', 1:M, 'YTick', 1:N,...
                'YDir', 'reverse', 'XAxisLocation', 'top');
     axes(hA11)
    title(strcat(num2str(Num),'. ',logText))
    fprintf(strcat(num2str(Num),'. ',logText,'\n'))
    for i=1:N
        for j=1:M
            C = Matrix(i, j);
            if (C >= 0)
                hT(i, j) = text(j, i, num2str(C), 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k');        
            end
            fprintf('%d ',C);
        end
        fprintf('\n');
    end
    fprintf('\n');    
end
 
function PrintMatrixWithSelection(Matrix, SelectionM, logText, outFigure, Num)
    [N, M] = size(Matrix);
    hT = zeros(N, M);
    K = 6;
    Y = 1.7 - (fix(Num / K)*0.6 + 0.1);
    X = (Num - K*fix(Num/K))*0.55 + 0.1;
    Height = 0.4;
    Width = 0.5;
    figure(outFigure);
    hA11 = axes('Position', [X Y Width Height], 'Units', 'centimeters',...
                'XLim', [0 M+1], 'YLim', [0 N+1], 'XTick', 1:M, 'YTick', 1:N,...
                'YDir', 'reverse', 'XAxisLocation', 'top');
 
    axes(hA11)
    title(strcat(num2str(Num),'. ',logText))
    for i=1:N
        for j=1:M
            C = Matrix(i, j);
            if (i == 1 || j == 1)
                hT(i, j) = text(j, i, num2str(C), 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k');
            else
                if (SelectionM(i,j) < 0)
                    hT(i, j) = text(j, i, '-', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k') ;
                else
                    if (SelectionM(i,j) == 1)
                        if (C < 0)
                            hT(i, j) = text(j, i, '+w', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');                            
                        else
                            hT(i, j) = text(j, i, strcat(num2str(C),'+w'), 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
                        end
                    else
                        if (SelectionM(i,j) == 2)
                            if (C < 0)
                                hT(i, j) = text(j, i, '-w', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'b');
                            else
                                hT(i, j) = text(j, i, strcat(num2str(C),'-w'), 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'b');
                            end
                        end
                    end
                end
            end
        end
    end    
end
 
function OutputText(textMessage, logText, outFigure, Num, N)
    K = 6;
    Y = 1.7 - (fix(Num / K)*0.6 + 0.1);
    X = (Num - K*fix(Num/K))*0.55 + 0.1;
    Height = 0.4;
    Width = 0.5;
    figure(outFigure);    
    hA11 = axes('Position', [X Y Width Height], 'Units', 'centimeters',...
                'Visible', 'off', 'YDir', 'reverse', 'XAxisLocation', 'top');
    axes(hA11)
    title(strcat(num2str(Num),'. ',logText))
    
    hT1 = text(0.5, 0.5, textMessage, 'HorizontalAlignment', 'center', 'FontSize', 15, 'Color', 'k');
    set(hT1);
    fprintf(textMessage);
    fprintf('\n');
end

function Matrix = CreateEmptyMatrix(WeightMatrix)
    [N, M] = size(WeightMatrix);
    Matrix = zeros(N,M);
    for i=1:N
        Matrix(i,1) = WeightMatrix(i,1);
    end
    for j=1:M
        Matrix(1,j) = WeightMatrix(1,j);
    end
    for i=2:N
        for j=2:M
            Matrix(i,j) = -1;
        end
    end
end
 
function SelectionMatrix = SelectStartSolution(SelectionMatrix)
    [N, M] = size(SelectionMatrix);
    iIndex = 2;
    jIndex = 2;
    
    while (iIndex <= N && jIndex <= M)
        if (SelectionMatrix(iIndex,1) < SelectionMatrix(1,jIndex))
            SelectionMatrix(iIndex,jIndex) = SelectionMatrix(iIndex,1);
            SelectionMatrix(1,jIndex) = SelectionMatrix(1,jIndex) - SelectionMatrix(iIndex,1);
            SelectionMatrix(iIndex,1) = 0;
            iIndex = iIndex + 1;
        else
            SelectionMatrix(iIndex,jIndex) = SelectionMatrix(1,jIndex);
            SelectionMatrix(iIndex,1) = SelectionMatrix(iIndex,1) - SelectionMatrix(1,jIndex);
            SelectionMatrix(1,jIndex) = 0;
            jIndex = jIndex + 1;
        end        
    end
end
 
function F = CountF(SelectionMatrix, WeightMatrix)
    [N, M] = size(SelectionMatrix);
    F = 0;
    for i=2:N
        for j=2:M
            if (SelectionMatrix(i,j) > 0)
                F = F + SelectionMatrix(i,j)*WeightMatrix(i,j);
            end
        end
    end
end
 
function ImproveResult = ImproveSolution(SelectionMatrix,WeightMatrix, outFigure, Num)
    ImproveResult = struct('Result',-1,'Matrix',SelectionMatrix,'Count',Num);
 
    [N, M] = size(SelectionMatrix);
    UVMatrix = FullUVMatrix(SelectionMatrix,WeightMatrix);
    UVMatrix = SolveUsingGauss(UVMatrix);
    UArray = GetUArray(UVMatrix,N-1);
    VArray = GetVArray(UVMatrix,M-1);
    MinStruct = FindMinD(SelectionMatrix,WeightMatrix,UArray,VArray);
    fprintf('Min = %d\n',MinStruct.Min);
    if (MinStruct.Min < 0)
        CycleMatrix = FindCycle(SelectionMatrix,MinStruct.iIndex,MinStruct.jIndex);
        PrintMatrixWithSelection(SelectionMatrix,CycleMatrix,'Матрица с циклом',outFigure,Num);
        MinElement = FindMinToBeDecreasedElementInCycle(SelectionMatrix,CycleMatrix);
        SelectionMatrix = ModifyMatrix(SelectionMatrix,CycleMatrix,MinElement);
        ImproveResult.Result = 0;
        ImproveResult.Matrix = SelectionMatrix;
        ImproveResult.Count = Num + 1;
    end
end
 
function UVMatrix = FullUVMatrix(SelectionMatrix,WeightMatrix)
    [N, M] = size(SelectionMatrix);
    N = N - 1;
    M = M - 1;
    UVRows = N+M;
    UVColumns = N+M+1;
    UVMatrix = zeros(UVRows,UVColumns);
    UVMatrix(1,1) = 1;
    UVMatrix(1,N+M+1) = 0;
    index = 2;
    for i=2:(N+1)
        for j=2:(M+1)
            if (SelectionMatrix(i,j) >= 0)
                UVMatrix(index,i-1) = 1;
                UVMatrix(index,N+j-1) = 1;
                UVMatrix(index,UVColumns) = WeightMatrix(i,j);
                index = index + 1;
            end
        end
    end
end
 
function UVMatrix = SolveUsingGauss(SelectionMatrix)
    [N, M] = size(SelectionMatrix);    
    for K=1:N
        Max = abs(SelectionMatrix(K,K));
        iMax = K;
        for i=(K+1):N
            if abs(SelectionMatrix(i,K)) > Max
                Max = abs(SelectionMatrix(i,K));
                iMax = i;
            end
        end
        if (iMax ~= K)
            for j=K:M
                temp = SelectionMatrix(K,j);
                SelectionMatrix(K,j) = SelectionMatrix(iMax,j);
                SelectionMatrix(iMax,j) = temp;
            end;
        end
        Koff = SelectionMatrix(K,K);
        SelectionMatrix(K,K) = 1;
        for j=(K+1):M
            SelectionMatrix(K,j) = SelectionMatrix(K,j) / Koff;
        end
        for i=(K+1):N
            Koff = SelectionMatrix(i,K);
            SelectionMatrix(i,K) = 0;
            for j=(K+1):M
                SelectionMatrix(i,j) = SelectionMatrix(i,j) - (Koff/SelectionMatrix(K,K)) * SelectionMatrix(K,j);
            end
        end
    end
    SelectionMatrix(N,M) = SelectionMatrix(N,M) / SelectionMatrix(N,N);
    K = N - 1;
    while (K >= 1)
        L = K + 1;
        R = SelectionMatrix(K,L)*SelectionMatrix(L,M);
        L = L + 1;
        while (L <= N)
            R = R + SelectionMatrix(K,L)*SelectionMatrix(L,M);
            L = L + 1;
        end
        SelectionMatrix(K,M) = (SelectionMatrix(K,M) - R) / SelectionMatrix(K,K);
        K = K - 1;
    end;
    UVMatrix = SelectionMatrix;
end
 
function UArray = GetUArray(UVMatrix,UCount)
    [~, M] = size(UVMatrix);
    UArray = zeros(UCount);
    for i=1:UCount
        UArray(i) = UVMatrix(i,M);
    end    
end
 
function VArray = GetVArray(UVMatrix,VCount)
    [N, M] = size(UVMatrix);
    VArray = zeros(VCount);
    delta = N-VCount;
    for i=1:VCount
        VArray(i) = UVMatrix(i+delta,M);
    end    
end
 
function MinStruct = FindMinD(SelectionMatrix,WeightMatrix,UArray,VArray)
    [N, M] = size(SelectionMatrix);    
    MinStruct = struct('Min',1,'iIndex',-1,'jIndex',-1);
    for i=2:N
        for j=2:M
            if (SelectionMatrix(i,j) < 0)
               temp = WeightMatrix(i,j) - UArray(i-1) - VArray(j-1);
               if (temp < MinStruct.Min)
                   MinStruct.Min = temp;
                   MinStruct.iIndex = i;
                   MinStruct.jIndex = j;
               end
            end            
        end
    end
end
 
function CycleMatrix = FindCycle(SelectionMatrix,iIndex,jIndex)
    CycleMatrix = CreateEmptyMatrix(SelectionMatrix);
    CycleMatrix(iIndex,jIndex) = 1;
    Result = GetCycleIndex(SelectionMatrix,CycleMatrix,2,jIndex);
    CycleMatrix = Result.Matrix;
end
 
function Result = GetCycleIndex(SelectionMatrix,CycleMatrix,FlagRowOrColumn,Index)
    Result = struct('Flag',1,'Matrix',CycleMatrix);
    [N, M] = size(SelectionMatrix);
    if (FlagRowOrColumn == 1) % Row
        for j=1:M
            if (CycleMatrix(Index,j) == 1) % Ð•Ð°Ð¼Ñ‹ÐºÐ°Ð½Ð¸Ðµ Ñ†Ð¸ÐºÐ»Ð°
                Result.Flag = 0;
                Result.Matrix = CycleMatrix;
                return;
            end
        end        
        for j=1:M
            if (SelectionMatrix(Index,j) > -1 && CycleMatrix(Index,j) == -1)
                CycleMatrix(Index,j) = 1;
                Result = GetCycleIndex(SelectionMatrix,CycleMatrix,2,j);
                if (Result.Flag == 0) % ÑƒÑÐ¿ÐµÑ…
                    return;
                end
                CycleMatrix(Index,j) = -1;    
            end
        end
    else % Column
        for i=1:N
            if (SelectionMatrix(i,Index) > -1 && CycleMatrix(i,Index) == -1)
                CycleMatrix(i,Index) = 2;
                Result = GetCycleIndex(SelectionMatrix,CycleMatrix,1,i);
                if (Result.Flag == 0) % ÑƒÑÐ¿ÐµÑ…
                    return;
                end
                CycleMatrix(i,Index) = -1;    
            end
        end
    end
end
function MinElement = FindMinToBeDecreasedElementInCycle(SelectionMatrix,CycleMatrix)
    MinElement = SelectionMatrix(1,1);
    [N, M] = size(SelectionMatrix);
    
    for i=2:N
        for j=2:M
            if (CycleMatrix(i,j) == 2 && SelectionMatrix(i,j) < MinElement)
                MinElement = SelectionMatrix(i,j);
            end
        end
    end
end
function SelectionMatrix = ModifyMatrix(SelectionMatrix,CycleMatrix,MinElement)
    [N, M] = size(SelectionMatrix);
    DeleteFlag = 0;
    for i=2:N
        for j=2:M
            if (CycleMatrix(i,j) == 1)
                if (SelectionMatrix(i,j) == -1)
                    SelectionMatrix(i,j) = MinElement;
                else
                    SelectionMatrix(i,j) = SelectionMatrix(i,j) + MinElement;
                end
            else
                if (CycleMatrix(i,j) == 2)
                    SelectionMatrix(i,j) = SelectionMatrix(i,j) - MinElement;
                    if (SelectionMatrix(i,j) == 0 && DeleteFlag == 0)
                        DeleteFlag = 1;
                        SelectionMatrix(i,j) = -1;
                    end
                end
            end
        end
    end
end
