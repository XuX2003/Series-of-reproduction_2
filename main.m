% A = [0 0 0 0 0 322.38;
%     0.966 0 0 0 0 0;
%     0.013 0.01 0.125 0 0 3.448;
%     0.007 0 0.125 0.238 0 30.17;
%     0.008 0 0.038 0.245 0.167 0.862;
%     0 0 0 0.023 0.75 0];
A = [0.1 7.7 1.9 1.7;
     9.1 26.2 0 18.6;
     0.1 1.3 8.3 0.4;
     2.1 18.5 0 4]; % 程老师测试
% A = [0 0 0 0 0 322.38;
%     0.966 0 0 0 0 0.3279;
%     0.013 0.01 0.125 0 0 3.448;
%     0.007 0 0.125 0.578 0 13.17;
%     0.008 0 0.038 0.245 0.167 0.862;
%     0 0 5.689 0.023 0.75 0];
% A = [0	0	0	0	2.30000000000000	7.20000000000000	12.8000000000000	9
% 0.200000000000000	0	0	0	0	0	0	0
% 0	0.600000000000000	0	0	0	0	0	0
% 0	0	0.800000000000000	0	0	4.70000000000000	0	0
% 0	0	0.100000000000000	0.800000000000000	0	0	0	0
% 0	0	0	0	0.860000000000000	0	0	0
% 0	0	0	0	0	0.870000000000000	0	0
% 0	0	0	0	0	0	0.750000000000000	0];

%% 判断矩阵是否为方阵
[nr, nc] = size(A);
if nr ~= nc 
    error('矩阵 A 并不是一个方阵！请重新检查矩阵！')
end
Num_Nodes = nr;  % 节点数量

%% 判断矩阵是否不可约
AI = A;
for j = 1 : Num_Nodes
    AI(j, j) = AI(j, j) + 1;
end
if (nnz(AI^(Num_Nodes - 1)) < Num_Nodes^2)
    error('矩阵 A 是一个可约矩阵！请重新检查矩阵！');
end

%% 判断矩阵是否非负
if ~isempty(find(A < 0, 1))
    error('矩阵 A 并不是一个非负矩阵！请重新检查矩阵！');
end

%% 计算矩阵的特征值（Lambda） 、特征向量、敏感度矩阵和弹性矩阵
[V, D] = eig(full(A));   
[~, Lambda_index] = max(max(real(D))); 
Right_EV = V(:, Lambda_index); % 右特征矩阵
Right_EV = Right_EV / sum(Right_EV); % 归一化
 
[V, D] = eig(full(A'));  
[Lambda, Lambda_index] = max(max(real(D)));
Left_EV = V(:, Lambda_index); % 左特征矩阵
Left_EV = Left_EV / Left_EV(1, 1);

% 计算敏感度矩阵
Sensitivity = (Left_EV * Right_EV') / (sum(Left_EV .* Right_EV));

% 计算弹性矩阵
Elasticity_Mat = 100 * (Sensitivity .* A) / Lambda;

%% 从弹性矩阵（Elasticity_Mat）中提取非零元素的位置对应的弹性值.
nzA = find(A ~= 0); % 获取A中非零元素的索引
Num_Edges = length(nzA); % 获取有效边的数量
Elasticity_Vect = Elasticity_Mat(nzA); % st. gamma * x = e中的e
[Node_List, Loop_Mat] = FindLoops(A); 
Loop_Cell = struct2cell(Node_List);
maxsize = 0;
for i = 1: length(Loop_Cell)
    if (length(Loop_Cell{i}) > maxsize)
        maxsize = length(Loop_Cell{i});
    end
end
maxsize = maxsize + (maxsize - 1);
Loop_Lengths = zeros(length(Loop_Cell), 1);
temp = Loop_Cell{1};
l = length(temp);
Loop_Lengths(1, 1) = l - 1;
Loop_List = 'Loop 1:';
    for j = 1:l
        Loop_List = [Loop_List ' ' num2str(temp(j))];
    end    
for i = 2: length(Loop_Cell)
    temp = Loop_Cell{i};
    l = length(temp);
    Loop_Lengths(i, 1) = l - 1;
    Loop_Line = ['Loop ' num2str(i) ':'];
    for j = 1: l
        Loop_Line = [Loop_Line ' ' num2str(temp(j))];
    end
    Loop_List = char(Loop_List, Loop_Line);
end    
Num_Loops = length(Loop_Lengths);
Cycle_Space = Num_Edges - Num_Nodes + 1;  % 循环空间维度
Solution_Set = Num_Loops - Cycle_Space;  % 解空间维度
loop_number=(1: Num_Loops)'; 

%% 如果解空间维度为0，说明存在唯一的循环分解，直接计算特征弹性（x）
if (Solution_Set == 0)
    x = Loop_Mat \ Elasticity_Vect; 
    Unique_Decomp = [loop_number, x, x .* Loop_Lengths];
    for i = 1: length(Unique_Decomp)
        disp(Unique_Decomp(i));
    end
    return
end

%% 计算和存储循环-边关联矩阵 Loop_Mat 的零空间（null space）的基向量，并找出每个基向量中正负分量的位置。
B = null(Loop_Mat, 'r'); % 计算矩阵 Loop_Mat 的零空间，并返回一个正交基
for j = 1: Solution_Set
    B_pos = find(B(:, j) == 1)';
    B_neg = find(B(:, j) == -1)';
    TradeOffs(j) = struct('basis_vector', j, 'B_pos', B_pos, 'B_neg', B_neg);
end

%% 确定特征弹性的上下界(Table 6)和循环的弹性(Table 5)
lb = zeros(Num_Loops, 1); 
char_elast = zeros(Num_Loops, 2);  
loop_elast = zeros(Num_Loops, 2);  

for i = 1: Num_Loops
    f = zeros(Num_Loops, 1); % 0 0 Length(i) 0 0 0 0 0
    f(i) = Loop_Lengths(i);
    x = linprog(f, [], [], Loop_Mat, Elasticity_Vect, lb, []);
    char_elast(i, 1) = x(i);
    x = linprog(-f, [], [], Loop_Mat, Elasticity_Vect, lb, []);
    char_elast(i, 2) = x(i); 
end

loop_elast(:, 1) = char_elast(:, 1) .* Loop_Lengths; % 将特征弹性乘以循环长度得到最小循环弹性.
loop_elast(:, 2) = char_elast(:, 2) .* Loop_Lengths; % 将特征弹性乘以循环长度得到最大循环弹性
 
options = optimoptions('linprog','Algorithm','dual-simplex-highs');
f = zeros(length(Loop_Lengths), 1);
% f = Loop_Lengths';
f(1: 4) = Loop_Lengths(1: 4);
f = f'
[x, fval1, exitflag1] = linprog(f, [], [], Loop_Mat, Elasticity_Vect, lb, [], options);
optimal_f(:, 1) = x .* Loop_Lengths;
[y, fval2, exitflag2] = linprog(-f, [], [], Loop_Mat, Elasticity_Vect, lb, [], options);
optimal_f(:, 2) = y .* Loop_Lengths;

if ((exitflag1 ~= 1) || (exitflag2 ~= 1))
    disp('未收敛到最优解');
    Exit_Flags = [exitflag1, exitflag2];
    disp(Exit_Flags(1));
    disp(Exit_Flags(2));
end

Decomp = [loop_number, optimal_f];