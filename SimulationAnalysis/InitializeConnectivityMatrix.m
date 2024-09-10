function R = InitializeConnectivityMatrix(params)
% 模拟增强子和启动子之间的相互作用

% Construct the connectivity matrix
% 构造邻接矩阵
R = zeros(params.beads_num); % (100,100)

% 设置两个对角线上的值为-1
R(diag(true(1,params.beads_num-1),1))  = -1;  % super diagonal
R(diag(true(1,params.beads_num-1),-1)) = -1;  % sub-diagonal

% Multiply by the spring constant values for heterogeneous polymer
% 乘以异质聚合物的弹簧常数值
R = R.*params.spring_const;

A = zeros(params.beads_num); % (100,100)
% 设置增强子和启动子之间的吸引力
A(params.enhancer_index,params.promoter_index) = -1;
A(params.promoter_index,params.enhancer_index) = -1;
A = A.*(params.attraction_coef);

% Sum rows to get the diagonal element value
% 求和行以获得对角元素值
d = diag(true(1,params.beads_num));
R(d) = -sum(A+R,2);
R = R + A;
end
