function R = InitializeConnectivityMatrix(params)
% ģ����ǿ�Ӻ�������֮����໥����

% Construct the connectivity matrix
% �����ڽӾ���
R = zeros(params.beads_num); % (100,100)

% ���������Խ����ϵ�ֵΪ-1
R(diag(true(1,params.beads_num-1),1))  = -1;  % super diagonal
R(diag(true(1,params.beads_num-1),-1)) = -1;  % sub-diagonal

% Multiply by the spring constant values for heterogeneous polymer
% �������ʾۺ���ĵ��ɳ���ֵ
R = R.*params.spring_const;

A = zeros(params.beads_num); % (100,100)
% ������ǿ�Ӻ�������֮���������
A(params.enhancer_index,params.promoter_index) = -1;
A(params.promoter_index,params.enhancer_index) = -1;
A = A.*(params.attraction_coef);

% Sum rows to get the diagonal element value
% ������Ի�öԽ�Ԫ��ֵ
d = diag(true(1,params.beads_num));
R(d) = -sum(A+R,2);
R = R + A;
end
