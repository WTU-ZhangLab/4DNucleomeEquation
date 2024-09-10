function params = ParametersBurst(params)

params.attraction_coef = 1; % attraction interaction constant 吸引相互作用常数
params.enhancer_index = 25; % enhancer index number 增强子索引
params.promoter_index = 75; % promoter index number 启动子索引

params.simulated_on   = true;
params.EP_flag        = true;
params.compute_TR     = true;
params.dimension      = 3; % spatial dimension 空间维度
params.beads_num      = 100; % beads number 球数
params.dt             = 0.01;  % time step [sec] 时间步长
params.diffusion_const = 4e-3/params.friction_coef; % diffusion constant [mu m]^2 /[sec] 扩散常数
params.b              = 0.1; % STD of distance between adjacent monomers 相邻单体之间的距离的标准差
params.spring_const   = ones(params.beads_num); % spring constant 弹簧常数
params.factor         = sqrt(2*params.diffusion_const*params.dt);
params.distance_T     = 0.1; 
params.distance_05    = 0.2;
params.H              = 3;
params.encounter_dist = params.distance_T; % interaction distance 交互距离

params.k_on_max       = 1; % off to on max 从off态到on态的最大速度
params.k_on           = 0.06;  % off to on min 从off态到on态的最小速度
params.k_off          = 0.2; % on to off 从on态到off态的速度
params.mu_max         = 3; % generate mRNA 产生mRNA的最大速率
params.mu             = 0.5; % generate mRNA 产生mRNA的速率
params.delta          = 0.1; % mRNA degradation mRNA降解的速率
params.friction_coef  = 50; % friction coefficient 摩擦系数
params.simulation_num    = 32; % simulation number 模拟次数
params.simulation_reaction_step = 1000000; % simulation reaction step 模拟反应步数
params.simulation_time   = 600000; % [s] 模拟时间（秒）

params.snapshot_interval = 50; % [s]

params.result_base_folder = fullfile(pwd,'Results');
params.filename           = [];

%% Read the input options
intput_params = fieldnames(params);
for i = 1:length(intput_params)
    name  = intput_params{i};
    value = params.(name);
    if isfield(params, name)
        params.(name) = value;
    end
end

% attraction_coef, enhancer_index, promoter_index, simulated_on, EP_flag, compute_TR, dimension, beads_num, dt, diffusion_const, b, spring_const, factor, distance_T, distance_05, H, encounter_dist, k_on_max, k_on, k_off, mu_max, mu, delta, friction_coef, simulation_num, simulation_reaction_step, simulation_time, snapshot_interval, result_base_folder, filename
% 吸引相互作用常数, 增强子索引, 启动子索引, 模拟, EP标志, 计算TR, 空间维度, 球数, 时间步长, 扩散常数, b, 弹簧常数, 因子, 距离T, 距离05, H, 相遇距离, k_on_max, k_on, k_off, mu_max, mu, 代尔塔, 摩擦系数, 模拟次数, 模拟反应步数, 模拟时间, 快照间隔, 结果基础文件夹, 文件名
% 1, 25, 75, true, true, true, 3, 100, 0.01, 0.004, 0.1, [100x100 double], 0.0632, 0.1, 0.2, 3, 0.1, 1, 0.06, 0.2, 3, 0.5, 0.1, 50, 50, 32, 1000000, 600000, 50, 'C:\Users\lyc\Documents\GitHub\BurstAnalysis\Results', []
end


