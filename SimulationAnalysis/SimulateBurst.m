function result = SimulateBurst(params,s_idx)
% This function simulates the dynamic evolution of E-P distance and transcriptional burst.
% 这个函数模拟E-P距离和转录爆发的动态演变。

%% Initialization
% pre-allocating transcriptional variable and Initialization
% 初始化转录变量

% recording the state after each reaction
% 记录每个反应后的状态
transcription_state = zeros(params.simulation_reaction_step,3); % (1000000,3)
transcription_state(1,1) = 1;

% recording the time at each reaction
% 记录每个反应的时间
transcription_time  = zeros(params.simulation_reaction_step,1); % (1000000,1)
transcription_time(1)  = 0;

% 记录每个反应之间的时间间隔
tau_array   = zeros(params.simulation_reaction_step,1); % (1000000,1)

% 记录每个反应的类型
reaction_array = zeros(params.simulation_reaction_step,1); % (1000000,1)

% 反应速率矩阵
r_mu  = [-1  1  0
    1   -1  0
    0   0  1
    0   0  -1] ;

% pre-allocating enhancer-promoter variable and Initialization
% 初始化增强子-启动子变量
if params.EP_flag

    % 初始化E-P距离
    current_position = zeros(params.beads_num, 3); % 每个单体的当前位置 (100,3)
    rand_array       = params.b.*randn(params.beads_num,params.dimension); % 生成随机位置 0.1*randn(100,3)
    % 从第二个单体开始，到100个单体，每个单体的位置等于前一个单体的位置加上随机数组
    for ps_idx = 2:params.beads_num
        current_position(ps_idx,:) = current_position(ps_idx-1,:) + rand_array(ps_idx,:);
    end

    % initialize rouse chain, harmonic potential between consecutive monomers
    % 初始化Rouse链，相邻单体之间的谐波势
    interaction_matrix = InitializeConnectivityMatrix(params);
    % an equilibratio of 1e4 steps was run before taking samples
    % 1e4步的平衡运行之后才取样
    randn_array = randn(params.beads_num,params.dimension,params.simulation_reaction_step).*params.factor;
    for index = 1:params.simulation_reaction_step
        current_position = current_position -(1./params.friction_coef).*...
            interaction_matrix*(current_position*params.dt) + randn_array(:,:,index);
    end
    EP_distance = pdist([current_position(params.enhancer_index,:);...
        current_position(params.promoter_index,:)]);
    i = 1;
    randn_array = randn(params.beads_num,params.dimension,params.simulation_reaction_step).*params.factor;
end
%% Initialization end

%% Main simulation loop
tic;
current_idx = 2; % 从第二个时间点开始
record_steps = 0; % 记录步数
while current_idx <= params.simulation_reaction_step && transcription_time(current_idx - 1) < params.simulation_time
    koff = params.k_off; delta = params.delta;
    if params.EP_flag
        u1 = rand(1); H = 1; steps = 0;
        while H > u1
            kon =  (EP_distance <= params.distance_T)*params.k_on_max + ...
                (EP_distance > params.distance_T)*(params.k_on + (params.k_on_max...
                - params.k_on)/(1+((EP_distance-params.distance_T)/(params.distance_05...
                -params.distance_T))^params.H));
            mu =  (EP_distance <= params.distance_T)*params.mu_max + ...
                (EP_distance > params.distance_T)*(params.mu + (params.mu_max...
                - params.mu)/(1+((EP_distance-params.distance_T)/(params.distance_05...
                -params.distance_T))^params.H));
            f_mu = [kon*transcription_state(current_idx - 1,1), koff*transcription_state(current_idx - 1,2),...
                mu*transcription_state(current_idx - 1,2), delta*transcription_state(current_idx - 1,3)*...
                (transcription_state(current_idx - 1,3)>0)];
            a_tot = sum(f_mu);
            H = H - a_tot*H*params.dt;
            current_position = current_position - (1./params.friction_coef).*...
                interaction_matrix*(current_position*params.dt) + randn_array(:,:,i);
            steps = steps + 1; i = i + 1; record_steps = record_steps + 1;
            if i == 100000; i = 1;randn_array = randn(params.beads_num,params.dimension,100000).*params.factor; end
            EP_distance = pdist([current_position(params.enhancer_index,:);current_position(params.promoter_index,:)]);
        end
        tau = steps*params.dt;
    else
        kon = params.k_on; mu = params.mu;
        f_mu = [kon*transcription_state(current_idx - 1,1), koff*transcription_state(current_idx - 1,2),...
            mu*transcription_state(current_idx - 1,2), delta*transcription_state(current_idx - 1,3)*...
            (transcription_state(current_idx - 1,3)>0)];
        a_tot = sum(f_mu);
        u1 = rand(1);
        tau = roundn((1/a_tot)*log(1/u1),-2);
    end
    u2 = rand(1)*a_tot;
    next_react = find(cumsum(f_mu) >= u2,1,'first');
    transcription_time(current_idx) = transcription_time(current_idx - 1) + tau;
    tau_array(current_idx) = tau;
    gain = r_mu(next_react,:);
    transcription_state(current_idx,:) = transcription_state(current_idx - 1,:) + gain;
    reaction_array(current_idx) = next_react;
    current_idx = current_idx + 1;
end % end loops

timerVal = toc;
disp(['Time taken for simulation',num2str(s_idx),': ',num2str(timerVal)]);

%% Taking sample
temp_snap_size = floor(transcription_time(current_idx - 1)/params.snapshot_interval);
mRNA_level = zeros(1, temp_snap_size);
for k = 1:temp_snap_size
    index = find(transcription_time > k*params.snapshot_interval, 1, 'first');
    snap_index = index - 1;
    mRNA_level(k) = transcription_state(snap_index,3);
end

%% Save simulation data
record_size = params.simulation_reaction_step;
if current_idx < params.simulation_reaction_step; record_size = current_idx - 1; end
result.params              = params;
result.record_size         = record_size;
result.transcription_state = transcription_state(1:record_size,:);
result.transcription_time  = transcription_time(1:record_size);
result.tau_array           = tau_array(1:record_size);
result.snap_size           = temp_snap_size;
result.mRNA_level          = mRNA_level;

filename = sprintf('//%d.mat', s_idx);
save([params.result_base_folder,params.filename,filename],'result');
end





