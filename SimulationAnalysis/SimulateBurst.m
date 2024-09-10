function result = SimulateBurst(params,s_idx)
% This function simulates the dynamic evolution of E-P distance and transcriptional burst.
% �������ģ��E-P�����ת¼�����Ķ�̬�ݱ䡣

%% Initialization
% pre-allocating transcriptional variable and Initialization
% ��ʼ��ת¼����

% recording the state after each reaction
% ��¼ÿ����Ӧ���״̬
transcription_state = zeros(params.simulation_reaction_step,3); % (1000000,3)
transcription_state(1,1) = 1;

% recording the time at each reaction
% ��¼ÿ����Ӧ��ʱ��
transcription_time  = zeros(params.simulation_reaction_step,1); % (1000000,1)
transcription_time(1)  = 0;

% ��¼ÿ����Ӧ֮���ʱ����
tau_array   = zeros(params.simulation_reaction_step,1); % (1000000,1)

% ��¼ÿ����Ӧ������
reaction_array = zeros(params.simulation_reaction_step,1); % (1000000,1)

% ��Ӧ���ʾ���
r_mu  = [-1  1  0
    1   -1  0
    0   0  1
    0   0  -1] ;

% pre-allocating enhancer-promoter variable and Initialization
% ��ʼ����ǿ��-�����ӱ���
if params.EP_flag

    % ��ʼ��E-P����
    current_position = zeros(params.beads_num, 3); % ÿ������ĵ�ǰλ�� (100,3)
    rand_array       = params.b.*randn(params.beads_num,params.dimension); % �������λ�� 0.1*randn(100,3)
    % �ӵڶ������忪ʼ����100�����壬ÿ�������λ�õ���ǰһ�������λ�ü����������
    for ps_idx = 2:params.beads_num
        current_position(ps_idx,:) = current_position(ps_idx-1,:) + rand_array(ps_idx,:);
    end

    % initialize rouse chain, harmonic potential between consecutive monomers
    % ��ʼ��Rouse�������ڵ���֮���г����
    interaction_matrix = InitializeConnectivityMatrix(params);
    % an equilibratio of 1e4 steps was run before taking samples
    % 1e4����ƽ������֮���ȡ��
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
current_idx = 2; % �ӵڶ���ʱ��㿪ʼ
record_steps = 0; % ��¼����
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





