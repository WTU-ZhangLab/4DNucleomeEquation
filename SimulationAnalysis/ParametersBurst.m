function params = ParametersBurst(params)

params.attraction_coef = 1; % attraction interaction constant �����໥���ó���
params.enhancer_index = 25; % enhancer index number ��ǿ������
params.promoter_index = 75; % promoter index number ����������

params.simulated_on   = true;
params.EP_flag        = true;
params.compute_TR     = true;
params.dimension      = 3; % spatial dimension �ռ�ά��
params.beads_num      = 100; % beads number ����
params.dt             = 0.01;  % time step [sec] ʱ�䲽��
params.diffusion_const = 4e-3/params.friction_coef; % diffusion constant [mu m]^2 /[sec] ��ɢ����
params.b              = 0.1; % STD of distance between adjacent monomers ���ڵ���֮��ľ���ı�׼��
params.spring_const   = ones(params.beads_num); % spring constant ���ɳ���
params.factor         = sqrt(2*params.diffusion_const*params.dt);
params.distance_T     = 0.1; 
params.distance_05    = 0.2;
params.H              = 3;
params.encounter_dist = params.distance_T; % interaction distance ��������

params.k_on_max       = 1; % off to on max ��off̬��on̬������ٶ�
params.k_on           = 0.06;  % off to on min ��off̬��on̬����С�ٶ�
params.k_off          = 0.2; % on to off ��on̬��off̬���ٶ�
params.mu_max         = 3; % generate mRNA ����mRNA���������
params.mu             = 0.5; % generate mRNA ����mRNA������
params.delta          = 0.1; % mRNA degradation mRNA���������
params.friction_coef  = 50; % friction coefficient Ħ��ϵ��
params.simulation_num    = 32; % simulation number ģ�����
params.simulation_reaction_step = 1000000; % simulation reaction step ģ�ⷴӦ����
params.simulation_time   = 600000; % [s] ģ��ʱ�䣨�룩

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
% �����໥���ó���, ��ǿ������, ����������, ģ��, EP��־, ����TR, �ռ�ά��, ����, ʱ�䲽��, ��ɢ����, b, ���ɳ���, ����, ����T, ����05, H, ��������, k_on_max, k_on, k_off, mu_max, mu, ������, Ħ��ϵ��, ģ�����, ģ�ⷴӦ����, ģ��ʱ��, ���ռ��, ��������ļ���, �ļ���
% 1, 25, 75, true, true, true, 3, 100, 0.01, 0.004, 0.1, [100x100 double], 0.0632, 0.1, 0.2, 3, 0.1, 1, 0.06, 0.2, 3, 0.5, 0.1, 50, 50, 32, 1000000, 600000, 50, 'C:\Users\lyc\Documents\GitHub\BurstAnalysis\Results', []
end


