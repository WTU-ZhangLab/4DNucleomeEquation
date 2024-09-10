clear;clc;
close all;
if isempty(gcp('nocreate'))
    numCores = feature('numcores');
    parpool(numCores);
end

% 定义增强子和启动子的索引
enhancer_idx = [51,55,60,65,70,75,80,85,90,95];
promoter_idx = [49,45,40,35,30,25,20,15,10,5];

input_options.attraction_coef   = 0.1; % K, E-P communication coefficient K，E-P通信系数
input_options.result_base_folder = fullfile(pwd, 'ResultDG');

% 循环遍历每对增强子和启动子
for idx = 1:length(enhancer_idx)
    input_options.enhancer_index    = enhancer_idx(idx); % enhancer index number 增强子索引
    input_options.promoter_index    = promoter_idx(idx); % promoter index number 启动子索引
    
    % 文件名，包括增强子索引、启动子索引、k_on_max、k_on、k_off、mu_max、mu、delta、friction_coef
    input_options.filename = sprintf('//E_%d_P_%d_%f_%f_%f_%f_%f_%f_%f',[input_options.enhancer_index,...
        input_options.promoter_index,input_options.k_on_max,input_options.k_on,input_options.k_off,...
        input_options.mu_max,input_options.mu,input_options.delta,input_options.friction_coef]);
    result_folder = fullfile(input_options.result_base_folder, input_options.filename);
    if ~exist(result_folder, 'dir')
        mkdir(result_folder);
    end

    params = ParametersBurst(input_options);

    if input_options.simulated_on
        tic;
        parfor s_idx = 1:input_options.simulation_num
            result = SimulateBurst(params,s_idx);
        end
        timerVal = toc;
        disp(['Total simulation time:', num2str(timerVal)]);
    end

    % Analysing Burst 分析爆发
    % Only parameters are needed, and data are loaded from the .mat files 只需要参数，数据从.mat文件中加载
    tic;
    results = AnalyseBurst(params);
    filename = sprintf('//E_%d_P_%d_%f_%f_%f_%f_%f_%f_%f.mat',[input_options.enhancer_index,...
        input_options.promoter_index,input_options.k_on_max,input_options.k_on,input_options.k_off,...
        input_options.mu_max,input_options.mu,input_options.delta,input_options.friction_coef]);
    save([input_options.result_base_folder,filename],'results');
    timerVal = toc;
    disp(['Analysing time:', num2str(timerVal)]);

    fig = figure;
    hold on
    if results.params.simulated_on == true
        h = histogram(results.mRNA_total,'BinEdges',0:max(results.mRNA_total),"Normalization",'probability','DisplayName','simulation');
    end
    plot(0.5+(0:1:length(results.mRNA_Prob)-1),results.mRNA_Prob,'LineWidth',1,'DisplayName','mixed')
%      plot(0.5+(0:1:length(results.mRNA_Prob)-1),results.PDF.mRNA_Prob_v1,'LineWidth',1,'DisplayName','fast')
%      plot(0.5+(0:1:length(results.mRNA_Prob)-1),results.PDF.mRNA_Prob_v2,'LineWidth',1,'DisplayName','slow')
    set(fig,'position',[300 400 280 190]);
    set(gca,'TickLength',[0.02,0.025]);
    legend1 = legend(gca,'show');
    box on
    xlim([0 200])
    title(['E-P ',num2str(abs(enhancer_idx(idx)-promoter_idx(idx))),...
        '  weight ',num2str(results.weight(1)) ]);
    % Save the figure
    saveas(fig, fullfile(input_options.result_base_folder, input_options.filename, 'figure.png'));
    data.mRNAdistri{1,idx} = results.mRNA_Prob;
    data.EPgenomedistance(idx) = abs(enhancer_idx(idx)-promoter_idx(idx));
    data.mRNAdistribin{1,idx} = 0:1:length(results.mRNA_Prob)-1;
end
