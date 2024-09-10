clear;clc;
close all;
if isempty(gcp('nocreate'))
    numCores = feature('numcores');
    parpool(numCores);
end

K_coefficient =  0.05:0.1:0.95; % K，E-P通信系数
input_options.enhancer_index    = 25; % enhancer index number 增强子索引
input_options.promoter_index    = 75; % promoter index number 启动子索引
input_options.result_base_folder = fullfile(pwd, 'ResultKEP');

% 循环遍历每个K系数
for idx =  1:length(K_coefficient)
    input_options.attraction_coef   = K_coefficient(idx); % K, E-P communication coefficient K，E-P通信系数

    input_options.filename = sprintf('//%f_%f_%f_%f_%f_%f_%f_%f',[input_options.attraction_coef,...
        input_options.k_on_max,input_options.k_on,input_options.k_off,input_options.mu_max...
        input_options.mu,input_options.delta,input_options.friction_coef]);
    if exist([input_options.result_base_folder,input_options.filename],'file') == 0
        mkdir(input_options.result_base_folder,input_options.filename);
    end

    params = ParametersBurst(input_options);
    
    if input_options.simulated_on
        tic;
        parfor s_idx = 1:input_options.simulation_num
            result = SimulateBurst(params,s_idx);
        end
        timerVal = toc;
        disp(['Total simulation time:',num2str(timerVal)]);
    end

    % Analyzing Burst 分析爆发
    % Only parameters are needed, and data are loaded from the .mat files 只需要参数，数据从.mat文件中加载
    tic;
    results = AnalyseBurst(params);
    filename = sprintf('//%f_%f_%f_%f_%f_%f_%f_%f.mat',[input_options.attraction_coef,...
        input_options.k_on_max,input_options.k_on,input_options.k_off,input_options.mu_max...
        input_options.mu,input_options.delta,input_options.friction_coef]);
    save([input_options.result_base_folder,filename],'results');
    timerVal = toc;
    disp(['Analysing time:',num2str(timerVal)]);
    
    fig = figure; 
    hold on
    if results.params.simulated_on == true
        h = histogram(results.mRNA_total,'BinEdges',0:max(results.mRNA_total),"Normalization",'probability','DisplayName','simulation');
        scatter(0.5+(0:3:length(h.BinEdges)-2),h.Values(1:3:end),'MarkerEdgeColor',[0.0392156876623631 0.658823549747467 0.250980406999588],...
    'LineWidth',1,...
    'MarkerFaceColor',[0.894117653369904 0.941176474094391 0.901960790157318]);
    end
    plot(0.5+(0:1:length(results.mRNA_Prob)-1),results.mRNA_Prob,'LineWidth',1,'DisplayName','mixed')
%     plot(0.5+(0:1:length(results.mRNA_Prob)-1),results.PDF.mRNA_Prob_v1,'LineWidth',1,'DisplayName','fast')
%     plot(0.5+(0:1:length(results.mRNA_Prob)-1),results.PDF.mRNA_Prob_v2,'LineWidth',1,'DisplayName','slow')
    set(fig,'position',[300 400 280 190]);
    set(gca,'TickLength',[0.02,0.025]);
    legend1 = legend(gca,'show');
    box on
%     title(['kep ',num2str(K_coefficient(idx)),'  weight ',num2str(results.weight) ]);
    % Save the figure
    saveas(fig, fullfile(input_options.result_base_folder, input_options.filename, 'figure.png'));
end
