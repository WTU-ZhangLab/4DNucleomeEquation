# 4DNucleomeEquation

4DNucleomeEquation simulation analysis and statistical inference of gene expression regulated by enhancer-promoter interaction across space and time described in the paper "4D nucleome equation predicts gene expression controlled by long-range enhancer-promoter interaction".

## **Environment setup**

First let’s set up the environment to run the simulations. The simulated system is complex and requires parallel computing. We can put the simulated code on the server for computing, or test a small demo on the local machine. First of all, download the repository as a zip file clicking the green button in the top right of this page to your local machine.

### **Download the scripts**

First of all, download the repository as a zip file clicking the green button in the top right of this page to your local machine.

### **Requirements**

The software prerequisites for the package to work are Matlab R2018a or later version of Matlab. You can start parallel pool (parpool) using the command `parpool`. Then, you can obtain the number of workers.

```matlab
%% current parallel pool
% By default, all available cores are used. If you want to specify the number of workers, 
% make sure the specified number is less than or equal to the total number of cores.
if isempty(gcp('nocreate'))
    numCores = feature('numcores'); % Get the number of available cores in the system
    parpool(numCores); % `numCores` is the number of workers you want to use, default is the number of cores in your computer
end
```

### **Loading the scripts**

The next step is to manually add the scripts in the same folder.

## Directories

There are two folders, and each folder contains several .m files. Next, we give detailed explanation each script:

### SimulationAnalysis

#### testdG.m, testKEP.m

The script can simulate and analyze different values of enhancer-promoter (E-P) genome distance (or E-P interaction strength). We simulate multiple times in parallel and store the simulated data in the folder generated according to the parameters.

#### ParametersBurst.m

This function generates a parameter structure to be passed to the simulation framework. The parameters are saved fields of a structure called `params`. We need to pass in the user defined parameters and set some default parameter values.

#### SimulateBurst.m

This .m files is core simulation framework. This script includes pre-allocating variable, initializing model, updating model and saving results.

#### InitializeConnectivityMatrix.m

This function is used to initialize generalized Rouse model, and was designed for initialized params structure input format. We supply this file to generate connectivity matrix to represent the connection of linear monomers. This file should not be used for other purposes.

#### AnalyseBurst.m

This .m file is the core analysis process. It contains two parts: Statistical analysis of simulation data and theoretical calculation of statistical indicators. We can obtain the results numerically and theoretically containing: mean mRNA level, CV, etc.

#### Poissbeta.m

This .m file can theoretically calculate the probability density function of ON-OFF model.

#### drawMeanVardG.m, drawMeanVarKEP.m

This script can plot the mean and CV changes of E-P genome distance (or E-P interaction strength) based on simulation results.

#### drawBC.m

This script theoretically computes the bimodal coefficient and determine the peak number boundaries.

#### drawDGpeakinformation.m, drawKEPpeakinformation.m

This script theoretically calculate the peak information for E-P genome distance (E-P interaction strength) under a fixed E-P interaction strength (E-P genome distance). We show the evolutionary process of the peak numbers and peak probabilities.

### StatisticalInference

#### dataToFit.mat, dGvsEncounterProb.mat

`dataToFit.mat`  provides the experimental data we are going to fit. The data includes the E-P genome distance, mRNA distribution, E-P encounter prob and mRNA mean level for different cell lines.  `dGvsEncounterProb.mat`  provides the encounter prob and corresponding E-P genomic distance.

#### fitmRNADistribution.m

This script  the master program for estimating the  E-P interaction dynamics and gene expression dynamics parameters from smFISH data and E-P genomic distance data.

#### fitEPEncounterProb.m

This script compares theoretical results based on statistically inferred parameters with experimental data. We show the relationship between E-P genome distance and contact probability. Mean eGFP mRNA level plotted against contact probability between the ectopic Sox2 promoter and SCR insertions. CV of eGFP level against contact probability.

#### CEInferenceBurst.m

This script is the core code of statical inference. We use the minimum cross entropy method to estimate the parameters.

#### CalculateTheoryProb.m

This .m file theoretical calculate of mRNA distribution.

## **Running the scripts**

先运行SimulationAnalysis文件夹中的`testdG.m`和`testKEP.m`（这俩不分先后），然后依次运行StatisticalInference文件夹中的`fitmRNADistribution.m`和`fitEPEncounterProb.m`。

## 与论文的匹配

### 4D核体方程

论文描述：论文提出了一个理论框架，即4D核体方程，整合了染色质上的E-P相互作用和基因转录的生化反应。
代码实现：SimulateBurst.m可能包含模拟4D核体方程的代码，负责模拟基因表达动态。

### 染色质动态模型

论文描述：染色质被建模为一个聚合物，通过哈密顿量描述其动态。
代码实现：InitializeConnectivityMatrix.m可能用于初始化广义Rouse模型，代表线性单体之间的连接。

### 基因表达动态模型

论文描述：基因表达过程被建模为一个两态模型，包括活跃（ON）和沉默（OFF）状态。
代码实现：ParametersBurst.m可能用于生成参数结构，这些参数可能包括基因表达的激活率和转录率。

### 链接函数

论文描述：链接函数将上游的E-P拓扑结构与下游的基因表达联系起来。
代码实现：SimulateBurst.m可能包含实现链接函数的代码，模拟E-P空间距离对转录相关反应率的影响。

### mRNA分布模型

论文描述：给定E-P空间距离时，mRNA丰度的稳态概率分布是泊松-贝塔分布。
代码实现：Poissbeta.m可能用于理论上计算ON-OFF模型的概率密度函数。

### 时间尺度分离方法

论文描述：通过时间尺度分离方法，作者推导出了mRNA的稳态分布。
代码实现：AnalyseBurst.m可能包含分析模拟数据和理论计算统计指标的代码，如平均mRNA水平和CV。

## 实验数据与代码文件对应

### 实验数据

论文描述：论文中使用了小鼠胚胎干细胞的smRNA-FISH数据和E-P基因组距离数据。
代码实现：dataToFit.mat提供了实验数据，用于拟合模型。

### 模型选择和参数推断

论文描述：4D核体方程允许模型选择和参数推断。
代码实现：fitmRNADistribution.m和fitEPEncounterProb.m用于从实验数据中估计E-P相互作用动态和基因表达动态参数。

### 双峰和多峰mRNA分布

论文描述：研究发现，长程E-P相互作用可以诱导双峰和三峰mRNA分布。
代码实现：drawBC.m用于计算双峰系数并确定峰值数边界；drawDGpeakinformation.m和drawKEPpeakinformation.m用于计算在固定E-P相互作用强度或基因组距离下的峰值信息。

### 统计推断

论文描述：进一步的统计推断表明E-P相互作用更倾向于通过控制启动子激活和转录起始率来调节mRNA水平。
代码实现：CEInferenceBurst.m使用最小交叉熵方法估计参数；CalculateTheoryProb.m理论上计算mRNA分布。