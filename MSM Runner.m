% data2 = csvread('SP500_Full_Real.csv');
% 
% data3=xlsread('data_demo.xls');

SP500dat = csvread('SP500_VIX.csv',1,1);
VIXdat = csvread('VIX_Daily_Returns.csv',1,1);

SP100dat = csvread('sp100.csv',1,1);
VXOdat = csvread('VXO.csv',1,1);


LL_sp500 = zeros(10,1);
LL_VIX = zeros(10,1);
LL_sp100 = zeros(10,1);
LL_VXO = zeros(10,1);

parameters_sp500 = zeros(10,4);
parameters_sp100 = zeros(10,4);
parameters_vix = zeros(10,4);
parameters_vxo = zeros(10,4);

LLs_sp5 = repmat(zeros(length(SP500dat)),10,1);
LLs_sp1 = repmat(zeros(length(SP100dat)),10,1);
LLs_vxo = repmat(zeros(length(VXOdat)),10,1);
LLs_vix = repmat(zeros(length(VIXdat)),10,1);


for kay = 8:10
%kay = 10
        [parameters,LL,LLs] = MSM_modified(SP500dat,kay);
        parameters_sp500(kay,:) = parameters;
        LLs_sp5(kay,:) = LLs;
        LL_sp500(kay,1) = LL;
        
        [parameters,LL,LLs] = MSM_modified(VIXdat,kay);
        parameters_vix(kay,:) = parameters;
        LLs_vix(kay,:) = LLs;
        LL_VIX(kay,1) = LL;
        
        [parameters,LL,LLs] = MSM_modified(SP100dat,kay);
        parameters_sp100(kay,:) = parameters;
        LLs_sp1(kay,:) = LLs;
        LL_sp100(kay,1) = LL;
        
        [parameters,LL,LLs] = MSM_modified(VXOdat,kay);
        parameters_vxo(kay,:) = parameters;
        LLs_vxo(kay,:) = LLs;
        LL_VXO(kay,1) = LL;
        
end



%%%%%%%% Now, Test the MLE against the Sim Series %%%%%%%%%%
  
%startingvals = [b, m0, gamma_kbar, sigma];
% startingvals = [];
% 
% %params_mat = zeros(4,1);
% 
%     

%params_mat(:) = parameters

%    for kay = 1:8
%        [parameters,LL] = MSM_modified(sim_xt,8)
%        LL_tester(kay,1) = -LL
%         if kay == kbar
%             parameter_matrix(:,simulation) = parameters
%         end
%    end
    
%    LL_tester

%-sum(LLs)