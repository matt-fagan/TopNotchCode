

vxo_clarketest = zeros(10,1);
vix_clarketest = zeros(10,1);
sp1_clarketest = zeros(10,1);
sp5_clarketest = zeros(10,1);

Bsp5 = zeros(10,1);
Bsp1 = zeros(10,1);
Bvix = zeros(10,1);
Bvxo = zeros(10,1);

for kay = 1:10
    
    di_sp5 = LLs_sp5(7,:) - LLs_sp5(kay,:);
    Bsp5(kay) = length(di_sp5(di_sp5>0))
    nsp5 = length(di_sp5)
    
    sp5_clarketest(kay,1)=binopdf(Bsp5(kay),nsp5,.5);

    di_sp1 = LLs_sp1(9,:) - LLs_sp1(kay,:);
    Bsp1(kay) = length(di_sp1(di_sp1>0));
    nsp1 = length(di_sp1);
    
    sp1_clarketest(kay,1)=binopdf(Bsp1(kay),nsp1,.5);
    
    di_vix = LLs_vix(10,:) - LLs_vix(kay,:);
    Bvix(kay) = length(di_vix(di_vix>0));
    nvix = length(di_vix);
    
    vix_clarketest(kay,1)=binopdf(Bvix(kay),nvix,.5);
    
    di_vxo = LLs_vxo(9,:) - LLs_vxo(kay,:);
    Bvxo(kay) = length(di_vxo(di_vxo>0));
    nvxo = length(di_vxo);
    
    vxo_clarketest(kay,1)=binopdf(Bvxo(kay),nvxo,.5);
    
end

%Model_Selection_Clarketests = [sp5_clarketest, sp1_clarketest, vix_clarketest, vxo_clarketest]

%write(array2table(Model_Selection_Clarketests),'Model_Selection_Clarketests.csv');