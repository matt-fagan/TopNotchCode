lrat_sp1 = zeros(10,1);
omega2 = zeros(10,1);
sp1Vscore = zeros(10,1);
sp1Vtest = zeros(10,1);

lrat_sp5 = zeros(10,1);
omega2 = zeros(10,1);
sp5Vscore = zeros(10,1);
sp5Vtest = zeros(10,1);

lrat_vix = zeros(10,1);
omega2 = zeros(10,1);
vixVscore = zeros(10,1);
vixVtest = zeros(10,1);

lrat_vxo = zeros(10,1);
omega2 = zeros(10,1);
vxoVscore = zeros(10,1);
vxoVtest = zeros(10,1);

for kay = 1:10

    lrat_vxo(kay,1) = LL_vxo_Me(9) - LL_vxo_Me(kay);
    omega2(kay,1) = var((LLs_vxo(9,:) - LLs_vxo(kay,:)));
    
    vxoVscore(kay,1) = sqrt(length(LLs_vxo(9,:)))*(lrat_vxo(kay)/length(LLs_vxo(kay,:)));
    vxoVtest(kay,1) = normcdf(vxoVscore(kay),0,sqrt(omega2(kay)));



    lrat_vix(kay,1) = LL_vix_Me(10) - LL_vix_Me(kay);
    omega2(kay,1) = var((LLs_vix(10,:) - LLs_vix(kay,:)));
    
    vixVscore(kay,1) = sqrt(length(LLs_vix(10,:)))*(lrat_vix(kay)/length(LLs_vix(kay,:)));
    vixVtest(kay,1) = normcdf(vixVscore(kay),0,sqrt(omega2(kay)));


    lrat_sp5(kay,1) = LL_sp5_Me(7) - LL_sp5_Me(kay);
    omega2(kay,1) = var((LLs_sp5(7,:) - LLs_sp5(kay,:)));
    
    sp5Vscore(kay,1) = sqrt(length(LLs_sp5(7,:)))*(lrat_sp5(kay)/length(LLs_sp5(kay,:)));
    sp5Vtest(kay,1) = normcdf(sp5Vscore(kay),0,sqrt(omega2(kay)));


    lrat_sp1(kay,1) = LL_sp1_Me(9) - LL_sp1_Me(kay);
    omega2(kay,1) = var((LLs_sp1(9,:) - LLs_sp1(kay,:)));
    
    sp1Vscore(kay,1) = sqrt(length(LLs_sp1(9,:)))*(lrat_sp1(kay)/length(LLs_sp1(kay,:)));
    sp1Vtest(kay,1) = normcdf(sp1Vscore(kay),0,sqrt(omega2(kay)));
    
end

Model_Selection_Vtests = [1-sp5Vtest, 1-sp1Vtest, 1-vixVtest, 1-vxoVtest]