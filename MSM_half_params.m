SP500dathalf = SP500dat(1:(length(SP500dat)/2));
SP100dathalf = SP100dat(1:((length(SP100dat)/2)+.5));
VIXdathalf = VIXdat(1:((length(VIXdat)/2)+.5));
VXOdathalf = VXOdat(1:((length(VXOdat)/2)+.5));

startingvals=[];

kay=7
    
A_template = zeros((2^kay),(2^kay));
    for i =0:2^kay-1;       
        for j = i:(2^kay-1)-i;  
            A_template(i+1,j+1) = bitxor(i,j);
        end
    end
    
    %%%%% SP500 Run %%%%%

[startingvals_sp5,sp5LLs,sp5ordered_parameters] = ...
    MSM_starting_values(SP500dathalf,startingvals,kay,A_template)

startingvals_sp5 = [startingvals_sp5(1,1),startingvals_sp5(3,1),startingvals_sp5(2,1),startingvals_sp5(4,1)]
    
[data, kay, startingvals_sp5, LB, UB, options] = ...
            MSM_parameter_check(SP500dathalf, kay, startingvals_sp5);

[sp5parameters,LL,exitflag,output]=fmincon('MSM_likelihood',startingvals_sp5,...
    [],[],[],[],LB,UB,[],options,kay,data,A_template);

%[sp5LL,sp5LLs] = MSM_likelihood(sp5parameters,kay,data,A_template);


%%%%%%% SP100

kay = 9

A_template = zeros((2^kay),(2^kay));
    for i =0:2^kay-1;       
        for j = i:(2^kay-1)-i;  
            A_template(i+1,j+1) = bitxor(i,j);
        end
end

[startingvals_sp1,sp1LLs,sp1ordered_parameters] = ...
    MSM_starting_values(SP100dathalf,startingvals,kay,A_template)

startingvals_sp1 = [startingvals_sp1(1,1),startingvals_sp1(3,1),startingvals_sp1(2,1),startingvals_sp1(4,1)]

[sp1data, kay, startingvals_sp1, LB, UB, options] = ...
            MSM_parameter_check(SP100dathalf, kay, startingvals_sp1);

[sp1parameters,LL,exitflag,output]=fmincon('MSM_likelihood',startingvals_sp1,...
    [],[],[],[],LB,UB,[],options,kay,sp1data,A_template);

%[sp1LL,sp1LLs] = MSM_likelihood(sp1parameters,kay,sp1data,A_template);


%%%%%%%%% VIX

kay = 10

A_template = zeros((2^kay),(2^kay));
    for i =0:2^kay-1;       
        for j = i:(2^kay-1)-i;  
            A_template(i+1,j+1) = bitxor(i,j);
        end
end

[startingvals_vix,vixLLs,vix_ordered_parameters] = ...
    MSM_starting_values(VIXdathalf,startingvals,kay,A_template)

startingvals_vix = [startingvals_vix(1,1),startingvals_vix(3,1),startingvals_vix(2,1),startingvals_vix(4,1)]
    
[vixdata, kay, startingvals_vix, LB, UB, options] = ...
            MSM_parameter_check(VIXdathalf, kay, startingvals_vix);

[vixparameters,LL,exitflag,output]=fmincon('MSM_likelihood',startingvals_vix,...
    [],[],[],[],LB,UB,[],options,kay,vixdata,A_template);

%[vixLL,vixLLs] = MSM_likelihood(vixparameters,kay,vixdata,A_template);


%%%%%%%%%%%%%%% VXO

kay = 9

A_template = zeros((2^kay),(2^kay));
    for i =0:2^kay-1;       
        for j = i:(2^kay-1)-i;  
            A_template(i+1,j+1) = bitxor(i,j);
        end
end

[startingvals_vxo,vxoLLs,vxo_ordered_parameters] = ...
    MSM_starting_values(VXOdathalf,startingvals,kay,A_template)

startingvals_vxo = [startingvals_vxo(1,1),startingvals_vxo(3,1),startingvals_vxo(2,1),startingvals_vxo(4,1)]
    
[vxodata, kay, startingvals_vxo, LB, UB, options] = ...
            MSM_parameter_check(VXOdathalf, kay, startingvals_vxo);

[vxoparameters,LL,exitflag,output]=fmincon('MSM_likelihood',startingvals_vxo,...
    [],[],[],[],LB,UB,[],options,kay,vxodata,A_template);

%[vxoLL,vxoLLs] = MSM_likelihood(vxoparameters,kay,vxodata,A_template);
