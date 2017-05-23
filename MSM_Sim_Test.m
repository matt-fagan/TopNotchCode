
%%%%%% Parameter Values %%%%%%%%%

b = 3 ;
m0 = 1.5;
gamma_kbar = 0.95; 
sigma = .25;

kbar=5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = zeros((2^kbar),(2^kbar));
    for i =0:2^kbar-1       
        for j = i:(2^kbar-1)-i  
            A(i+1,j+1) = bitxor(i,j);
        end
    end
    
    gamma = zeros(kbar,1);                          
    gamma(1) = 1-(1-gamma_kbar)^(1/(b^(kbar-1)));
    for i = 2:(kbar)
        gamma(i,1) = 1-(1-gamma(1,1))^(b^(i-1));
    end
    gamma = gamma*0.5;
    gamma(:,2)=gamma(:,1);
    gamma(:,1) = 1 - gamma(:,1);  
    kbar1 = kbar +1;
    kbar2 = 2^kbar;
    prob = ones(kbar2,1);    
    
    for i=0:2^kbar-1    %Works out probability associated with each XOR number
        for m = 1:kbar  
            prob(i+1,1) = prob(i+1,1) * gamma(kbar1-m, (bitget(i,m)+1));
    end
    end
    for i =0:2^(kbar-1)-1   %Copies probabilities to the transition matrix
        for j = i:(2^(kbar-1)-1)  
            A(kbar2-i,j+1) = prob(kbar2-A(i+1,j+1),1);%Copies each probability to the other 8 symmetrical locations
            A(kbar2-j,i+1) =  A(kbar2-i,j+1);
            A(j+1,kbar2-i) =  A(kbar2-i,j+1);
            A(i+1,kbar2-j) =  A(kbar2-i,j+1);    
            A(i+1,j+1) = prob(A(i+1,j+1)+1,1);
            A(j+1,i+1) = A(i+1,j+1);
            A(kbar2-j,kbar2-i) = A(i+1,j+1);
            A(kbar2-i,kbar2-j) = A(i+1,j+1);
        end
    end 


% calculate all of the possible volatility states 
    m1=2-m0;
    kbar2 = 2^kbar;
    g_m1 = [0:(kbar2-1)];

    for i = 1:(kbar2)
        g=1;
        for j = 0:(kbar-1)       
            if(bitand(g_m1(i),(2^j))~=0)    %
                g=g*m1;
            else g=g*m0;
            end
        end
        g_m1(i)=g;
    end
    
    g_m=sqrt(g_m1);
    
 numsims = 400;
 
 parameter_matrix = zeros(4,400);
 LL_matrix = zeros(8,400);
 
% LL_tester = zeros(8,1);
 
for simulation = 141:numsims    
    
    
    %%%%% Simulated Series %%%%%%
    
    sample_state = 2^kbar; %randsample(1:32,1);
    
    sims = 7500;

  sim_xt = []; 
  sample_state_vect = sample_state;
  sample_gm = []; %zeros(sims+1,1);
  epsilon_vect = []; %zeros(sims+1,1);
  sample_gm_vect = []; % zeros(sims+1,1);

  for sim=1:sims
    
      sample_state = randsample(1:2^kbar, 1, true, A(sample_state,:));
    
    epsilon = normrnd(0,1,1,1);
    
    sample_gm = g_m(sample_state);
    
    sim_xt = cat(1, sim_xt , sigma*sample_gm*epsilon );
    sample_state_vect = cat(1, sample_state_vect,sample_state);
    epsilon_vect = cat(1, epsilon_vect, epsilon);
    sample_gm_vect=cat(1, sample_gm_vect, sample_gm);

  end
  
  
     for kay = 1:8
       [parameters,LL] = MSM_modified(sim_xt,kay)
       LL_matrix(kay,simulation) = LL;
        if kay == kbar
            parameter_matrix(:,simulation) = parameters;
        end
     end
% 
%   sim_xt_ts = timeseries(sim_xt);
%   plot(sim_xt_ts) %,type="l", main=paste("initial state =",sample.state.vect[1]))
  

end

