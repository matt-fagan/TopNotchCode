function [parameters,LL,LLs] = MSM_modified(data,kbar,startingvals, options)

% -------------------------------------------------------------------------
%                   Markov Switching Multifractal (MSM)                   
%                      Maximum likelihood estimation                       
%                                  v1.0                                      
%                Copyright © 2010 Multifractal-finance.com                 
% -------------------------------------------------------------------------
%
% USAGE
%    [PARAMETERS] = MSM(DATA, K)
%    [PARAMETERS, LL, LLs, DIAGNOSTICS] = MSM(DATA, K, STARTING_VALUES, OPTIONS)
%
% INPUTS:
%    DATA                -   A column (or row) of mean zero data
%    KBAR                -   The number of frequency components
%    STARTINGVALS        -   [OPTIONAL] Starting values for optimization
%                            [b, m0, gamma_k, sigma]
%                               b       - (1,inf) 
%                               m0      - (1,2] 
%                               gamma_k - (0,1) 
%                               sigma   - [0,inf)
%    OPTIONS             -   {OPTIONAL] User provided options structure
%
% OUTPUTS:
%    PARAMETERS          -   A 4x1 row vector of parameters 
%                            [b, m0, gamma_k, sigma]
%    LL                  -   The log-likelihood at the optimum
%    LLs                 -   Individual daily log-likelihoods at optimum
%    diagnostics         -   Structure of optimization output information.
%                            Useful for checking convergence problems
%
% ASSOCIATED FILES:
%    MSM_likelihood.m, MSM_parameter_check.m, MSM_starting_values.m
%
% REFERENCES:
%    [1] Calvet, L., Adlai Fisher (2004). "How to Forecast long-run 
%        volatility: regime-switching and the estimation of multifractal 
%        processes". Journal of Financial Econometrics 2: 49–83.
%    [2] Calvet, L., Adlai Fisher (2008). "Multifractal Volatility: Theory, 
%        Forecasting and Pricing". Elsevier - Academic Press..
% ------------------------------------------------------------------------- 

help MSM_modified

switch nargin
    case 2
        [data, kbar, startingvals, LB, UB, options] = ...
            MSM_parameter_check(data, kbar);
    case 3
        [data, kbar, startingvals, LB, UB, options] = ...
            MSM_parameter_check(data, kbar, startingvals);
    case 4
        [data, kbar, startingvals, LB, UB, options] = ...
            MSM_parameter_check(data, kbar, startingvals, options);
    otherwise
        error('Number of inputs must be between 2 and 4');
end

% Set up template for transition matrix to save initializing each iteration
A_template = T_mat_template(kbar);

% Get starting values (if none are supplied by user)
[startingvals,LLs,ordered_parameters] = ...
    MSM_starting_values(data,startingvals,kbar,A_template);

% % % Start seraching for minimum negative likelihood value
[parameters,LL,exitflag,output]=fmincon('MSM_likelihood',startingvals,...
    [],[],[],[],LB,UB,[],options,kbar,data,A_template);

% Compute log-likelihood and daily likelihood values
if nargout>1
    [LL,LLs] = MSM_likelihood(parameters,kbar,data,A_template);
    LL = -LL;
end

% Set up diagnostics for output
diagnostics.EXITFLAG=exitflag;
diagnostics.ITERATIONS=output.iterations;
diagnostics.FUNCCOUNT=output.funcCount;
diagnostics.MESSAGE=output.message;

disp(' ')
disp(' Estimated parameters ');
disp(' ')
disp(sprintf('  b = %3.4f ',parameters(1)))
disp(sprintf('  m0 = %3.4f ',parameters(2)))
disp(sprintf('  gamma_k = %3.4f ',parameters(3)))
disp(sprintf('  sigma = %3.4f ',parameters(4)))
disp(' ')
disp(sprintf('  Log-likelihood = %3.4f ',LL))

disp(' Diagnostic ')
disp(diagnostics);
disp(' ')
disp(' ')

    
%    PARAMETERS, LL, LLs, DIAGNOSTICS

% Save a template for the transition matrix
function A = T_mat_template(kbar)  
    A = zeros((2^kbar),(2^kbar));
    for i =0:2^kbar-1;       
        for j = i:(2^kbar-1)-i;  
            A(i+1,j+1) = bitxor(i,j);
        end
    end

