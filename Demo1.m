% Note: Before using this example, you need to first download the CasADi toolkit
% and add it to your MATLAB path.

% How to Set Up the CasADi Path
%
% 1. Visit the CasADi official website [CasADi Download Page](https://github.com/casadi/casadi/releases).
% 2. Download the CasADi release suitable for your MATLAB version (e.g., casadi-3.6.6-windows64-matlab2018b.zip).
% 3. Extract the downloaded file.
% 4. Open MATLAB and add the extracted CasADi folder path to the MATLAB path using the following command:
%    addpath('path/to/casadi-folder');
%    savepath; % Save the path for future use
%
% If you encounter any issues, please refer to the CasADi documentation or contact technical support.

%% Objective function
myfun = @(x) x^2;

%% Set parameters
X0 = 10.0; 
n = numel(X0);
R = 100*eye(n);  
MaxIter = 100;      % Maximum number of iterations
U = 0.01*ones(n,MaxIter);

%% Method1
Method = 'Method1'; % Choose Method1
Lambda = 1e-4;      % Step size parameter, scalar

% Call OMBOC function
[X1, FVAL1, EXITFLAG1, GRAD1] = OMBOC(myfun, X0, R, Method, MaxIter, Lambda,U);

%% Method2
Method = 'Method2'; % Choose Method2
Lambda = 0.01*eye(MaxIter*n); % Step size parameter, matrix

% Call OMBOC function
[X2, FVAL2, EXITFLAG2, GRAD2] = OMBOC(myfun, X0, R, Method, MaxIter, Lambda,U);
