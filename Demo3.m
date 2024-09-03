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
myfun = @(x)x-4*x^2+0.2*x^3 +2*x^4;

%% Set parameters
X0 = -10; % Initial point
n = numel(X0);
Method = 'Method2'; % Choose Method2
MaxIter = 100;      % Maximum number of iterations
Lambda = 0.01*eye(MaxIter*n); % Step size parameter

%% Find the first local minimum point
R = 0.01*eye(n);  
[X1, FVAL1, EXITFLAG1, GRAD1] = OMBOC(myfun, X0, R, Method, MaxIter, Lambda);

%% Find the second local minimum point
R = 100*eye(n);  
[X2, FVAL2, EXITFLAG2, GRAD2] = OMBOC(myfun, X0, R, Method, MaxIter, Lambda);