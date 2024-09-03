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
myfun = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;

%% Set parameters
X0 = [5;2];
n = numel(X0);
R = 0.1*eye(n);   

%% Method2
Method = 'Method2'; % Choose Method2, descrease the computation time
MaxIter = 10;      % Maximum number of iterations
Lambda = 0.01*eye(MaxIter*n); % Step size parameterm, matrix

% Call OMBOC function
[X, FVAL, EXITFLAG, GRAD] = OMBOC(myfun, X0, R, Method, MaxIter, Lambda);
