function [X,FVAL,EXITFLAG,GRAD] = OMBOC(FUN,X0,R,Method,MaxIter,Lambda,U,TolX,ShowGraph,varargin)
% OMBOC finds an unconstrained minimum of a non-convex function of several variables .
%   OMBOC attempts to solve problems of the form:
%    min F(X)
%     X
%
%   X = OBOMBOCOMBOC(FUN,X0,R,Method,MaxIter,Lambda,U,TolX,ShowGraph,varargin)
%   starts at X0 and finds a minimum X to the function FUN. FUN a accepts input 
%   X and returns a scalar function value F evaluated at X . X0 may be a
%   scalar, vector, or matrix. R is the positive definite control weight
%   matrix.
%   Method is a flag to select the optimazation methods listed below.
%     'Method1'  Algorithm in [1] .(default) 
%     'Method2'  Algorithm in (9)-(10) of [2](R=Lambda). If you want to
%                descrease the computation time (especially in 
%                high-dimensional cases), it is recommended this method.
%   Lambda is an iteration step size in Method1 or Method2.
%   MaxIter is the maximum number of iterations allowed. U is the guess of the
%   initial control sequence. TolX is the termination tolerance on x (1e-3 for 
%   default). ShowGraph is a {'on'}|{'off'} flag to show the iteration process 
%   or not.
%
%   [X,FVAL] = OMBOC(FUN,X0,...) returns the value of the objective 
%   function FUN at the solution X .
%
%   [X,FVAL,EXITFLAG] = OMBOC(FUN,X0,...) returns an EXITFLAG that 
%   describes the exit condition of OMBOC . Possible values of EXITFLAG 
%   and the corresponding exit conditions are listed below .
%   
%   All algorithms:
%     1  First order optimality conditions satisfied to the specified 
%        tolerance .
%     0  Maximum number of function evaluations or iterations reached .
%    -1  Preseved .
%
%   [X,FVAL,EXITFLAG,GRAD] = OMBOC(FUN,X0,...) returns the value of the
%   gradient of FUN at the solution X .
%
%   Here, myfun is the objective function that you want to minimize:
%
%       myfun = @(x)x-4*x^2+0.2*x^3 +2*x^4;
%
%
%   Finally, pass these anonymous functions to OMBOC:
%
%        x = OMBOC(@(x)myfun(x),x0,...,TolX,)
%                
%   See details in Demos .
%
%   OMBOC uses the Huanshui Zhang and Hongxia Wang optimization method based
%   on optimal oontrol (described in [1] and [2]), and is coded in MATLAB 2008a 
%   and tested in subsequent versions of MATLAB .
%
%   References:
%   [1] Yeming Xu, Ziyuan Guo, Hongxia Wang, Huanshui Zhang, "Optimization Method
%   Based on Optimal Control," 2023, https://arxiv.org/abs/2309.05280 .
%	[2] Huanshui Zhang, Hongxia Wang, "Optimization Methods Rooted in
%	Optimal Control," 2023, https://arxiv.org/abs/2312.01334.

%   See also FMINCON, FMINUNC, FMINBND, FMINSEARCH, @, FUNCTION_HANDLE .

%   Copyright 2024-2034
%   $Revision: 1.0.0.0.0.1 $  $Date: 2024/09/02 20:11:42 $ $Coder: Yeming Xu&Chuanzhi Lv&Kai Peng$
%   $Tester: Yeming Xu$ $Reviewer: Hongxia Wang$


if nargin <2
    error('OMBOC:NotEnoughInputs',...
        'OMBOC requires at least TWO input arguments');
end

% Check for non-double inputs
X0 = X0(:);
if ~isa(X0,'double')
    error('OMBOC:NonDoubleInput', ...
        'OMBOC only accepts inputs of data type double . ')
end
n = numel(X0);

% Check if R is a positive definite matrix
if  nargin < 3 || isempty(R)
    R = eye(n);
else
    nR = size(R);
    if length(nR) > 2 || length(nR) < 2 || nR (1) ~= nR(2) || nR (1) ~= n || 1 ~= check_matrix_definiteness(R)
        error('OMBOC:OptRNotPositiveDefiniteMatrix',...
            'Option ''R'' must be an n-by-n positive definite matrix where n is equal to deminsion of ''X0''')

    end
end

% Check if Method is selected properly
if nargin < 4 || isempty(Method)
    Method = 'Method1';  % Default to Method1 if not provided
elseif ~strcmpi(Method, 'Method1') && ~strcmpi(Method, 'Method2')
    error('OCP:OptMethodInputError',...
        'Option ''Method'' must be selected as either ''Method1'' or ''Method2''.');
end

% Check N is a positive integer value
if  nargin < 5 || isempty(MaxIter)
    MaxIter = 50 * n;
elseif ~isscalar(MaxIter) || MaxIter < 1
    error('OMBOC:OptMaxIterNotInteger',...
        'Option ''MaxIter'' must be a positive integer value if not the default . ')
end

% Check Lambda for Method1 or Method2
if nargin < 6 || isempty(Lambda)
    if strcmpi(Method, 'Method1')
        Lambda = 1e-6;  % Set Lambda for Method1
    elseif strcmpi(Method, 'Method2')
        Lambda = 0.1 * eye(n * MaxIter);  % Set Lambda for Method2
    end
elseif (strcmpi(Method, 'Method1') && (~isscalar(Lambda) || Lambda > 1)) || ...
       (strcmpi(Method, 'Method2') && (~ismatrix(Lambda) || any(size(Lambda) ~= [n * MaxIter, n * MaxIter]) || ...
       ~isequal(Lambda, Lambda') || any(eig(Lambda) < 0)))
    error('OMBOC:OptLambdaNotValid',...
        ['Option ''Lambda'' must be a positive scalar less than 1 for Method1 or ', ...
        'a positive semi-definite matrix of size n*MaxIter-by-n*MaxIter for Method2.']);
end

% Check U is an n-by-MaxIter matrix
if  nargin < 7 || isempty(U)
    U = ones(n,MaxIter);
elseif ~ismatrix(U) || size (U, 1) ~= n || size (U, 2) ~= MaxIter
    error('OMBOC:OptUNotControlSequence',...
        'Option ''U'' must be an n-by-MaxIter matrix if not the default . ')
end

% Check if TolX is a positive scalar
if  nargin < 8 || isempty(TolX)
    TolX = 1e-3;
elseif ~isscalar(TolX) || TolX <=0
    error('OMBOC:OptTolXNotPositiveScalar',...
        'Option ''TolX'' must be a positive scalar if not the default . ')
end

% Check if ShowGraph is on
if  nargin < 9 || isempty(ShowGraph)
    ShowGraph = 'on';
elseif ~strcmpi(ShowGraph,'on') && ~strcmpi(ShowGraph,'off')
    error('OMBOC:OptShowGraphInputError',...
        'Option ''ShowGraph'' must be selected in {''on''} | ''off'' . ')
end

% Check if Fun is scalar
if ~isempty(FUN),
    F0 = checkfun(X0,FUN,varargin,'1');
else
    error('OMBOC:OptFUNInputError','Option ''FUN'' can''t be []');
end

% ---------------------------------------------------
% Start the Iteration Precess
% ------------------------1--------------------------
disp('Iteration is started . ')
X = X0;
hf = 0;
EXITFLAG = 1;
iterCount = 0;  % Initialize iteration counter
maxIterations = 10^6;
if strcmpi(ShowGraph,'on'),callAllOptimPlotFcns('init'),end

switch lower(Method)
    case 'method1'      % Algorithm in [1]
        
        x=zeros(n,MaxIter);
        x(:, 1) = X0;
        H_u(:,1)=ones(n,1);
        % Define the gradient function using GradFcn
        gradStruct = Gradfcn(FUN,n);          % Call GradFcn to get the structure
        f_gradients = gradStruct.F_gradients; % Extract the gradient function
        
        % Main optimization loop
        while norm(H_u, 'fro') > TolX
            
             iterCount = iterCount + 1;  % Update iteration counter
             
            % Check for NaN or Inf in U, x, and H_u
            if any(isnan(x), 'all') || any(isinf(x), 'all')
                error('OMBOC:xContainsNaN', 'Iteration trajectory x generates NaN or Inf values, please decrease Lambda');
            end
            
            % Solving Forward-Backward Differential Equations (FBDEs)
            for i = 1:MaxIter
                x(:, i+1) = x(:, i) + U(:, i);
            end
            
            % Compute gradient at the final state
            l(:, MaxIter+1) = full(f_gradients(x(:, MaxIter+1)));  % Convert CasADi result to full matrix
            
            % Backward computation for gradients
            for i = MaxIter:-1:1
                l(:, i) = l(:, i+1) + full(f_gradients(x(:, i)));  % Convert CasADi result to full matrix
            end

            % Calculate H_u
            for i = 1:MaxIter
                H_u(:, i) = R * U(:, i) + l(:, i+1);
            end
            
            % Update control inputs U using gradient descent
            deltau = -Lambda * H_u;
            U = U + deltau;

             % Check for maximum iterations
            if iterCount >= maxIterations
                disp('Maximum number of iterations reached.');
                EXITFLAG = 0;  % Set EXITFLAG to 0 for maximum iterations
                break;
            end
            
        end
       
    case 'method2'      % Algorithm in [2]
        
        % Using casadi to get gradient and heesians in Method2
        param=reshape(U, MaxIter * n, 1);
        result = a_high_dimension_problem(MaxIter, n, FUN);

        while true
            grad = result.Loss_gradients(param, X, R);
            hessi = result.Loss_hessians_Inv(param, X, R, Lambda);  % Note: Lambda as a matrix
            g = hessi * grad;
       
            for j = 2:iterCount
                g = hessi * (grad + Lambda * g);
            end
            
            % Check for NaN or Inf in g
            if any(isnan(full(g)), 'all') || any(isinf(full(g)), 'all')
                error('OMBOC:gContainsNaN', 'Iteration trajectory x contains NaN or Inf values. Please adjust Lambda.');
            end
         
            % Compute the Frobenius norm of g using CasADi functions
            if full(norm(grad, 'fro')) < TolX
                break;
            end
          
            param = param - g;
            iterCount = iterCount + 1;
            
            % Check for maximum iterations
            if iterCount >= maxIterations
                disp('Maximum number of iterations reached.');
                EXITFLAG = 0;  % Set EXITFLAG to 0 for maximum iterations
                break;
            end
            
        end
        
       %% Calculate State and Performance Metrics 
        deal_param = reshape(param, n, MaxIter);
        u = full(deal_param);
        x = zeros(n, MaxIter+1);
        x(:, 1) = X0;
        for k = 1:MaxIter
            x(:, k+1) = x(:, k) + u(:, k);
        end
      
    otherwise
        error('OCP:OptMethodError','Unknown method.')
end

In = eye(n); tmpX = X; tmpF = F0;
% Assuming GradFcn is the function structure obtained from the GradFcn function
gradStruct = Gradfcn(FUN,n);  % Call GradFcn to get the structure containing CasADi functions
f_gradients = gradStruct.F_gradients; % Extract the gradient function

for j = 1:MaxIter
    X = x(:, j);
   % If graphical output is enabled, call plotting functions
   if ShowGraph,callAllOptimPlotFcns('iter');end
end

FVAL = feval(FUN, X, varargin{:});
GRAD = full(f_gradients(X));  % Use CasADi gradient function for final gradient

% Handle graphical output
if ShowGraph && ishandle(hf)
    figure(hf);
    % Any additional plotting code
end


disp('Iteration is Finished . ')
% pause
% if ishandle(hf),close(hf),end

% ----------------------------
% Nested Fun for PLOT
% ============================
function callAllOptimPlotFcns(state)
    switch state
        case 'init'
            hf = figure('numbertitle','off','Name',['Iteration Process'],'Visible','on');
            clf
            subplot(211),hold on,box on
            ylabel('F')
            subplot(212),hold on,box on
            xlabel('Iteration Number')
            ylabel('X')
            ShowGraph = true;
        case 'iter'
            F0 = feval(FUN, X, varargin{:});
            set(0,'CurrentFigure',hf)
            subplot(211),hold on,box on
            plot([j,j+1],[tmpF,F0],'r-d','MarkerEdgeColor','r',...
                'MarkerFaceColor','r',...
                'MarkerSize',5)
            subplot(212)
            plot([j,j+1],[tmpX,X]','-d','MarkerSize',5)
            tmpX = X; tmpF = F0;
        otherwise
            error('OMBOC:callAllOptimPlotFcns:state','Unkown state . '); 
    end
end % end of callAllOptimPlotFcns

end % end of OMBOC

function matrix_type = check_matrix_definiteness(A)
    % input:
    % A - the matrix to be checked if positive definite
    % outputs:
    % matrix_type -: 1-'positive_definite',0- 'not positive_definite'

    if ~isequal(A, A')
        error('Matrix must be square . ');
    end
    eigenvalues = eig(A);

    num_positive = sum(eigenvalues > 0);
    n = length(eigenvalues);

    if num_positive == n
        matrix_type = 1;
    else
        matrix_type = 0;
    end
end

function f = checkfun(x,userfcn,varargin,flag)
    % CHECKFUN checks for complex or NaN results from userfcn .

    f = feval(userfcn, x, varargin{:});
    % Note: we do not check for Inf as OMBOC handles it naturally .
    if isnan(f)
        error('OMBOC:checkfun:NaNFval', ...
            'User function ''% s'' returned NaN when evaluated;\n OMBOC cannot continue . ', ...
            localChar(userfcn));  
    elseif ~isreal(f)
        error('OMBOC:checkfun:ComplexFval', ...
            'User function ''% s'' returned a complex value when evaluated;\n OMBOC cannot continue . ', ...
            localChar(userfcn));  
    end
    switch flag
        case '1'
            if ~isscalar(f),
                error('OMBOC:checkfun:ScalarFval', ...
                    'User function ''% s'' returned a non-scalar value when evaluated;\n OMBOC cannot continue . ', ...
                    localChar(userfcn)); 
            end

        case '2'
            m = size(f);
            if length (m) ~= 2 || length (f(:)) ~= length(x(:)) || numel (f)~= size(f,1)
                error('OMBOC:checkfun:VectorGradFval', ...
                    'User function ''% s'' did not return a column vector with dimension not matching x when evaluated;\n OMBOC cannot continue . ', ...
                    localChar(userfcn)); 
            end

        case '3'
            n = length(x(:));
            m = size(f);
            if length (m) ~= 2 || m (1) ~= n || m (2) ~= n
                error('OMBOC:checkfun:ScalarFval', ...
                    'User function ''% s'' did not return a square matrix with dimension not matching x when evaluated;\n OMBOC cannot continue . ', ...
                    localChar(userfcn)); 
            end
        otherwise
    end
end
function strfcn = localChar(fcn)
    % Convert the fcn to a string for printing

    if ischar(fcn)
        strfcn = fcn;
    elseif isa(fcn,'inline')
        strfcn = char(fcn);
    elseif isa(fcn,'function_handle')
        strfcn = func2str(fcn);
    else
        try
            strfcn = char(fcn);
        catch
            strfcn = '(name not printable)';
        end
    end
end
function result = a_high_dimension_problem(N, dim, f)
    % Helper function to set up CasADi optimization problem
    % Inputs:
    %   N    - number of control steps
    %   dim  - dimensionality of the state space
    %   f    - function handle for the state-dependent part of the loss function
    %
    % Outputs:
    %   result - structure containing CasADi functions for gradients and Hessian inverses
    
    import casadi.*
    
    % Define symbolic variables
    r = SX.sym('r',dim,dim);
    RR = SX.sym('RR', dim*N, dim*N);  % Ensure RR is a matrix of the correct size
    UU = SX.sym('UU', dim, N);  % Control variables for each dimension
    XX = SX.sym('XX', dim, N+1);  % State variables for each dimension
    P = SX.sym('P', dim, 1);  % Parameterized initial state
    
    % State transition
    XX(:, 1) = P;
    for k = 1:N
        XX(:, k+1) = XX(:, k) + UU(:, k);
    end
    
    % Loss function
    Loss = 0.0;
    for k = 1:N
        sta = XX(:, k);
        con = UU(:, k);
        Loss = Loss + f(sta) + con' * r * con;
    end
    sta = XX(:, N+1);
    Loss = Loss + f(sta);
    
    % Compute gradients and Hessians
    Parameters = reshape(UU, N * dim, 1);
    gradients = gradient(Loss, Parameters);
    hessians = hessian(Loss, Parameters);
    Inverse = inv(RR + hessians);  % Ensure RR is defined and invertible
    
    % CasADi functions
    result.Loss_gradients = Function('Loss_gradients', {Parameters, P, r}, {gradients});
    result.Loss_hessians_Inv = Function('Loss_hessians_Inv', {Parameters, P, r, RR}, {Inverse});

end
function result = Gradfcn(f,n)
    % Helper function to set up CasADi optimization problem
    % Inputs:
    %   f    - function handle for the state-dependent part of the loss function
    %
    % Outputs:
    %   result - structure containing CasADi function for the gradient
    
    import casadi.*  % Import CasADi functions
        
    % Create symbolic variables
    x = SX.sym('x', n);
    
    % Define function handle for CasADi
    func = @(x) feval(f, x);  % Convert function handle to be used with CasADi
    
    % Convert MATLAB function handle to CasADi function
    f_casadi = Function('f', {x}, {func(x)});
    
    % Compute gradient of function 'f' with respect to 'x'
    gradients = gradient(f_casadi(x), x);
    
    % Create CasADi function for the gradient
    F_gradients = Function('gradients', {x}, {gradients});
   
    % Return structure containing CasADi function for the gradient
    result.F_gradients = F_gradients;
end


