function [X,FVAL,EXITFLAG,GRAD,HESSIAN] = OCP(FUN,GradFcn,HessFcn,X0,R,MaxIter,TolX,Method,ShowGraph,varargin)
%OCP finds a unconstrained minimum of a function of several variables .
%   OCP attempts to solve problems of the form:
%    min F(X)
%     X
%
%   X = OCP(FUN,GradFcn,HessFcn,X0,R,MaxIter,TolX,Method,ShowGraph) starts
%   at X0 and finds a minimum X to the function FUN with user-supplied
%   Gradient GradFcn and Hessian FUN_HESS. FUN accepts input X and returns
%   a scalar function value F evaluated at X. X0 may be a scalar, vector,
%   or matrix. R is the positive scalar or positive definite matrix to
%   adjust iteration step size. MaxIter is the maximum number of iterations
%   allowed. TolX is the termination tolerance on x (1e-6 for default).
%   Method is a flag to select the optimazation methods listed below.
%     'Method1'  Algorithm in (9)-(10) of [1] or in (15) of
%                [2].(default) 
%     'Method2'  Algorithm described in (22) of [2] (M = R). If you want to
%                avoid the inverse of the Hessian matrix (especially in 
%                high-dimensional cases), it is recommended this method.
%   ShowGraph is a {'on'}|'off' flag to show the iteration process or not.
%
%   [X,FVAL] = OCP(FUN,X0,...) returns the value of the objective 
%   function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = OCP(FUN,X0,...) returns an EXITFLAG that 
%   describes the exit condition of OCP. Possible values of EXITFLAG 
%   and the corresponding exit conditions are listed below.
%   
%   All algorithms:
%     1  First order optimality conditions satisfied to the specified 
%         tolerance.
%     0  Maximum number of function evaluations or iterations reached.
%    -1  Preseved.
%
%   [X,FVAL,EXITFLAG,GRAD] = OCP(FUN,X0,...) returns the value of the
%   gradient of FUN at the solution X.
%
%   [X,FVAL,EXITFLAG,GRAD,HESSIAN] = OCP(FUN,X0,...) returns the value of
%   the Hessian at solution X.
%
%   FUN, GradFcn and HessFcn must be parameterized, you can use anonymous
%   functions to capture the problem-dependent parameters. Suppose you want
%   to minimize the objective given in the function myfun with Gradient
%   mygradfun and Hessian myhessfun, where these three functions are
%   parameterized by their second argument a1, a2, and a3, respectively.
%   Here myfun, mygradfun and myhessfun are M-file functions such as
%
%        function f = myfun(x,a1)
%        f = x(1)^3 + a1*x(2)^3;
%
%        function fg = mygradfun(x,a2)
%        fg = [3*x(1)^2;
%              3*x(2)^2 * a2];
%
%        function fh = myhessfun(x,a3)
%        fg = [6*x(1), 0;
%              0, 6*x(2) * a3];
%
%   To optimize for specific values of a1, a2 and a3, first assign the
%   values to these three parameters. Then create three one-argument
%   anonymous functions that capture the values of a1, a2 and a2, and call
%   myfun, mygradfun and myhessfun with two arguments. Finally, pass these
%   anonymous functions to OCP:
%
%        a1 = 2; a2 = 2; a3 = 2; % define parameters first
%
%        x = OCP(@(x)myfun(x,a1),@(x)mygradfun(x,a2),@(x)myhessfun(x,a3),...
%                	x0,R,MaxIter,TolX,Method)
%   Also, you may tranfer problem-dependent parameters to FUN, GradFcn and
%   HessFcn by varargin. See details in Demos.
%
%   OCP uses the Huanshui Zhang and Hongxia Wang optimization methods
%   rooted in optimal control( described in [1]), and is coded in MATLAB
%   2008a and tested in subsequent versions of MATLAB.
%
%   References:
%	[1] Huanshui Zhang, Hongxia Wang, "Optimization Methods Rooted in
%	Optimal Control," 2023, https://arxiv.org/abs/2312.01334.
%	[2] Hongxia Wang, Yeming Xu, Ziyuan Guo, and Huanshui Zhang,
%	"Optimization Algorithms with Superlinear Convergence Rate," 2024,
%	https://arxiv.org/abs/2403.11115.
%
%   See also FMINCON, FMINUNC, FMINBND, FMINSEARCH, @, FUNCTION_HANDLE.

%   Copyright 2024-2034
%   $Revision: 1.0.0.0.0.1 $  $Date: 2024/08/15 20:11:42 $ $Coder: Kai Peng$
%   $Tester: Ziyuan Guo$ $Reviewer: Hongxia Wang$


if nargin < 4
    error('OCP:NotEnoughInputs',...
        'OCP requires at least FOUR input arguments');
end

% Check for non-double inputs
X0 = X0(:);
if ~isa(X0,'double')
    error('OCP:NonDoubleInput', ...
        'OCP only accepts inputs of data type double.')
end
n = numel(X0);
numberOfVariables = n;

% Check if R is a positive definite matrix
if  nargin < 5 || isempty(R)
    R = eye(n);
else
    nR = size(R);
    if length(nR) > 2 || length(nR) < 2 || nR(1) ~= nR(2) || nR(1) ~= n || 1 ~= check_matrix_definiteness(R)
        error('OCP:OptRNotPositiveDefiniteMatrix',...
            'Option ''R'' must be an n-by-n positive definite matrix where n is equal to deminsion of ''X0''')

    end
end

% Check if R is a positive definite matrix
if  nargin < 6 || isempty(MaxIter)
    MaxIter = 50 * numberOfVariables;
elseif ~isscalar(MaxIter) || MaxIter < 1
    error('OCP:OptMaxIterNotInteger',...
        'Option ''MaxIter'' must be a positive integer value if not the default.')
end

% Check if TolX is a positive scalar
if  nargin < 7 || isempty(TolX)
    TolX = 1e-6;
elseif ~isscalar(TolX) || TolX <=0
    error('OCP:OptTolXNotPositiveScalar',...
        'Option ''TolX'' must be a positive scalar if not the default.')
end

% Check if Method is selected proper
if  nargin < 7 || isempty(Method)
    Method = 'Method1';
elseif ~strcmpi(Method,'Method1') && ~strcmpi(Method,'Method2')
    error('OCP:OptMethodInputError',...
        'Option ''Method'' must be selected in {''Method1''} | ''Method2''.')
end

% Check if ShowGraph is on
if  nargin < 8 || isempty(ShowGraph)
    ShowGraph = true;
elseif ~strcmpi(ShowGraph,'on') && ~strcmpi(ShowGraph,'off')
    error('OCP:OptShowGraphInputError',...
        'Option ''ShowGraph'' must be selected in {''on''} | ''off''.')
end

% Check if Fun is scalar && GradFcn is Vector with dimension = n && HessFcn
% is nSQUARE MATRIX with dimension = n
if ~isempty(FUN),
    F0 = checkfun(X0,FUN,varargin,'1');
else
    error('OCP:OptFUNInputError','Option ''FUN'' can''t be []');
end
if ~isempty(GradFcn),
    F1 = checkfun(X0,GradFcn,varargin,'2');
else
    error('OCP:OptGradFcnInputError','Option ''GradFcn'' can''t be []');
end
if ~isempty(HessFcn),
    F2 = checkfun(X0,HessFcn,varargin,'3');
else
    error('OCP:OptHessFcnInputError','Option ''HessFcn'' can''t be []');
end

% ---------------------------------------------------
% Start the Iteration Precess
% ------------------------1--------------------------
disp('Iteration is started.')
X = X0;
hf = 0;
EXITFLAG = 0;
if strcmpi(ShowGraph,'on'),callAllOptimPlotFcns('init'),end

In = eye(n); tmpX = X; tmpF = F0;
switch lower(Method)
    case 'method1'      % Algorithm in (9)-(10) of [1] or in (15) of [2]
        for j = 0:MaxIter
            g = (R + F2)\F1;
            if norm(g) < TolX,EXITFLAG = 1;break;end
            for i=1:j
                g = (R + F2) \ (F1 + R * g);
            end
                        
            X = X - g;
            F1 = feval(GradFcn, X, varargin{:});
            F2 = feval(HessFcn, X, varargin{:});
            
            if ShowGraph,callAllOptimPlotFcns('iter');end
        end
    case 'method2'      % Algorithm described in (22) of [2]
        for j = 0:MaxIter
            g = R * F1;
            if norm(g) < TolX,EXITFLAG = 1;break;end
            for i=1:j;
                g = R * F1 + (In - R * F2) * g;
            end

            X = X - g;
            F1 = feval(GradFcn, X, varargin{:});
            F2 = feval(HessFcn, X, varargin{:});
            
            if ShowGraph,callAllOptimPlotFcns('iter');end
        end
    otherwise
        error('OCP:OptMethodError','Unknown method.')
end
FVAL = feval(FUN, X, varargin{:});
GRAD = F1;
HESSIAN = F2;

if ShowGraph && ishandle(hf), figure(hf);end

disp('Iteration is Finished.')
% pause
% if ishandle(hf),close(hf),end

% ----------------------------
% Nested Fun for PLOT
% ============================
function callAllOptimPlotFcns(state)
    switch state
        case 'init'
            hf = figure('numbertitle','off','Name',['Iteration Process of ',Method,' with OCP'],'Visible','on');
            clf
            subplot(211),hold on,box on
            %     subplot(2,1,1,'YScale','log','YMinorTick','on'),hold on,box on
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
            error('OCP:callAllOptimPlotFcns:state','Unkown state.'); 
    end
end % end of callAllOptimPlotFcns

end % end of OCP

function matrix_type = check_matrix_definiteness(A)
    % input:
    % A - the matrix to be checked if positive definite
    % outputs:
    % matrix_type -: 1-'positive_definite',0- 'not positive_definite'

    if ~isequal(A, A')
        error('Matrix must be square.');
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
    % CHECKFUN checks for complex or NaN results from userfcn.

    f = feval(userfcn, x, varargin{:});
    % Note: we do not check for Inf as OCP handles it naturally.
    if isnan(f)
        error('OCP:checkfun:NaNFval', ...
            'User function ''%s'' returned NaN when evaluated;\n OCP cannot continue.', ...
            localChar(userfcn));  
    elseif ~isreal(f)
        error('OCP:checkfun:ComplexFval', ...
            'User function ''%s'' returned a complex value when evaluated;\n OCP cannot continue.', ...
            localChar(userfcn));  
    end
    switch flag
        case '1'
            if ~isscalar(f),
                error('OCP:checkfun:ScalarFval', ...
                    'User function ''%s'' returned a non-scalar value when evaluated;\n OCP cannot continue.', ...
                    localChar(userfcn)); 
            end

        case '2'
            m = size(f);
            if length(m) ~= 2 || length(f(:)) ~= length(x(:)) || numel(f)~= size(f,1)
                error('OCP:checkfun:VectorGradFval', ...
                    'User function ''%s'' did not return a column vector with dimension not matching x when evaluated;\n OCP cannot continue.', ...
                    localChar(userfcn)); 
            end

        case '3'
            n = length(x(:));
            m = size(f);
            if length(m) ~= 2 || m(1) ~= n || m(2) ~= n
                error('OCP:checkfun:ScalarFval', ...
                    'User function ''%s'' did not return a square matrix with dimension not matching x when evaluated;\n OCP cannot continue.', ...
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