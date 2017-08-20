function Output = refine(Mdl,Y,params0,varargin)

%REFINE Refine initial parameters to aid estimation of state-space models
%
% Syntax:
%
%   Output = refine(Mdl,Y,params0)
%   Output = refine(Mdl,Y,params0,name,value,...)
%   refine(...)
%
% Description:
%
%   For observation vector y(t), state vector x(t), and uncorrelated, unit-
%   variance white noise vector processes u(t) and e(t), maximum likelihood 
%   parameter estimation of the following state-space model:
%
%   State equation:       x(t) = A(t) * x(t-1) + B(t) * u(t)
%   Observation equation: y(t) = C(t) * x(t)   + D(t) * e(t)
%
%   is often sensitive to the initial parameters specified by the user.
%
%   In the event a preliminary estimation fails to converge, or converges to
%   an unsatisfactory solution, this utility function refines the original 
%   initial parameters in an attempt to improve model estimation results.
%
%   In the model above, the length of x(t), y(t), u(t), and e(t) is m, n, k, 
%   and h, respectively.
%
% Input Arguments:
%
%   Mdl - A state-space model with unknown parameters to estimate, created 
%     by the SSM constructor.
%
%   Y - Observed response data to which the model is fit. For time-invariant 
%     models in which the length of each observation vector (n) is the same, 
%     Y is a T-by-n matrix. For time-varying models in which the length of 
%     the observation vector changes, Y is a T-by-1 cell array in which 
%     each element contains a time-varying n-element vector of observations, 
%     y(t), associated with the corresponding period. The last observation 
%     is the most recent.
%
%   params0 - A vector containing the initial values of unknown 
%     parameters associated with model coefficients A, B, C, and D, and 
%     optionally the mean vector (Mean0) and covariance matrix (Cov0) of 
%     initial states x(0), estimated by maximum likelihood. For models created 
%     explicitly, parameters mapped to NaN values are found by a column-wise 
%     search of A, followed by B, then C, then D, and finally Mean0 and Cov0. 
%     For models created implicitly, the parameter function ParamMap is 
%     solely responsible for mapping the initial parameter vector into model 
%     coefficients A, B, C, and D, as well as additional information 
%     regarding initial states and types if necessary. 
%
% Optional Input Name/Value Pairs:
%
%   'Predictors'  T-by-d matrix of common predictor variables used to
%                 include a regression component in the observation equation. 
%                 Observations at time t are deflated such that
%
%                 [y(t) - z(t)*b] = C * x(t) + D * e(t)
%
%                 where z(t) is a vector of predictor variables and b is 
%                 the regression coefficient vector (see below). The default
%                 is an empty matrix (no regression component)
%
%   'Beta0'       d-by-n matrix of initial values of regression coefficients 
%                 associated with predictors (see above). If the model contains
%                 a regression component, coefficients are estimated along 
%                 with other unknown parameters in A, B, C, D, Mean0, and 
%                 Cov0; the default initial values are obtained by ordinary 
%                 least squares (OLS) by regressing Y on the explanatory 
%                 variables.
%
% Output Argument:
%
%   Output - Output structure array in which each element has the following 
%     fields:
%
%     o Description    A brief description of the refinement algorithm.
%                      Current algorithms are:
%
%                      'Quasi-Newton'
%                      'Nelder-Mead simplex'
%                      'Loose bound interior point'
%                      'Starting value perturbation'
%                      'Starting value shrinkage'
%
%     o Parameters:    Vector of initial parameter values.
%
%     o LogLikelihood: Log-likelihood associated with the initial vector.
%
% Notes:
%
%  o When called with no output argument, a tabular display of summary
%    information related to the resulting vector of initial parameter values
%    associated with each algorithm is printed to the screen. If the Output
%    structure array is requested, then no summary display is printed.
%
%  o The various refined initial parameter vectors may look similar to each 
%    other and to the original (see params0 input above). Manually select a
%    candidate vector of initial parameter values that makes economic sense
%    and is associated with a relatively large log-likelihood value.
%
%  o If a refinement attempt is unsuccessful, its corresponding log-likelihood 
%    is set to -Inf and its parameter vector is empty. Error messages are 
%    displayed on the screen.
% 
% See also SSM, DSSM, ESTIMATE, FILTER, SMOOTH, FORECAST, SIMULATE, SIMSMOOTH.

% Copyright 2015 The MathWorks, Inc.


isDiffuseMdl = isa(Mdl,'dssm');

% Parse regressors
if isDiffuseMdl
    callerName = 'dssm.refine';
else
    callerName = 'ssm.refine';
end
parseObj = inputParser;
addParameter(parseObj,'Predictors',[],@(x)validateattributes(x,{'numeric','logical'},{'2d'},callerName));
addParameter(parseObj,'Beta0',[],@(x)validateattributes(x,{'numeric','logical'},{'2d','finite'},callerName));
parse(parseObj,varargin{:});
Predictors = parseObj.Results.Predictors;
Beta0 = parseObj.Results.Beta0;
displayRefine = nargout == 0;

% Other name-value pair settings
unitFlag = 1;
sqrtFlag = 0;
tol = 1e-8;
displaySet = 'off';

wstate = warning;
cleanUp = onCleanup(@()warning(wstate));
warning('off','econ:statespace:estimate:SeparateParams');
warning('off','econ:statespace:estimate:NotDiagonalRepair');
warning('off','econ:statespace:estimate:InfLogLt');
warning('off','econ:statespace:estimate:AttemptFminsearch');
warning('off','econ:statespace:estimate:NonConvergence');
warning('off','econ:statespace:estimate:BadEstMdl');
warning('off','econ:statespace:estimate:BindingConstraints');
warning('off','econ:statespace:estimate:SingularCov');
warning('off','econ:statespace:estimate:NotPSDCov');
warning('off','econ:statespace:estimate:LargeSE');
warning('off','econ:statespace:statespace:TooManyConstant');
warning('off','econ:statespace:statespace:MismatchedTime');
warning('off','econ:statespace:statespace:NotPSD');
warning('off','econ:statespace:statespace:NotSymmetry')
warning('off','econ:statespace:statespace:ConstOneMean');
warning('off','econ:statespace:statespace:ConstOneCov');
warning('off','econ:statespace:statespace:TransposedData');
warning('off','econ:statespace:statespace:SquareRootInactive');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

% Determine whether it is suitable for univariate treatment
[A,B,C,D] = getABCD(Mdl,params0);
timeInvariant = ~iscell(A) && ~iscell(B) && ~iscell(C) && ~iscell(D);
if ~timeInvariant && ~ssm.isDiag(D)
    unitFlag = 0;
end

if displayRefine
    disp('Trying to refine starting values...')
end
Output = struct('Description',{},'LogLikelihood',{},'Parameters',{});

% Quasi-newton algorithm
Output(1).Description = 'Quasi-Newton';
try
    options = optimoptions('fminunc','Algorithm','quasi-newton','Display','off','MaxIter',5000);
    [~,Output(1).Parameters,~,Output(1).LogLikelihood] = estimate(Mdl,Y,params0,...
        'Univariate',unitFlag,'SquareRoot',sqrtFlag,'Tolerance',tol,...
        'Predictors',Predictors,'Beta0',Beta0,'Options',options,...
        'Display',displaySet);
    if displayRefine
        fprintf('logL = %4.3e, params = ',Output(1).LogLikelihood);
        internal.econ.tableprint(Output(1).Parameters(:)','numDigits',3);
    end
catch Exception
    disp(Exception.message)
    Output(1).LogLikelihood = -Inf;
end


% Nelder-Mead simplex search 
Output(2).Description = 'Nelder-Mead simplex';
try
    options = struct('NelderMeadSimplex',1,'BHHH',1);
    [~,Output(2).Parameters,~,Output(2).LogLikelihood] = estimate(Mdl,Y,params0,...
        'Univariate',unitFlag,'SquareRoot',sqrtFlag,'Tolerance',tol,...
        'Predictors',Predictors,'Beta0',Beta0,'Options',options,...
        'Display',displaySet);
    if displayRefine
        fprintf('logL = %4.3e, params = ',Output(2).LogLikelihood);
        internal.econ.tableprint(Output(2).Parameters(:)','numDigits',3);
    end
catch Exception
    disp(Exception.message)
    Output(2).LogLikelihood = -Inf;
end

% Loose bound interior point
Output(3).Description = 'Loose bound interior point';
try    
    if ~isempty(Predictors)
        Beta0Use = Beta0;
        if isempty(Beta0Use)
            fullData = all(~isnan([Y,Predictors]),2);
            PredictorsUse = Predictors(fullData,:);
            Yuse = Y(fullData,:);
            Beta0Use = PredictorsUse \ Yuse;
        end
        params0Big = [params0(:);Beta0Use(:)];
    else
        params0Big = params0;
    end    
    plusMask = params0Big > 0;
    lb = params0Big;
    lb(plusMask) = lb(plusMask) * 0.01;
    lb(~plusMask) = lb(~plusMask) * 10;
    options = optimoptions('fmincon','Algorithm','interior-point','Display','off','MaxIter',5000);    
    [~,Output(3).Parameters,~,Output(3).LogLikelihood] = estimate(Mdl,Y,params0,...
        'Univariate',unitFlag,'SquareRoot',sqrtFlag,'Tolerance',tol,...
        'Predictors',Predictors,'Beta0',Beta0,'Options',options,...
        'Display',displaySet,'lb',lb);
    if displayRefine
        fprintf('logL = %4.3e, params = ',Output(3).LogLikelihood);
        internal.econ.tableprint(Output(3).Parameters(:)','numDigits',3);
    end
catch Exception
    % Handle the case of Beta0 in params0
    if strcmp(Exception.identifier,'econ:statespace:estimate:MismatchedConstraint')
        try
            [~,Output(3).Parameters,~,Output(3).LogLikelihood] = estimate(Mdl,Y,params0,...
                'Univariate',unitFlag,'SquareRoot',sqrtFlag,'Tolerance',tol,...
                'Predictors',Predictors,'Beta0',Beta0,'Options',options,...
                'Display',displaySet,'lb',lb(1:length(params0)));
            if displayRefine
                fprintf('logL = %4.3e, params = ',Output(3).LogLikelihood);
                internal.econ.tableprint(Output(3).Parameters(:)','numDigits',3);
            end
        catch Exception
            disp(Exception.message)
            Output(3).LogLikelihood = -Inf;
        end
    else
        disp(Exception.message)
        Output(3).LogLikelihood = -Inf;
    end    
end

% Starting value perturbation
Output(4).Description = 'Starting value perturbation';
try
    cutoff = round(numel(params0)/3);
    params0New = params0;
    params0New(1:cutoff) = params0New(1:cutoff) * 1.01 - sqrt(eps);
    params0New(cutoff+1:end) = params0New(cutoff+1:end) * 0.99 + sqrt(eps);    
    options = optimoptions('fminunc','Algorithm','quasi-newton','Display','off','MaxIter',5000);
    [~,Output(4).Parameters,~,Output(4).LogLikelihood] = estimate(Mdl,Y,params0New,...
        'Univariate',unitFlag,'SquareRoot',sqrtFlag,'Tolerance',tol,...
        'Predictors',Predictors,'Beta0',Beta0,'Options',options,...
        'Display',displaySet);
    if displayRefine
        fprintf('logL = %4.3e, params = ',Output(4).LogLikelihood);
        internal.econ.tableprint(Output(4).Parameters(:)','numDigits',3);
    end
catch Exception
    disp(Exception.message)
    Output(4).LogLikelihood = -Inf;
end

% Starting value shrinkage
Output(5).Description = 'Starting value shrinkage';
try
    params0New = params0;
    params0New(1) = params0New(1) * 1.01 - sqrt(eps);
    params0New(end) = params0New(end) * 0.99 + sqrt(eps);
    params0New = params0New ./ 9.8;
    options = optimoptions('fminunc','Algorithm','quasi-newton','Display','off','MaxIter',5000);
    [~,Output(5).Parameters,~,Output(5).LogLikelihood] = estimate(Mdl,Y,params0New,...
        'Univariate',unitFlag,'SquareRoot',sqrtFlag,'Tolerance',tol,...
        'Predictors',Predictors,'Beta0',Beta0,'Options',options,...
        'Display',displaySet);
    if displayRefine
        fprintf('logL = %4.3e, params = ',Output(5).LogLikelihood);
        internal.econ.tableprint(Output(5).Parameters(:)','numDigits',3);
    end
catch Exception
    disp(Exception.message)
    Output(5).LogLikelihood = -Inf;
end



