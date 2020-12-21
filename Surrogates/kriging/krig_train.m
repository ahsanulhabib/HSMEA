function [Krigmodel,Options] = krig_train(xtr, ytr, varargin)
%% Check for other user inputs for option
i = 1;
while i <= length(varargin)
    if strcmpi(varargin{i},'FIT_Fn') || strcmpi(varargin{i},'fitfunc') || strcmpi(varargin{i},'fit_func')
        
        FIT_Fn = varargin{i+1};
        i = i+2;
    elseif strcmpi(varargin{i},'KRG_RegressionModel') || strcmpi(varargin{i},'RegressionModel')...
            || strcmpi(varargin{i},'Regression_Model') || strcmpi(varargin{i},'Regr_Model') || strcmpi(varargin{i},'Reg_Model')...
            || strcmpi(varargin{i},'RegModel') || strcmpi(varargin{i},'RegrModel')
        
        RegressionModel = varargin{i+1};
        i = i+2;
    elseif strcmpi(varargin{i},'KRG_CorrelationModel') || strcmpi(varargin{i},'CorrelationModel')...
            || strcmpi(varargin{i},'Correlation_Model') || strcmpi(varargin{i},'Corr_Model') || strcmpi(varargin{i},'Cor_Model')...
            || strcmpi(varargin{i},'CorModel') || strcmpi(varargin{i},'CorrModel')
        
        CorrelationModel = varargin{i+1};
        i = i+2;
    elseif strcmpi(varargin{i},'KRG_Theta0') || strcmpi(varargin{i},'Theta0')...
            || strcmpi(varargin{i},'Theta_0') || strcmpi(varargin{i},'Theta_init') || strcmpi(varargin{i},'Init_theta')
        
        Theta0 = varargin{i+1};
        i = i+2;
    elseif strcmpi(varargin{i},'KRG_LowerBound') || strcmpi(varargin{i},'ThetaLB') || strcmpi(varargin{i},'Theta_LB')...
            || strcmpi(varargin{i},'Theta_LowerBound') || strcmpi(varargin{i},'Theta_Lower_Bound')
        
        Theta_LB = varargin{i+1};
        i = i+2;
    elseif strcmpi(varargin{i},'KRG_UpperBound') || strcmpi(varargin{i},'ThetaUB') || strcmpi(varargin{i},'Theta_UB')...
            || strcmpi(varargin{i},'Theta_UpperBound') || strcmpi(varargin{i},'Theta_Upper_Bound')
        
        Theta_UB = varargin{i+1};
        i = i+2;
    else
        i = i+1;
    end
end

if ~exist('FIT_Fn','var')
    FIT_Fn = [];
end
if ~exist('RegressionModel','var')
    RegressionModel = [];   
end
if ~exist('CorrelationModel','var')
    CorrelationModel = [];
end
if ~exist('Theta0','var')
    Theta0 = [];
end
if ~exist('Theta_LB','var')
    Theta_LB = [];
end
if ~exist('Theta_UB','var')
    Theta_UB = [];
end

Options = KrigSetOptions(xtr, ytr, FIT_Fn, RegressionModel, CorrelationModel, Theta0, Theta_LB, Theta_UB);
[Kriging,State] = KrigFitModel(Options);
Krigmodel = Kriging.model;
Krigmodel.perf = State.Perf;
return

function Option = KrigSetOptions(X, Y, varargin)
% Load options
% Initialize options
if nargin < 2 % Return empty if either X or Y missing
    Option.X = [];
    Option.Y = [];
    
    Option.FIT_Fn = [];
    
    Option.RegressionModel  = [];
    Option.CorrelationModel = [];
    Option.Theta0           = [];
    Option.Theta_LB         = [];
    Option.Theta_UB         = [];
    
else          % Load default values
    [npoints,nvariables] = size(X);
    
    Option.X = X;
    Option.Y = Y;
    
    Option.FIT_Fn = @krig_fit;
    
    Option.RegressionModel  = @krig_regpoly0;
    Option.CorrelationModel = @krig_corrgauss;
    Option.Theta0           = (npoints^(-1/nvariables))*ones(1, nvariables);
    Option.Theta_LB         = 1e-5*ones(1, nvariables);
    Option.Theta_UB         = 100*ones(1, nvariables);
end

for i = 1:numel(varargin)
    if ~isempty(varargin{i})
        if i == 1
            Option.FIT_Fn = varargin{i};
        elseif i == 2
            Option.RegressionModel = varargin{i};
        elseif i == 3
            Option.CorrelationModel = varargin{i};
        elseif i == 4
            Option.Theta0 = varargin{i};
        elseif i == 5
            Option.Theta_LB = varargin{i};
        elseif i == 6
            Option.Theta_UB = varargin{i};
        end
    end
end
return

function [Surrogate, State] = KrigFitModel(Option)
State = srgtsFitCreateState(Option);
if isempty(Option.Theta_LB) || isempty(Option.Theta_LB) % no optimization for theta
    [Surrogate.model, State.Perf, State.FITFuncVal] = krig_fit(...
        Option.X, ...
        Option.Y, ...
        Option.RegressionModel, ...
        Option.CorrelationModel, ...
        Option.Theta0,[],[]);
else
    [Surrogate.model, State.Perf, State.FITFuncVal] = krig_fit(...
        Option.X, ...
        Option.Y, ...
        Option.RegressionModel, ...
        Option.CorrelationModel, ...
        Option.Theta0, ...
        Option.Theta_LB, ...
        Option.Theta_UB);
end
Surrogate.model.surro_type = 'krig';
return

function srgtSTT = srgtsFitCreateState(srgtOPT)
srgtSTT.FIT_Fn    = srgtOPT.FIT_Fn;
srgtSTT.FIT_FnVal = NaN;
return
