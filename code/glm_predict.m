function [X_new, y_pred, slope, intercept] = glm_predict(X, y, X_new, distribution)
% DESCRIPTION:
%   General linear regression model prediction
%
% USAGE: 
%   X = rand(100, 1); y = rand(100, 1);
%   X_new = rand(100, 1);
%   distribution = 'normal';
%   [X_new, y_pred, slope, intercept] = glm_predict(X, y, X_new, distribution)
%
% Inputs:   X, y: features and target
%
%           X_new: features to predict
%
%           distribution: distribution of the response variable
%
% Outputs:  X_new: features to predict
%
%           y_pred: predicted target 
%
%           slope, intercept: result parameters of GLM model
    if size(X, 1) < size(X, 2)
        X = X';
    end
    if size(X_new, 1) < size(X_new, 2)
        X_new = X_new';
    end

    if ~exist('distribution', 'var')
        distribution = 'normal';
    end

    mdl = fitglm(X, y, 'Distribution', distribution, 'Link', 'identity');
    
    
    y_pred = predict(mdl, X_new);
    
    
    coefficients = mdl.Coefficients.Estimate;
    intercept = coefficients(1);
    slope = coefficients(2); 

end
