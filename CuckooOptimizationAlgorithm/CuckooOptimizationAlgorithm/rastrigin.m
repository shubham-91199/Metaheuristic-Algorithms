function y = rastrigin(x)

% global FunctionCallNum


y = 10.0 * size(x,2) + sum(x .^2 - 10.0 * cos(2 * pi .* x),2);

% FunctionCallNum = floor(0.8*FunctionCallNum);
% FunctionCallNum = FunctionCallNum + 1;