function [Ypred, Yvar] = bcmfwd(net,Xtest,querySize,verbosity)
% bcmfwd - Forward propagation in Bayesian Committee Machine
%
% Synopsis:
%   Ypred = bcmfwd(net,Xtest)
%   [Ypred,Yvar] = bcmfwd(net,Xtest,querySize,verbosity)
%   
% Arguments:
%   net: BCM structure
%   Xtest: [Q d] matrix of test data, Q points in d dimensions
%   querySize: Size of query set (prediction is based on blocks of test
%       points of size querySize). Default value, if omitted: 500.
%   verbosity: (optional) Use a value >0 to display progress information
%   
% Returns:
%   Ypred: [Q 1] vector of predictions (predictive mean) for each test point
%   Yvar: [Q 1] vector of predictive variances (error bars) for each test
%       point 
%   
% Description:
%   Forward propagation in Bayesian Committee Machine. The test data is
%   split up into blocks of size querySize. For each block, all GP
%   modules in the BCM make their prediction, the prediction is then
%   weighted by the inverse predictive covariance, summed and normalized
%   to give the BCM output.
%   Typically, the performance of the BCM increases as querySize
%   increases.
%   Instead of passing querySize as a parameter, it can also be set in a
%   field 'querySize' of the BCM structure.
%   
% Examples:
%   Build a BCM with modules that contain 500 training points each:
%     gp1 = gp(5, 'sqexp');
%     bcm1 = bcm(gp1);
%     bcm1.querySize = 500;
%     bcm1 = bcminit(bcm1, Xtrain, Ytrain);
%   Train the BCM, by maximizing the training data marginal likelihood
%   for each module individually:
%     bcm1 = bcmtrain(bcm1,'individual','scg');
%   Compare the predictions of the BCM with different query set size:
%     pred1 = bcmfwd(bcm1, Xtest, 10);
%     pred2 = bcmfwd(bcm1, Xtest, 800);
%   
%   
% See also: bcm,bcminit,bcmtrain,bcmprepare
% 

% Author(s): Anton Schwaighofer, Nov 2004
% $Id: bcmfwd.m,v 1.2 2004/11/23 23:23:58 anton Exp $

error(nargchk(2, 4, nargin));
error(consist(net, 'bcm', Xtest));
if nargin<3,
  querySize = [];
end
if isempty(querySize),
  if isfield(net, 'querySize'),
    querySize = net.querySize;
  else
    querySize = 500;
  end
end
if nargin<4,
  verbosity=0;
end

if isempty(net.invPrior) | isempty(net.weights),
  net = bcmprepare(net);
end
P = size(Xtest, 1);
% Number of query sets of maximum size net.querySize
nQueries = ceil(P/querySize);
nModules = length(net.module);
Ypred = zeros([P 1]);
Yvar = zeros([P 1]);

if verbosity>0,
  fprintf('\nStarting forward propagation (%i query sets).\n', nQueries);
end
if verbosity==1,
  fprintf('Query set ');
end
t1 = cputime;
for j = 1:nQueries,
  if verbosity==1,
    fprintf('%i ', j);
  end
  if verbosity==2,
    fprintf('Query set %i: ', j);
  end
  ind1 = (1+(j-1)*querySize):min(P, j*querySize);
  Xtest1 = Xtest(ind1, :);
  % A small regularization matrix for inversions
  smallEye = eps^(2/3)*speye(length(ind1));
  % Prediction for the current query set
  Ypred1 = zeros([length(ind1) 1]);
  % Overall covariance matrix for current query set
  Ycov1 = 0;
  % The original BCM where all modules share the same hyperparameters:
% $$$   K11 = gpcovarp(net.module(1), Xtest1, Xtest1);
% $$$   Ycov1 = -(nModules-1)*inv(K11+smallEye);
  startInd = 1;
  for i = 1:length(net.module),
    netI = net.module(i);
    K11 = gpcovarp(netI, Xtest1, Xtest1);
    K12 = gpcovarp(netI, Xtest1, netI.tr_in);
    % Prediction of current module
    Ypred2 = K12*net.weight{i};
    % Covariance of current module
    Ycov2 = K11-K12*net.invPrior{i}*K12';
    invYcov2 = inv(Ycov2+smallEye);
    % Add weighted prediction of the current module
    Ypred1 = Ypred1+invYcov2*Ypred2;
    % Update overall covariance matrix
    Ycov1 = Ycov1+invYcov2;
    % Instead of the above (M-1)*inv(K11) term: this has one contribution for
    % the prior covariance for each module but the last/first one. The
    % last module is usually the smallest, drop this one
    if i~=length(net.module),
      Ycov1 = Ycov1 - inv(K11+smallEye);
    end
    if verbosity==2,
      fprintf('.');
    end
  end
  % Ycov1 is the *inverse* covariance of the overall prediction
  Ycov1 = inv(Ycov1+smallEye);
  % Rescale the sum of the modules' predictions and write into result
  Ypred(ind1) = Ycov1*Ypred1;
  Yvar(ind1) = diag(Ycov1);
  if verbosity>0,
    fprintf('\n');
  end
end
