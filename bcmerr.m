function [e,edata,eprior] = bcmerr(net, x, t, exactCov, Xtest)
% bcmerr - Error function for Bayesian Committee Machine
%
% Synopsis:
%   [e,edata,eprior] = bcmerr(net)
%   e = bcmerr(net, [], [], Xtest) (use with care, see below)
%   
% Arguments:
%   net: BCM structure
%   Xtest: [Q net.nin] matrix of test data
%   
% Returns:
%   e: Value of the error function (marginal likelihood)
%   edata: Data contribution to e
%   eprior: Prior contribution to e
%   
% Description:
%   This function returns the sum of the error functions of each module,
%   that is, the unnormalized negative log-likelihood.  Error function is
%   computed on the basis of the pre-initialized data in each GP module,
%   thus no data is required as input. Still, to be compatible with the
%   standard Netlab error functions, bcmerr.m accepts input arguments in
%   the form bcmerr(net, x, t).
%
%   In the second calling syntax, bcmerr(net, [], [], Xtest), the exact
%   BCM evidence is returned. This is given by
%     -1/2 log det C - 1/2 t' C^{-1} t + const
%   with C given by
%      C = BD[K] + sigma^2 I + (K_c K_t^{-1} K_c' - BD[K_c K_t^{-1} K_c']
%   where BD[...] denotes a block-diagonal approximation of the argument,
%   K_c is the kernel matrix of all training points versus test points,
%   K_t is the test point kernel matrix, and t is a vector of all
%   training targets.
%   C is a matrix of size [N N] for N training data, thus the exact BCM
%   evidence can only be computed for cases where also an exact GP
%   solution can be found. Use this feature only with moderately sized
%   problems.
%   
%   
% See also: bcm,bcmtrain,bcminit
% 

% Author(s): Anton Schwaighofer, Nov 2004
% $Id: bcmerr.m,v 1.2 2004/11/23 21:43:51 anton Exp $

% Input arguments x and t are effectively ignored

error(nargchk(3, 4, nargin));
if nargin<4,
  Xtest = [];
end

if isempty(Xtest),
  % No test data given: Default case of summing up individual modules' evidence
  e = 0;
  edata = 0;
  eprior = 0;
  for i = 1:length(net.module),
    netI = net.module(i);
    [a, b, c] = gperr(netI, netI.tr_in, netI.tr_targets);
    e = e+a;
    edata = edata+b;
    eprior = eprior+c;
  end
else
  % Compute BCM error with the actual BCM covariance matrix. This matrix
  % has size [N N] for N training points, thus we can typically not hold
  % it in memory. For this, knowledge of the test data is required.

  % For the block diagonal approximation, we need to know the number of
  % data in each module:
  modSize = zeros(1, length(net.module));
  for i = 1:length(net.module),
    modSize(i) = length(net.module(i).tr_targets);
  end
  % Reconstruct the full training data
  N = sum(modSize);
  Xtrain = zeros(N, net.nin);
  ind = 1;
  for i = 1:length(net.module),
    netI = net.module(i);
    Xtrain(ind:(ind+length(netI.tr_targets)-1),:) = netI.tr_in;
    ind = ind+length(netI.tr_targets);
  end
  % Major part of the overall kernel matrix is a form of Schur complement:
  Kt = gpcovarp(net.gpnet, Xtest, Xtest);
  Kc = gpcovarp(net.gpnet, Xtrain, Xtest);
  smallEye = eps^(2/3)*speye(size(Kt));
  C = Kc*inv(Kt+smallEye)*Kc';
  % Overwrite diagonal blocks with exact covariance matrix, meaning that
  % the kernel matrix is exact for points within the same module
  startInd = 1;
  for i = 1:length(net.module),
    ind = startInd:(startInd+modSize(i)-1);
    netI = net.module(i);
    % Use gpcovar here, so that the contribution of the noise variance is
    % already taken into account
    C(ind,ind) = gpcovar(net.gpnet, Xtrain(ind,:));
    startInd = startInd+modSize(i);
  end
  % With this matrix C, we can compute evidence as usual:
  C(isnan(C)) = realmax;
  C(isinf(C)&(C<0)) = -realmax;
  C(isinf(C)&(C>0)) = realmax;
  eigC = eig(C, 'nobalance');
  % Guard against eventual tiny negative eigenvalues (eg. in the Matern
  % kernel with large values of nu)
  if any(eigC<=0),
    warning('Skipping some negative eigenvalues. Results may be inaccurate');
  end
  edata = 0.5*(sum(log(eigC(eigC>0)))+t'*inv(C)*t);
  eprior = 0;
  e = edata+eprior;
end
