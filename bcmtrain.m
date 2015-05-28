function net = bcmtrain(net,method,alg,options)
% bcmtrain - Kernel parameter optimization for Bayesian Committee Machine
%
% Synopsis:
%   net = bcmtrain(net,method,alg,options)
%   
% Arguments:
%   net: Initialized BCM structure, as output by bcminit.m
%   method: String, one of 'shared', 'individual'. If method=='shared', all
%       modules share the same hyperparameters, chosen such that the sum of all
%       module marginal likelihoods is maximized. if method=='individual', kernel
%       parameters are optimized for each modules, such that each
%       modules' marginal likelihood is maximal.
%   alg: Optimization routine to use for kernel parameters, e.g. 'scg'
%   options: Options vector for the optimization routine
%   
% Returns:
%   net: Modified BCM structure
%   
% Description:
%   For reasons of efficiency, BCM hyperparameter selection is only
%   implemented as heuristics, where the marginal likelihood of
%   individual BCM modules is considered. Two strategies are available:
%   'shared': All BCM modules share the same hyperparameters (e.g.,
%       kernel parameters or noise variance. Training is done by
%       maximizing the sum of marginal likelihoods in each module.
%   'individual': Each BCM module has its distinct set of
%       hyperparameters. Training is done by maximizing marginal
%       likelihood in each module separately.
%   For most cases, it seems that shared hyperparameters
%   (method=='shared') leads to better performance than individual
%   hyperparameters. 
%   
%   
% See also: bcm,bcminit,bcmprepare,bcmfwd
% 

% Author(s): Anton Schwaighofer, Nov 2004
% $Id: bcmtrain.m,v 1.2 2004/11/23 22:49:37 anton Exp $

error(nargchk(2, 4, nargin));
if nargin<4,
  options = zeros([1 18]);
  options(1) = 0;
  options(2) = 1e-4;
  options(3) = 1e-6;
  options(9) = 0;
  options(14) = 50;
end
if nargin<3,
  alg = 'scg';
end

% Invalidate all eventually computed prior matrices
net.invPrior = {};
net.weights = {};
  

net.method = method;
switch method
  case 'shared'
    % Use default netopt. bcmerr and bcmgrad have been adapted such that
    % they are appropriate for this type of kernel parameter
    % optimization: bcmerr computes the sum of the individual module
    % likelihoods, bcmgrad computes the sum of the gradients
    net = netopt(net, options, [], [], 'scg');
  case 'individual'
    % This is easy... just loop over all modules, train them via standard
    % evidence maximization
    for i = 1:length(net.module),
      netI = net.module(i);
      netI = netopt(netI, options, netI.tr_in, netI.tr_targets, alg);
      net.module(i) = netI;
    end
  otherwise
    error('Invalid value for parameter ''method''');
end
