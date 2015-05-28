function net = bcmprepare(net, verbosity)
% bcmprepare - Pre-compute prior matrices for Bayesian Committee Machine (BCM)
%
% Synopsis:
%   net = bcmprepare(net)
%   net = bcmprepare(net,verbosity)
%   
% Arguments:
%   net: Initialized BCM structure, as output by bcminit.m (training data must
%       already be assigned to each module)
%   verbosity: (optional) Use a value >0 to display progress information
%   
% Returns:
%   net: Modified BCM structure, where now the fields net.invPrior and
%       net.weight are computed
%   
% Description:
%   Pre-compute the matrices that are repeatedly used in BCM forward
%   propagation, that are the inverse covariance matrix of each module,
%   and the weight vector for GP predictions.
%   
%   
% See also: bcm,bcminit
% 

% Author(s): Anton Schwaighofer, Nov 2004
% $Id: bcmprepare.m,v 1.1 2004/11/18 21:20:55 anton Exp $

error(nargchk(1, 2, nargin));
error(consist(net, 'bcm'));
if nargin<2,
  verbosity=0;
end

if verbosity>0,
  fprintf('Pre-computing prior matrices for %i modules ', nbModules);
end
for i = 1:length(net.module),
  netI = net.module(i);
  % gpcovar computes the kernel matrix of the given points, and also adds
  % the measurement noise.
  Kprior = gpcovar(netI, netI.tr_in);
  net.invPrior{i} = inv(Kprior);
  % Measurement noise is restricted to a minimum value of 1e-8 in the
  % Netlab routines. Thus, the matrices should be so well conditioned
  % that we can solve the linear system by inversion, instead of mldivide
  net.weight{i} = net.invPrior{i} * netI.tr_targets;
  if verbosity==2,
    fprintf('.');
  end
end
if verbosity==2,
  fprintf('\n');
end
