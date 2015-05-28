function g = bcmgrad(net, x, t)
% bcmgrad - Error gradient for Bayesian Committee Machine
%
% Synopsis:
%   g = bcmgrad(net)
%   
% Arguments:
%   net: BCM structure
%   
% Returns:
%   g: Gradient of the error function (marginal likelihood) with respect
%       to the kernel parameters
%   
% Description:
%   Error function and gradient are computed on the basis of the
%   pre-initialized data in each GP module, thus no data is required as
%   input.
%   
%   
% See also: bcm,bcmtrain,bcminit,bcmerr
% 

% Author(s): Anton Schwaighofer, Nov 2004
% $Id: bcmgrad.m,v 1.1 2004/11/18 21:19:46 anton Exp $

g = 0;
for i = 1:length(net.module),
  netI = net.module(i);
  gI = gpgrad(netI, netI.tr_in, netI.tr_targets);
  g = g+gI;
end
