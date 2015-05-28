function net = bcmunpak(net,w)
% bcmunpak - Copy kernel parameters for BCM from vector
%
% Synopsis:
%   net = bcmunpak(net,w)
%   
% Arguments:
%   net: BCM structure
%   w: Vector of GP parameter
%   
% Returns:
%   net: Modified BCM structure, where the GP parameters for each module are set
%       to the values given in vector w
%   
% Description:
%   This routine is only meant to be used as a subroutine of
%   bcmtrain.m. The routine copies the gien GP parameters into each GP
%   module, as well as into the GP given in net.gpnet.
%
%   
% See also: bcmpak
% 

% Author(s): Anton Schwaighofer, Nov 2004
% $Id: bcmunpak.m,v 1.1 2004/11/18 21:21:11 anton Exp $

error(nargchk(2, 2, nargin));
error(consist(net, 'bcm'));

net.gpnet = gpunpak(net.gpnet, w);
for i = 1:length(net.module),
  net.module(i) = gpunpak(net.module(i), w);
end
