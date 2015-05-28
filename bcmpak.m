function w = bcmpak(net)
% bcmpak - Combine kernel parameters of BCM into vector
%
% Synopsis:
%   w = bcmpak(net)
%   
% Arguments:
%   net: BCM structure
%   
% Returns:
%   w: Vector of GP parameters taken from field net.gpnet
%   
% Description:
%   This routine is only meant to be used as a subroutine of
%   bcmtrain.m. The routine returns the GP parameters taken from
%   net.gpnet, which are assumed to be equal to all of the BCM's module
%   parameters.
%   
%   
% See also: bcm,bcmtrain,bcmunpak
% 

% Author(s): Anton Schwaighofer, Nov 2004
% $Id: bcmpak.m,v 1.1 2004/11/18 21:20:47 anton Exp $

error(nargchk(1, 1, nargin));
error(consist(net, 'bcm'));
w = gppak(net.gpnet);
