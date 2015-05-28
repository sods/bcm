function net = bcm(gpnet)
% bcm - Bayesian Committee Machine
%
% Synopsis:
%   net = bcm(gpnet)
%   
% Arguments:
%   gpnet: A Gaussian process template for BCM modules, as output by Netlab's
%       function gp.m. Each module of the BCM will inherit its initial
%       parameters from gpnet.
%   
% Returns:
%   net: Structure describing the BCM
%   
% Description:
%   The Bayesian Committee Machine (BCM) is an approximation method for
%   large-scale Gaussian process regression. The training data is split
%   into a number of blocks, for which individual Gaussian process
%   predictors ("modules") are trained. The prediction of a BCM is a
%   weighted combination of the predictions of individual modules on the
%   test data. Also, test data is processed in blocks, which leads to
%   improved performance.
%   The code here is a wrapper routine for Gaussian process routines
%   provided by the Netlab toolbox. Netlab is thus required for this code
%   to run.
%
% Examples:
%   Building a BCM for 7-dimensional input, where each module is a GP
%   with squared-exponential kernel:
%       gpnet = gp(7, 'sqexp');
%       net = bcm(gpnet);
%   Equip the BCM with its training data, split up into modules of size
%   500: 
%       net = bcminit(net, Xtrain, Ytrain, 500);
%   Fit each module's hyperparameters, and pre-compute a few matrices:
%       net = bcmtrain(net, 'individual');
%       net = bcmprepare(net);
%   For increased performance: cluster the training data beforehand (10
%   clusters in the example below) then assign clusters to modules:
%       options = [1 1e-5 1e-4 0 0 0 0 0 0 0 0 0 0 30];
%       r = randperm(size(Xtrain,1));
%       [centres,opt,post] = kmeans(Xtrain(r(1:10)),Xtrain,options);
%       [m,assignment] = max(post,[],2);
%       net = bcminit(net, Xtrain, Ytrain, assignment);
%       net = bcmprepare(net);
%   Now can do prediction:
%       [pred, errorBar] = bcmfwd(net, Xtest, 400);
%   
% See also: bcminit,bcmprepare,bcmtrain,bcmfwd,bcmerr,bcmgrad,bcmpak,bcmunpak
% 

% Author(s): Anton Schwaighofer, Nov 2004
% $Id: bcm.m,v 1.1 2004/11/18 21:18:24 anton Exp $

error(nargchk(1, 1, nargin));

net = struct('type', 'bcm', 'gpnet', gpnet);
net.nin = gpnet.nin;
net.nout = 1;
net.module = [];
net.invPrior = {};
net.weights = {};
