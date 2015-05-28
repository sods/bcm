function net = bcminit(net, Xtrain, Ytrain, assignment)
% bcminit - Initialization for Bayesian Committee Machine (BCM)
%
% Synopsis:
%   net = bcminit(net,Xtrain,Ytrain,assignment)
%   
% Arguments:
%   net: BCM structure, as output by bcm.m
%   Xtrain: [N d] matrix of training data, N points in d dimensions
%   Ytrain: [N 1] vector of training targets
%   assignment: Scalar or [N 1] vector. Number of training data that are
%       assigned to each module. If assignment is a scalar K, each module is
%       assigned K points, eventually the last module is given fewer
%       points. If assignment is a vector of length N, module I will be
%       assigned all points J for which assignment(J)==I.
%   
% Returns:
%   net: BCM structure. net has the following newly added fields:
%       .module: Structure array, containing the Netlab-like GP
%           description for module I in net.module(I)
%       .invPrior: Cell array, with the inverse covariance of module I's 
%           training data in net.invPrior{I} (empty)
%       .weight: Cell array, with the GP weight vector of module I in
%           net.weight{I} (empty)
%   
% Description:
%   Purpose of this routine is to split up the training data into data
%   for the individual GP modules. Each module is a replica of the
%   template GP, given in bcm.gpnet 
%   
%   
% Examples:
%   Give each module an equal share of 500 training data:
%       net = bcminit(net, Xtrain, Ytrain, 500);
%   For improved performance, use kmeans clustering to get modules that
%   are spatially separated. Intialize kmeans with 5 random centres:
%       options = [1 1e-5 1e-4 0 0 0 0 0 0 0 0 0 0 30];
%       r = randperm(size(Xtrain,1));
%       centres = Xtrain(r(1:10));
%       [centres,opt,post] = kmeans(centres,Xtrain,options);
%   Extract the point/centres assignment and use directly in bcminit:
%       [m,assignment] = max(post,[],2);
%       net = bcminit(net, Xtrain, Ytrain, assignment);
%
%   
% See also: bcm,bcmprepare,bcmtrain,bcmfwd
% 

% Author(s): Anton Schwaighofer, Nov 2004
% $Id: bcminit.m,v 1.1 2004/11/18 21:19:53 anton Exp $

error(nargchk(4, 4, nargin));
error(consist(net, 'bcm', Xtrain, Ytrain));

[N, dim] = size(Xtrain);
% Assignment given as a scalar:
if prod(size(assignment))==1,
  modulesize = assignment;
  nModules = ceil(N/modulesize);
  r = rem(N, modulesize);
  % Each modules gets an equal share of the training data
  if r==0,
    modulesize = repmat(modulesize, [1 nModules]);
  else
    modulesize = [repmat(modulesize, [1 nModules-1]) r];
  end
  % Generate the assignment vector: assignment(j)==i if point j goes to
  % module i
  assignment = zeros([N 1]);
  start = 1;
  for i = 1:length(modulesize),
    ind = start:(start+modulesize(i)-1);
    assignment(ind) = i;
    start = ind(end)+1;
  end
else
  if length(assignment)~=N,
    error('Length of vector assignment match match the number of training data');
  end
  % Uniquify the whole thing, this gets us rid of any nonsense data,
  % wrong module numbers, and such
  [B, dummy, assignment] = unique(assignment);
  nModules = length(B);
end

% Initialize the GPs for each module with its data
net.module = net.gpnet;
for i = 1:nModules,
  netI = net.gpnet;
  ind = (assignment==i);
  netI = gpinit(netI, Xtrain(ind,:), Ytrain(ind,:));
  net.module(i) = netI;
end
% Initialize empty data for inverse prior matrices and weight vectors
net.invPrior = [];
net.weight = [];
