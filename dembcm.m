function dembcm()
% dembcm - Demo program for BCM approximation for large scale GP regression
%
% Synopsis:
%   dembcm;
%
% Description:
%   This routine demonstrates how the provided routines for the Bayesian
%   Committee Machine can be used for analyzing data
%   Basic steps are
%   - Generate a data set (linear combination of some random basis
%     functions) of 500 data points
%   - Split the data into modules of 100 points each
%   - Train Gaussian process models for each module
%   - Use the BCM approximation to obtain a prediction
%   
%   Also, the demo compares the prediction accuracy with
%   - A Gaussian Process model that is trained on all 500 points
%   - A Gaussian Process model trained on only 300 points
%
%   The Bayesian Committee Machine is used in two variants:
%   - In the standard form, training data are assigned modules at random
%   - Alternatively, use k-means clustering, and assign points the points
%     from each cluster to a module.
%
% See also: bcm,gp,bcminit,bcmtrain,bcmprepare
%   

% Author(s): Anton Schwaighofer, Nov 2004
% $Id: dembcm.m,v 1.1 2005/11/16 17:12:41 anton Exp $

% The demo also works for other random states, no worries ;-)
randstate = 1;

rand('state', randstate);
randn('state', randstate);

% ----------------------------------------------------------------------
fprintf('Generating training and test data...\n');
% 500 training data points, save basis points to later generate test
% data. Generate low noise data
noiselevel = 0.1;
[Xtrain, Ytrain, Xbasis, Ybasis, Ytrain0] = art_data(500, 5, 0, noiselevel);
% Generate 2000 test data from the same function. Use the "true" function
% values Ytest0 (not the ones corrupted by noise) for testing
[Xtest, Ytest, dummy1, dummy2, Ytest0] = art_data(2000, 5, 0, noiselevel, ...
                                                  5, Xbasis, Ybasis);

% Options for scg:
scgopt = foptions;
scgopt(1) = 1;
scgopt(2) = 1e-4;
scgopt(3) = 1e-4;
scgopt(14) = 15;

% ----------------------------------------------------------------------
fprintf('Training a full GP model on all training data...\n');

% Full gp model: Use the standard Netlab routines to train the thingy
% Rational quadratic kernel is much better than the squexp kernel
fullgp = gp(5, 'ratquad');
fullgp = gpinit(fullgp, Xtrain, Ytrain);
fullgp = netopt(fullgp, scgopt, Xtrain, Ytrain, 'scg');
[fullpred,fullvar] = gpfwd(fullgp, Xtest);

% ----------------------------------------------------------------------
fprintf('Training a full GP model on 300 (out of the 500) training data...\n');

% We also compare with a full GP trained on only 300 (out of 500) points
full1gp = gp(5, 'ratquad');
full1gp = gpinit(full1gp, Xtrain(1:300,:), Ytrain(1:300));
full1gp = netopt(full1gp, scgopt, Xtrain(1:300,:), Ytrain(1:300), 'scg');
[full1pred,full1var] = gpfwd(full1gp, Xtest);

% ----------------------------------------------------------------------
fprintf('Training the modules of the Bayesian Committee Machine...\n');

% Now build BCM model: start with defining a 'template' GP that is the
% basis for each BCM module
gp0 = gp(5, 'ratquad');
% Build a BCM model from the template
bcm0 = bcm(gp0);
% Give the BCM its data. The training data will be split up such that
% each  module gets 100 points
fprintf('BCM: Each module has 100 data points\n');
bcm0 = bcminit(bcm0, Xtrain, Ytrain, 100);
% Fit the BCM modules. Do this by optimizing evidence for each module
% with shared hyperparameters
bcm1 = bcmtrain(bcm0, 'shared', 'scg', scgopt);
bcm1 = bcmprepare(bcm1);
[bcm1pred, bcm1var] = bcmfwd(bcm1, Xtest);

% ----------------------------------------------------------------------
fprintf('Starting to cluster the training data...\n');

% Clustered BCM:
kmeansopt = [1 1e-5 1e-4 0 0 0 0 0 0 0 0 0 0 30];
r = randperm(size(Xtrain,1));
[centres,opt,post] = kmeans(Xtrain(r(1:5),:),Xtrain,kmeansopt);
[m,assignment] = max(post,[],2);
for i = 1:5,
  fprintf('Clustered BCM: Module %i has %i data points\n', i, nnz(assignment==i));
end
% ----------------------------------------------------------------------
fprintf('Training the modules of the clustered Bayesian Committee Machine...\n');
bcm5 = bcminit(bcm0, Xtrain, Ytrain, assignment);
bcm5a = bcmtrain(bcm5, 'shared', 'scg', scgopt);
bcm5a = bcmprepare(bcm5a);
[bcm5apred, bcm5avar] = bcmfwd(bcm5a, Xtest);


% ----------------------------------------------------------------------
fprintf('\nEvaluating all models in terms of\n');
fprintf('RMSE (root mean squared error)\n');
fprintf('logProb (negative log probability of test data under the predictive distribution\n\n');

loss_negLogProb = inline('0.5*(log(2*pi*var) + ((pred-label).^2)./var)','label','pred','var');

fprintf('Full GP model:\n');
pred = fullpred; var = fullvar;
fprintf('RMSE = %f, logProb = %f\n\n', sqrt(mean((pred-Ytest0).^2)), ...
        mean(loss_negLogProb(Ytest0, pred, var)));
fprintf('Full GP model on 300 data points:\n');
pred = full1pred; var = full1var;
fprintf('RMSE = %f, logProb = %f\n\n', sqrt(mean((pred-Ytest0).^2)), ...
        mean(loss_negLogProb(Ytest0, pred, var)));
fprintf('BCM model with shared hyperparams:\n');
pred = bcm1pred; var = bcm1var;
fprintf('RMSE = %f, logProb = %f\n\n', sqrt(mean((pred-Ytest0).^2)), ...
        mean(loss_negLogProb(Ytest0, pred, var)));
fprintf('Clustered BCM model with shared hyperparams:\n');
pred = bcm5apred; var = bcm5avar;
fprintf('RMSE = %f, logProb = %f\n\n', sqrt(mean((pred-Ytest0).^2)), ...
        mean(loss_negLogProb(Ytest0, pred, var)));

% $$$ figure(10);
% $$$ clf;
% $$$ val_errorbars(Ytest0', fullpred', sqrt(fullvar)');
% $$$ set(gcf, 'Name', 'Full GP with ratquad');
% $$$ figure(11);
% $$$ clf;
% $$$ val_errorbars(Ytest0', full1pred', sqrt(full1var)');
% $$$ set(gcf, 'Name', 'Full GP on subset of 300 points');
% $$$ figure(12);
% $$$ clf;
% $$$ val_errorbars(Ytest0', bcm1pred', sqrt(bcm1var)');
% $$$ set(gcf, 'Name', 'BCM shared');
% $$$ figure(15);
% $$$ clf;
% $$$ val_errorbars(Ytest0', bcm5apred', sqrt(bcm5avar)');
% $$$ set(gcf, 'Name', 'Clustered BCM shared');

return

function [X, Y, Xbasis, Ybasis, Ynoisefree] = art_data(npoints, ndim, classification, noise, nbasis, Xbasis, Ybasis)
% ART_DATA - Generate Volker's artificial data set
%   Set all random seeds to 0 before calling to reproduce the exact data set.
%

if nargin<6,
  Xbasis = [];
  Ybasis = [];
end
if nargin<5,
  nbasis = 5;
end
if nargin<4,
  noise = 0;
end
if nargin<3,
  classification = 0;
end
if nargin<2,
  ndim = 5;
end

[Nb, dimb] = size(Xbasis);
if (Nb>0) & (~all(size(Ybasis) == [Nb, 1])),
  error('Size of basis function matrices XBASIS YBASIS does not match');
end


if Nb==0,
  % No basis functions: generate new ones randomly
  % X-Prototypes in range [-1...+1]
  Xbasis = 2*rand(nbasis, ndim)-1;
  % Target values of the prototypes
  Ybasis = randn(nbasis, 1);
  if classification & (ndim==5) & (nbasis==5),
    % Volker's modification to generate a more balanced classification
    % data-set, works only with random seed 0
    Ybasis(2) = -Ybasis(2);
    Ybasis(3) = -Ybasis(3);
  end
end

adisa = mean(mean(sqrt(dist2(Xbasis, Xbasis))));
aids_sig = adisa/4;

X = 2*rand(npoints, ndim)-1;
ad  = sqrt(dist2(X, Xbasis));
Yex = exp(-1/(2*(aids_sig*aids_sig))*(ad.*ad));
Ynoisefree = (Yex * Ybasis) ./ (Yex * ones(nbasis, 1));
Y = Ynoisefree + noise * randn(npoints,1);
if classification,
  Y = sign(Y);
end
