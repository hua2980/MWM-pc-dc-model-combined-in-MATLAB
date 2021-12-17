% 20211216 choosing model: choose model and display Q before each step of movement.

function [which_mod, mod_prob] = model_selection (weights_mod, beta)
% softmax: computing probability of choosing each model
mod_prob = exp(weights_mod.^beta);
mod_prob = max(0, mod_prob)./sum(mod_prob); 

% Make decision (store index into the variable 'which_model') with the mod_prob
which_mod = randsample([1,-1], 1, true, mod_prob); % 1:pc-based model and -1 :dc-based









