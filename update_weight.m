% 20211216 choosing model: update weight for each model after each step.
function [weights_mod, tderr_mod] = update_weight (rewd, weights_mod,discf, prev_modQ, lrate)
% compute curr_modQ
curr_modQ = weights_mod;
% compute temporal difference error for each model and update weights_mod
tderr_mod = rewd + discf*curr_modQ - prev_modQ; %temporal difference error
weights_mod = max(0, weights_mod + lrate*tderr_mod); 