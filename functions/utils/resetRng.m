function resetRng( n_seed )
%% RESETRNG Resets the random number generation seed

    if (n_seed > 0)
        rng(n_seed);
    elseif (n_seed == 0)
        rng('default');
    end
end