function idx_selection = RouletteWheelSelection(probabilities)

r = rand*sum(probabilities);

c = cumsum(probabilities);

idx_selection = find(r <= c, 1, 'first');