Tierpsy_256_40 = {Tierpsy_256_40worms{2:end,2:end}};
Tierpsy_256_40 = reshape(Tierpsy_256_40,[numel(Tierpsy_256_40)/256,256]);
Tierpsy_256_5 = {Tierpsy_256_5worms{2:end,2:end}};
Tierpsy_256_5 = reshape(Tierpsy_256_5,[numel(Tierpsy_256_5)/256,256]);
[pc, score, ~, ~, explained] = pca(cell2mat(Tierpsy_256_5));