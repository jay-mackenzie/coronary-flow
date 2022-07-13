out100 = load("test2_100.2d");
out0 = load("test2_0.2d");
%%
for o = {out0, out100}
    out = o{1};
    for i = unique(out(:, 1))'
        tmp = out(out(:, 1) == i, 2:end);

        x = unique(tmp(:, 2))';
        t = unique(tmp(:, 1));
        T = max(t) - min(t) + t(2)-t(1);
        t = linspace(0, T, length(t));
        leg = {};
        c = 1;
        figure(floor(i)+1)
        for j = x(end) % space
            tmp2 = tmp(tmp(:, 2) == j, 3:end);
            for k = 1:3
                subplot(1, 3, k); hold on
                plot(t, tmp2(:,k), 'LineWidth',2)
            end
        end
    end
end
