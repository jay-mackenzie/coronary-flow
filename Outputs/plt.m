% load what you want to plot
out100 = load("850trif_100.2d");
%%
for o = {out100} % put the names of what you want to plot here
    plotter(o{1})
end

function plotter(out)
    for i = unique(out(:, 1))' % loop through vessel ids
        tmp = out(out(:, 1) == i, 2:end);

        x = unique(tmp(:, 2))'; % find spatial point in those vessels
        t = unique(tmp(:, 1)); % find time
        T = max(t) - min(t) + t(2)-t(1); % find period
        t = linspace(0, T, length(t)); % shift solves to compare sims
        leg = {}; % TO DO: fill to make a legend
        figure(floor(i)+1) % the (i+1)-th figure is for the i-th vessel
        for j = x % x(end) % x: all spatial points 
            tmp2 = tmp(tmp(:, 2) == j, 3:end);
            for k = 1:3 % plot flp, q, tmp
                subplot(1, 3, k); hold on
                plot(t, tmp2(:,k), 'LineWidth',2)
            end
        end
    end
end
