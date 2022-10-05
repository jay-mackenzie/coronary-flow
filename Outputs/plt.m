% load what you want to plot
% out = load('out.2d');
out = load('scott/out.2d');
% out = load('coronary-flow3/Outputs/one_0.2d');
% out = load('V16.2d');

% out100 = load("650trif_100.2d");
% out0 = load("650trif_0.2d");
% out100 = load("scott_100.2d");
% out0 = load("scott_0.2d");
% %%
% close all

T = 1.1;
for i = 0:2
    figure(i+1); hold on
    out0 = out(find(out(:, 1) == i), :);
    j = unique(out0(:, 3))';
    for k = j
        out1 = out0(find(out0(:, 3) == k), 4:end);
% t = linspace(0, T, length(out1));
% xlim([0 128])
% t = t(t-(max(t)-T) >0);
        plot(out1(:, 2));
% xlim([max(t)-T, max(t)])
        hold on
    end
shg
end
%%

% figure;
% for o = out%{out0,out100} % put the names of what you want to plot here
%     plotter(o)
% end

function plotter(out)
    for i = unique(out(:, 1))' % loop through vessel ids
        
        tmp = out(out(:, 1) == i, 2:end);

        x = unique(tmp(:, 2))'; % find spatial point in those vessels
        t = unique(tmp(:, 1)); % find time
        T = max(t) - min(t) + t(2)-t(1); % find period
        t = linspace(0, T, length(t)); % shift solves to compare sims
        leg = {}; % TO DO: fill to make a legend
%         figure(floor(i)+1) % the (i+1)-th figure is for the i-th vessel
%         figure;
        for j = x % x(end) % x: all spatial points 
            tmp2 = tmp(tmp(:, 2) == j, 3:end);
            for k = 1:3 % plot flp, q, tmp
                subplot(1, 3, k); hold on
                plot(t, tmp2(:,k), 'LineWidth',2)
                
            end
            break
        end
        break
    end
end
