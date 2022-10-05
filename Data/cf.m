%% load all the data and compare raw versions
close all
clear all

P = './';
S = dir(fullfile(P,'*.dat'));

for k = 1:numel(S)
    F = dir(fullfile(P,S(k).name));
    out{k} = load(F.name);
    out{k} = out{k}(1:2^16+1);
    l{k} = F.name;
end 

compare(out)



% [~, goal] = max(out{5})
%% normalise
figure;
c = 1;

for i = 1:6 % from S.name
    out{i} = out{i}-min(out{i})*ones(size(out{i}));
    out{i} = out{i}/max(out{i});
    if intersect(i, [2 3 5])
        out{i} = out{i}*36;
        out{i} = out{i}+65;
        [~, goal(i)] = max(out{i});
        ll{c} = l{i};
        subplot(1, 2, 1); hold on
        plot(out{i});
        c = c + 1;
    else
        out{i} = out{i}*100;
        subplot(1, 2, 2)
        hold on
        plot(out{i})
    end

end

%% shift peaks

goal = max(goal);
clc
figure(4); hold on
for i = [1:6]
    [~, curr(i)] = max(out{i});
        out{i} = shft(out{i}, goal, curr(i));
    if intersect(i,[2, 3, 5])
        subplot(1,2, 1); hold on; plot(out{i})
    else
        subplot(1, 2, 2); plot(out{i}); hold on
    end
end
legend(ll)

% compare(out)

%% save the data
names = {"JE", "JA", "NA", "NE", "SA", "SE"};

for i = 1:6
    p = out{i};
    a = sprintf("save -ascii %s.dat p", names{i});
    eval(a);
end


%% functions

function out = shft(in, goal, curr)
    
    pts_shift = goal-curr;
    l = 2^16+1;
    in = in(1:l);
    out = [in(1)*ones(1, pts_shift)'; in]; % shift right
    out(1:pts_shift)  = out(l+1:end); % copy end to start
    out = out(1:l); % trim end
end

function compare(out)
figure
subplot(1, 2, 1); hold on
for i = [2, 3, 5]
    plot(out{i})
end

subplot(1, 2, 2); hold on
for i = [1,4,6]
    plot(out{i})
end
pairs = [2, 1; 3 4; 5 6];

figure
for i = 1:3
    subplot(1, 3, i); hold on
    for j = pairs(i, :)
            plot(out{j})
    end
end
end
