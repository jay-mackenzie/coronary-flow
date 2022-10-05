load /home/jay/Documents/coronary-flow/Data/pres_sub3.dat
load /home/jay/Documents/coronary-flow/Data/ExtPLong.dat
load /home/jay/Documents/coronary-flow/Data/MynardSmolichLong.dat
load /home/jay/Documents/coronary-flow/Data/corPress_scott.dat
%%
close all
[~, goal] = max(corPress_scott)
[~, curr] = max(MynardSmolichLong)

l = 2^16+1;

% jay_ao = shft(MynardSmolichLong, goal, l);
% plot(jay_ao)
% hold on
% plot(corPress_scott)

shft2(MynardSmolichLong, goal, l, curr);


function out = shft2(in, goal, curr)
pts_shift = goal-curr;
l = length(in);
out = [in(1)*ones(1, pts_shift)'; in]; % shift right
out(1:pts_shift)  = out(l+1:end); % copy end to start
out = out(1:l); % trim end
end

function out = shft(in, goal, l)
out = [in; in];
[~, temp] = max(out);

while temp < goal
out = out(temp+1:end);
[~, temp] = max(out);
end
temp;
out = out([temp-goal:temp-goal+l]);
end

%%
% pts_shift = p1-p2
% 
% t = linspace(0, T, length(corPress_scott));
% dt = t(2) - t(1);
% 
% t_shift= dt*pts_shift;

%%
% jay = MynardSmolichLong;
% scott = corPress_scott;
% [~, p1] = max(pres_sub3);
% driving = shift(MynardSmolichLong, p1, p2, t_shift, t, T, dt);
% %%
% hold on
% plot(t(p1),corPress_scott(p1), 'o')
% 
% external = shift(ExtPLong*100, p1, p2, t_shift, t, T, dt);
% 
% 
% 
% figure;plot(corPress_scott);hold on;plot(pres_sub3)
% 
% figure;plot(MynardSmolichLong);hold on;plot((ExtPLong*15)+85)
    
function [out] = shift(in, p1, p2, t_shift, t, l, dt)
out = [in; in];
t_MS = linspace(t_shift-T, T+t_shift, length(out))';


out = out(t_MS>0);
t_MS = t_MS(t_MS>0);
out = out(t_MS<T);
t_MS = t_MS(t_MS<T);
t_MS = [t_MS;t_MS(end)+dt];
out = [out; out(1)];
length(out)
figure;
[~, p3] = max(out);
out = out;
end



