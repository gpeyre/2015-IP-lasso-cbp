function plot_homotopy(lambda_list, A, J, I, v)

% plot_homotopy - display a continuation path
%
%   plot_homotopy(lambda_list, A, J, I, v)
%
% Copyright (c) 2015 Gabriel Peyre

c = round(rescale(v)*255)+1;
CM = jet(256);

lw = 2;
hold on;
for i=1:size(A,1)
    style = '-';
    j = J(i);
    col = CM(c(i),:);
    if not(ismember(j,I))
        style = '--';
    end
    plot(lambda_list, A(i,:), style, 'LineWidth', lw, 'Color', col);
end
