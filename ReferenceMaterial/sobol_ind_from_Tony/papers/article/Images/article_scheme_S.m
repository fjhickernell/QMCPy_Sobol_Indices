    set(0,'defaultaxesfontsize',28,'defaulttextfontsize',28) %make font larger
    set(0,'defaultLineLineWidth',3) %thick lines
    set(0,'defaultTextInterpreter','latex') %latex axis labels
    set(0,'defaultLegendInterpreter','latex') %latex axis labels
    set(0,'defaultLineMarkerSize',40) %latex axis labels

S = @(x) exp(-(x-1).^2) + 1;

I = 0.90;
err = 0.2;
Ireal = I - err + 0.03;
Smax = S(1); xmax = 1;
Smin = S(I - err); xmin = I - err;


x = I - err - 0.06:0.02:I + err + 0.06;
ymin = min(S(x));
ymax = max(S(x));

plot(x, S(x), 'k')
ax = gca;
set(ax,'TickLabelInterpreter', 'Latex');
ax.XTick = [I-err Ireal I I+err];
ax.YTick = [Smin S(Ireal) (Smin+Smax)/2 S(I) Smax];
ax.XTickLabel = {'$\hat{I}-\varepsilon_{I}$','$I$','$\hat{I}$','$\hat{I}+\varepsilon_{I}$'};
ax.YTickLabel = {'$S_{\min}$', '$S(I)$','$(S_{\min}+S_{\max})/2$', '$S(\hat{I})$', '$S_{\max}$'};
hold on
area(I - err:0.02:I + err, S(I - err:0.02:I + err))
plot([xmin xmin], [ymin Smin], 'k--')
plot([x(1) xmin], [Smin Smin], 'k--')
plot([I I], [ymin S(I)], 'k--')
plot([x(1) I], [S(I) S(I)], 'k--')
plot([xmax xmax], [ymin Smax], 'k--')
plot([x(1) xmax], [Smax Smax], 'k--')
plot([Ireal Ireal], [ymin S(Ireal)], 'k--')
plot([x(1) Ireal], [S(Ireal) S(Ireal)], 'k--')
hold off
axis([x(1) x(end) ymin ymax*1.01])