function P_CORRE_BASIC(CORRE, contactIdx, f, ax)

hold off;
set(0, 'currentfigure', f);
set(f, 'currentaxes', ax);
hold on;

nele = size(CORRE,1);
CORRE = H_SMOOTHD1(CORRE');
CORRE = H_SMOOTHD1(CORRE');
imagesc(CORRE)

% tej = flipud(colormap('jet'));
% colormap(tej);

c1 = colorbar;
xlabel('Upper <- (Depth) -> Lower');
ylabel('Lower <- (Depth) -> Upper');

for ii = 1 : length(contactIdx)
    i = contactIdx(ii);
    if iseven(i)
        labels{ii} = num2str(i);
    else
        labels{ii} = [];
    end
end
% 
set(gca,'ydir', 'rev','ylim', [1 size(CORRE,1)], 'xlim', [1 size(CORRE,1)],...
    'ytick', linspace(1, size(CORRE,1), nele), 'yticklabel', labels, ...
    'xtick', linspace(1, size(CORRE,1), nele), 'xticklabel', labels)
xtickangle(90)

ylabel(c1, 'R')

end