function storeFigure( si_name )
%STOREFIGURE Stores a figure and exports to PNG and PDF files

    % save figure
    savefig(si_name);
    
    % also export to PNG and PDF
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6.25 7.5]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 6.25 7.5]);
    set(gcf, 'renderer', 'painters');
    print(gcf, '-dpdf', [si_name, '.pdf']);
    print(gcf, '-dpng', [si_name, '.png']);
%     print(gcf, '-depsc2', 'my-figure.eps');
end