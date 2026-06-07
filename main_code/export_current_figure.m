function export_current_figure(out_dir, base_name)
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    fig = gcf;
    set(fig, 'Color', 'w');
    drawnow;
    png_path = fullfile(out_dir, [base_name, '.png']);
    fig_path = fullfile(out_dir, [base_name, '.fig']);
    try
        exportgraphics(fig, png_path, 'Resolution', 300);
    catch
        saveas(fig, png_path);
    end
    try
        savefig(fig, fig_path);
    catch
    end
end
