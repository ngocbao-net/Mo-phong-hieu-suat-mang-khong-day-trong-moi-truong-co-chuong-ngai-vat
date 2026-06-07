function true_labels = compute_true_labels_from_shadow(total_shadow, mesh)
    in_shadow = false(size(mesh,1), 1);
    if ismethod(total_shadow, 'regions')
        shadow_regions = regions(total_shadow);
        for k = 1:length(shadow_regions)
            v = shadow_regions(k).Vertices;
            in_shadow = in_shadow | inpolygon(mesh(:,1), mesh(:,2), v(:,1), v(:,2));
        end
    else
        v = total_shadow.Vertices;
        in_shadow = inpolygon(mesh(:,1), mesh(:,2), v(:,1), v(:,2));
    end
    true_labels = ~in_shadow;
end
