function prob = probability_OK_cpp(points, lambda, Lmax, N_grid)
% probability_OK_cpp - Tính toán xác suất (phiên bản an toàn)

    % 1. Kiểm tra đầu vào rỗng
    if isempty(points) || size(points, 1) == 0 || Lmax <= 0
        prob = 1.0;
        return;
    end
    
    N_grid = round(N_grid);
    if N_grid < 1, N_grid = 1; end
    
    N1 = ceil(sqrt(N_grid));
    N2 = ceil(N_grid / N1);
    
    n = size(points, 1);
    total_area = 0.0;
    dl = Lmax / N1;
    dtheta = pi / N2; 
    
    xs = points(:, 1);
    ys = points(:, 2);
    
    % 2. Vòng lặp chính
    for n1 = 1:N1
        l = n1 * dl;
        sum_aux = 0.0;
        
        for n2 = 1:N2
            theta = n2 * dtheta;
            
            l_cos_theta = l * 0.5 * cos(theta);
            l_sin_theta = l * 0.5 * sin(theta);
            
            % Tính diện tích để lọc bỏ các đa giác suy biến (diện tích ~ 0)
            areas = abs((xs .* l_sin_theta) - (ys .* l_cos_theta)) * 2;
            valid_idx = find(areas > 1e-10);
            
            p_area = 0.0;
            
            if ~isempty(valid_idx)
                rects = [];
                % Tạo đa giác TỪNG CÁI một để tránh lỗi polyshape
                for k = 1:length(valid_idx)
                    idx = valid_idx(k);
                    x_curr = xs(idx);
                    y_curr = ys(idx);
                    
                    % 4 đỉnh của đa giác
                    X_pol = [-l_cos_theta; x_curr-l_cos_theta; x_curr+l_cos_theta; l_cos_theta];
                    Y_pol = [-l_sin_theta; y_curr-l_sin_theta; y_curr+l_sin_theta; l_sin_theta];
                    
                    % Tạo polyshape cho 1 đa giác
                    pol = polyshape(X_pol, Y_pol, 'Simplify', false);
                    
                    % Gộp vào chuỗi đa giác
                    if isempty(rects)
                        rects = pol;
                    else
                        rects = [rects; pol]; % Nối thêm vào danh sách
                    end
                end
                
                % Tính union của tất cả các đa giác hợp lệ
                union_shape = union(rects);
                p_area = area(union_shape);
            end
            
            sum_aux = sum_aux + p_area * dtheta;
        end
        
        total_area = total_area + sum_aux * dl;
    end
    
    prob = exp(-lambda * total_area / (pi * Lmax));
end