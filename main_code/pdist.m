function D = pdist(X)
% PDIST Replacement function for missing Statistics Toolbox
% Tính toán khoảng cách Euclidean giữa các cặp điểm.
% Input: X - Ma trận kich thuoc m x n (m diem, n chieu).
% Output: D - Vector hang chua khoang cach.

    % Kiem tra so chieu cua X
    [m, n] = size(X);
    
    % Neu chi co 1 diem thi tra ve rong
    if m < 2
        D = [];
        return;
    end

    % Tinh khoang cach Euclidean bang phep toan vector hoa (vectorized)
    % Cong thuc: |x-y|^2 = x^2 + y^2 - 2xy
    Sq = sum(X.^2, 2); % Tong binh phuong moi hang
    % Tinh ma tran khoang cach binh phuong
    D_sq = Sq + Sq.' - 2 * (X * X.');
    
    % Lay can bac 2 (dung gia tri nho de tranh loi so hoc)
    D_matrix = sqrt(max(0, D_sq));
    
    % Trich xuat cac phan tu duong cheo duoi (lower triangle) de tao vector output
    % Thu tu cua pdist la (2,1), (3,1), (3,2), ...
    D = D_matrix(logical(tril(ones(m), -1))).';
end