function Idx = knnsearch(X, Y, varargin)
% KNNSEARCH Custom implementation for K-Nearest Neighbors search
%   Finds the K nearest neighbors in X for each point in Y.
%
%   Inputs:
%       X   : Reference data matrix (N x D)
%       Y   : Query data matrix (M x D) or a single point vector
%       'K' : Number of neighbors to find (default is 1)
%
%   Output:
%       Idx : Indices of the nearest neighbors in X (M x K)

    % Xử lý tham số đầu vào
    p = inputParser;
    defaultK = 1;
    addParameter(p, 'K', defaultK, @isnumeric);
    parse(p, varargin{:});
    K = p.Results.K;

    % Đảm bảo Y là ma trận (dạng hàng)
    if isvector(Y)
        Y = Y(:)'; 
    end

    [N, D] = size(X); % N điểm tham chiếu
    M = size(Y, 1);   % M điểm truy vấn

    % Khởi tạo kết quả
    Idx = zeros(M, K);

    % Tính khoảng cách và tìm K điểm gần nhất
    % Sử dụng vectorization để tối ưu tốc độ
    for i = 1:M
        % Tính khoảng cách Euclidean bình phương từ điểm Y(i) đến tất cả điểm X
        % Công thức: |x - y|^2 = sum(x^2) + sum(y^2) - 2*x*y'
        % X là NxD, Y(i) là 1xD
        sq_dists = sum((X - Y(i, :)).^2, 2);
        
        % Sắp xếp để tìm chỉ số các điểm gần nhất
        [~, sorted_indices] = sort(sq_dists);
        
        % Lấy K chỉ số đầu tiên
        Idx(i, :) = sorted_indices(1:K)';
    end
end