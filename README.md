Đề tài: Mô phỏng hiệu suất mạng không dây trong môi trường có chướng ngại vật.



\# Giới thiệu



&#x20; Dự án này mô phỏng hiệu suất mạng không dây trong môi     trường có vật cản bằng MATLAB. Mô phỏng được thực hiện bằng MATLAB, nhằm phản ánh gần với điều kiện thực tế khi triển khai mạng, đặc biệt là sự khác biệt giữa:

&#x20;   LOS (Line-of-Sight)

&#x20;   NLOS (Non-Line-of-Sight)



\# Mục tiêu

&#x20; Xây dựng môi trường mô phỏng 2D có vật cản ngẫu nhiên

&#x20; Xác định trạng thái truyền: LOS và NLOS

&#x20; Áp dụng mô hình suy hao theo 3GPP TR 38.901

&#x20; Đánh giá hiệu suất mạng qua:

&#x20;   Suy hao đường truyền (Path loss) 

&#x20;   Tốc độ truyền dữ liệu (Throughput)

&#x20; So sánh độ chính xác của hai thuật toán:

    kNN

    kN-MAP



\# Môi trường mô phỏng



  Bán kính vùng: 150 m

  Người dùng phân bố ngẫu nhiên

  Cột phát sóng đặt tại (0,0)

  Vật cản dạng hình chữ nhật sinh ngẫu nhiên (PPP)



\# Cấu trúc thư mục



Mo-phong-hieu-suat-mang-khong-day-trong-moi-truong-co-chuong-ngai-vat-main/

 ├── main\_code/            # Code chính MATLAB

 └── paper\_simulations/     # Kết quả mô phỏng





\# Hướng dẫn chạy



1\. Chạy mô phỏng chính



&#x20; Bước 1: Mở MATLAB

&#x20; Bước 2: Di chuyển đến thư mục chứa dự án

&#x20; Bước 3: Chạy file main\_simulation.m



2\. Chạy khảo sát theo mật độ vật cản



&#x20; Chạy file: run\_sweep\_lambda.m

&#x20; Kết quả thu được: Đánh giá ảnh hưởng của mật độ vật cản đến độ chính xác của hai thuật toán kNN kNMAP, suy hao đường truyền và tốc độ truyền dữ liệu.



3\. Chạy khảo sát theo tần số 



&#x20; Chạy file: run\_sweep\_frequency.m

&#x20; Kết quả thu được: Đánh giá ảnh hưởng của tần số truyền dẫn đến độ chính xác của hai thuật toán kNN kNMAP, suy hao đường truyền và tốc độ truyền dữ liệu.



\# Tác giả



Sinh viên ngành Mạng máy tính \& Truyền thông dữ liệu



Nguyễn Tuấn Dương

Nguyễn Ngọc Bảo

Nguyễn Thị Thanh Trà

