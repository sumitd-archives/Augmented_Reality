%% Part 1: Camera Calibration using 3D calibration object
%% driver function which invokes other functions to solve the sub parts of the problem
function[] = part1()
    object_points = [2 2 2; -2 2 2; -2 2 -2; 2 2 -2; 2 -2 2; -2 -2 2; -2 -2 -2; 2 -2 -2];
    image_points = [422 323; 178 323; 118 483; 482 483; 438 73; 162 73; 78 117; 522 117];
    %Plot the image points
    draw_image_points(image_points);
    %Generate the p_matrix
    p_matrix = Generate_P_Matrix(object_points , image_points);     
    %get the M matrix from svd of p_matrix
    m_matrix = perform_svd(p_matrix);
    %calculate the translation matrix vector
    calculate_translation_matrix(m_matrix);
    %calculate the matrix M'
    m_dash = calculate_m_dash(m_matrix);
    %calculate the rotation vector Rx and matrix N
    n_matrix = calculate_R_x(m_dash);
    %calculate theta_z in degrees
    k_matrix = calculate_theta_z(n_matrix);
    %calculate the camera intrinsic params
    calculate_intrinsic_params(k_matrix)
end

%% This function plots the image points
function[] = draw_image_points(image_points)
    nPoints = length(image_points);
    x = image_points(:, 1);
    y = image_points(:, 2);
    scatter(x, y , 'filled');
end

%% This function returns the two rows of the P matrix
function [rows] = GetRows(object_point , image_point)
    %convert object_point to homogeneous coordinates
    object_point = [object_point 1];
    u = image_point(1);
    v = image_point(2);
    zeroes = [0 0 0 0];
    row1 = [object_point zeroes object_point.*(-u)];
    row2 = [zeroes object_point object_point.*(-v)];
    rows = [row1; row2];
end

%% This function generates the 16 rows and 12 columns of the P Matrix
function [p_matrix] = Generate_P_Matrix(object_points , image_points)
    p_matrix = [];
    [rows_n col_n] = size(object_points);
    for i = 1 : rows_n
        object_point = object_points(i , :);
        image_point = image_points(i , :);
        rows = GetRows(object_point , image_point);
        p_matrix = [p_matrix ; rows];
    end
    %%
    %display the p_matrix
    display(p_matrix);
end

%% This function will calculate the M Matrix by performing svd on P Matrix
function [m_matrix] = perform_svd(p_matrix)
    [u,s,v] = svd(p_matrix);
    last_column_v = v(: , end);
    transpose_m_matrix = reshape(last_column_v, [], 3);
    m_matrix = transpose(transpose_m_matrix);
    %%
    %display the projection matrix M
    m_matrix
end

%% This function will calculate euclidean coordinates of camera center in the frame of reference of cube
function [] = calculate_translation_matrix(m_matrix)
    [u,s,v] = svd(m_matrix);
    camera_center_homogeneous = v(: , end);
    %convert to Euclidean co-ordinates
    x = camera_center_homogeneous(1)/camera_center_homogeneous(4);
    y = camera_center_homogeneous(2)/camera_center_homogeneous(4);
    z = camera_center_homogeneous(3)/camera_center_homogeneous(4);
    camera_center_euclidean = [x y z];
    %%
    %print the euclidean coordinates of camera center
    camera_center_euclidean
end

%% This function calculates the matrix M' composed of first 3 columns of M matrix
function [m_dash] = calculate_m_dash(m_matrix)
    m_dash = m_matrix(:,1:3);
    scale_factor = m_dash(3,3);
    m_dash = m_dash./scale_factor;
    %%
    %print matrix M'
    m_dash
end

%% This function calculates the Rx , Theta_x and N
function [n_matrix] = calculate_R_x(m_dash)
    m33 = m_dash(3, 3);
    m32 = m_dash(3, 2);
    cos_theta_x = m33 / sqrtm(m33^2 + m32^2);
    sin_theta_x = -m32 / sqrtm(m33^2 + m32^2);
    Rx_matrix = [1 0 0; 0 cos_theta_x -sin_theta_x; 0 sin_theta_x cos_theta_x];
    n_matrix = m_dash*Rx_matrix;
    %%
    %print Rx , Theta_x and matrix N
    Rx_matrix     
    theta_x = atand(sin_theta_x / cos_theta_x)
    n_matrix
end

%% This function calculates the angle of rotation about z axis(theta_z) and returns the K matrix
function [k_matrix] = calculate_theta_z(n_matrix)
    n22 = n_matrix(2, 2);
    n21 = n_matrix(2, 1);
    cos_theta_z = n22 / sqrtm(n21^2 + n22^2);
    sin_theta_z = -n21 / sqrtm(n21^2 + n22^2);
    Rz_matrix = [cos_theta_z -sin_theta_z 0; sin_theta_z cos_theta_z 0; 0 0 1];
    %calculate the K matrix
    k_matrix = n_matrix*Rz_matrix;
    %%print theta_z
    theta_z = atand(sin_theta_z / cos_theta_z)
end    

%% This function uses the k matrix to calculate the focal lengths and center of the camera
function [] = calculate_intrinsic_params(k_matrix)
    %Rescale the k matrix
    scale_factor = k_matrix(3, 3);
    k_matrix = k_matrix./scale_factor;
    %%
    %print K matrix
    k_matrix
    %%
    %Calculate the focal lengths 
    output_x = ['focal_length_x = k_matrix(1,1) = ' , num2str(k_matrix(1, 1))];
    output_y = ['focal_length_y = k_matrix(2,2) = ' , num2str(k_matrix(2, 2))];
    focal_length_x = k_matrix(1, 1);
    %display the focal length along the x any y axes in pixels
    disp(output_x);
    disp(output_y);
    %%
    %Calculate the center coordinates in pixels
    u = k_matrix(1, 3);
    v = k_matrix(2, 3);
    center = [num2str(u) , ', ' , num2str(v)];
    output = ['center = (u, v) = ' , '(' , center , ')'];
    %display center
    disp(output);
end
    