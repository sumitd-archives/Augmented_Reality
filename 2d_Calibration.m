%% Part II: Camera Calibration using 2D calibration object
%% driver function which invokes other functions to solve the sub parts of the problem
function[homographies k_matrix r_matrices t_matrices] = part2()
    %images
    img1 = 'images2.png';
    img2 = 'images9.png';
    img3 = 'images12.png';
    img4 = 'images20.png';
    images = {img1, img2, img3, img4};
    
    %object corner co-ordinates
   % corner1 = [0, 0, 1];%top-left
   % corner2 = [0, 210, 1];%top-right
   % corner3 = [270, 0, 1];%bottom-left
   % corner4 = [270, 210, 1];%bottom-right
    corner1 = [0, 0, 1];%bottom-left
    corner2 = [270, 0, 1];%bottom-right
    corner3 = [0, 210, 1];%top-left
    corner4 = [270, 210, 1];%top-right
    corners = [corner1; corner2; corner3; corner4];
    
    %get the Homography for each image
    homographies = [];
    for i = 1 : 4
        homography = calculate_homography(images{i} , i);
        homographies(: , : , i) = homography;
    end
    
    %display the homographies for all the images
    %%display_homography(images , homographies);
    %Calculate the intrinsic parameters
    [a_matrix lambda] = get_intrinsic_params(homographies , length(images));
    %Display the extrinsic parameters
    display_extrinsic_parameters(a_matrix, lambda, homographies, images);
    %Improving Accuracy
    [homographies k_matrix r_matrices t_matrices] = improve_accuracy(corners, images, homographies);
end

%% This function calculates the homographies of an image
function[homography] = calculate_homography(image , number)
    %Calculate the corner homogeneous coordinates in the object
    %corner1 = [0, 0, 1];%top-left
    %corner2 = [0, 210, 1];%top-right
    %corner3 = [270, 0, 1];%bottom-left
    %corner4 = [270, 210, 1];%bottom-right
    %corners = [corner1; corner2; corner3; corner4];
    corner1 = [0, 0, 1];%bottom-left
    corner2 = [270, 0, 1];%bottom-right
    corner3 = [0, 210, 1];%top-left
    corner4 = [270, 210, 1];%top-right
    corners = [corner1; corner2; corner3; corner4];
    object_points = transpose(corners);
    %Calculate the corner homogeneous coordinates in the image
    %using the gipnput() function. The corners of the squares 
    %should be selected in the following order .
    %bottom-left , bottom-right , top-left , top-right
    heading = ['image : ' , num2str(number)];
    imshow(imread(image));
    title(heading);
    [x y] = ginput(4);
    npts = length(x);
    corners = [];
    for i = 1 : npts
        corner = [x(i) y(i) 1];
        corners = [corners; corner];
    end  
    image_points = transpose(corners);
    args = [object_points; image_points];
    homography = homography2d(args);
end

%% This function displays the homographies for all the images
function[] = display_homography(images, homographies)
    nImages = length(images);
    for i = 1 : nImages
        image =  ['image : ' , images{i}];
        disp(image);
        display('homography : ')
        disp(homographies(:, :, i));
    end
end
    
%% This function displays the intrinsic parameters and returns a matrix of the intrinsic params
function[a_matrix lambda] = get_intrinsic_params(homographies, nImages)
    V_matrix = [];
    for i = 1 : nImages
        v_11 = get_v_matrix(homographies(:, :, i)' , 1, 1);
        v_12 = get_v_matrix(homographies(:, :, i)' , 1, 2);
        v_22 = get_v_matrix(homographies(:, :, i)' , 2, 2);
        v_matrix = [v_12 ; v_11 - v_22];
        V_matrix = [V_matrix ; v_matrix];
    end  
    %get the SVD of V matrix
    [U, S, V] = svd(V_matrix);
    %disp('v matrix = ');
    %disp(v);
    v = V(:, end);
    B = [v(1) v(2) v(4); v(2) v(3) v(5); v(4) v(5) v(6)];
    disp('B Matrix = ');
    disp(B)
    %calculate the intrinsic parameters
    v_0 = (B(1, 2)*B(1, 3) - B(1, 1)*B(2, 3))/(B(1, 1)*B(2, 2)-B(1, 2)^2);
    lambda = B(3, 3) - [B(1, 3)^2 + v_0*(B(1, 2)*B(1, 3)-B(1, 1)*B(2, 3))]/B(1, 1);
    alpha = sqrtm(lambda / B(1, 1));
    beta = sqrtm(lambda*B(1, 1) / (B(1, 1)*B(2, 2) - B(1, 2)^2));
    gamma = -(B(1, 2)*alpha^2*beta/lambda);
    u_0 = gamma*v_0/alpha - B(1, 3)*alpha^2/lambda;
    a_matrix = [alpha gamma u_0; 0 beta v_0; 0 0 1];
    %%
    %Display the intrinsic parameters
    disp('u_0 = ');
    %disp(sprintf('%.8f' , u_0));
    disp((u_0));
    disp('v_0 = ');
    %disp(sprintf('%.8f' , v_0));
    disp((v_0));
    disp('alpha = ');
    %disp(sprintf('%.8f' , alpha));
    disp((alpha));
    disp('beta = ');
    %disp(sprintf('%.8f' , beta));
    disp((beta));
    disp('gamma = ');
    %disp(sprintf('%.8f' , gamma));
    disp((gamma));
end

%%Get The v matrix
function[v_ij]  = get_v_matrix(h, i, j)
    v_ij = [h(i, 1)*h(j, 1), h(i, 1)*h(j, 2) + h(i, 2)*h(j, 1), h(i, 2)*h(j, 2), h(i, 3)*h(j, 1) + h(i, 1)*h(j, 3), h(i, 3)*h(j, 2) + h(i, 2)*h(j, 3), h(i, 3)*h(j, 3)];
end

%% This function calculates and displays the extrinsic parameters for each image
function[] = display_extrinsic_parameters(A_matrix, lambda, homographies, images)
    A_inverse = inv(A_matrix);
    r_matrices = [];
    t_matrices = [];
    nImages = length(images);
    %calculate the extrinsic parameters for each image
    for i = 1 : nImages
        r_1 = (A_inverse*homographies(:, 1, i)).*lambda;
        r_2 = (A_inverse*homographies(:, 2, i)).*lambda;
        r_3 = cross(r_1, r_2);
        t = (A_inverse*homographies(:, 3, i)).*lambda;
        r_matrix = [r_1 r_2 r_3];
        r_matrices(:, :, i) = r_matrix;
        t_matrices(:, :, i) = t;
    end
    %print the R matrix and the T matrix 
    for i = 1 : nImages              
        image = ['image : ' , images{i}];
        disp(image);
        sprintf('\n');
        R = r_matrices(:, :, i);
        R_T_R = transpose(R)*R;
        disp('R = ');
        %disp(sprintf('%.4f' ,real(R)));
        disp((R));
        disp('T = ');
        %disp(sprintf('%.4f' ,real(t_matrices(:, :, i))));
        disp((t_matrices(:, :, i)));
        disp('R_transpose_R = ');
        %disp(sprintf('%.4f' ,real(R_T_R)));
        disp((R_T_R));
        %%
    end    
    %%
    disp(sprintf('Below are the R and R_Transpose_R after enforcing R to be a rotation matrix\n\n'));
    for i = 1 : nImages
        R = r_matrices(:, :, i);
        [u s v] = svd(R);
        R_modified = u*transpose(v);
        unit_matrix = transpose(R_modified)*R_modified;
        image = ['image : ' , images{i}];
        disp(image);
        disp('R_modified = ');
        %disp(sprintf('%.4f' ,real(R_modified)));
        disp(real(R_modified));
        disp('R_transpose_R = ');
        %disp(sprintf('%.4f' ,real(unit_matrix)));
        disp(real(unit_matrix));
    end
end

%% PART II : Improving accuracy
%% driver function which invokes other functions to solve the sub parts of the improving accuracy problem 
function[new_homographies k_matrix r_matrices t_matrices] = improve_accuracy(corners, images, homographies)
    nImages = length(images);
    %call the function to display the projected grid corners
    p_approx = project_grid_corners(corners, images, homographies, 1);
    %display the Harris Corners for the images
    all_harris_corners = display_harris_corners(images);
    %display the closest Harris Corner
    p_correct = display_closest_harris_corners(images, p_approx, all_harris_corners);
    %Compute the new homography for each image
    new_homographies = calculate_new_homographies(corners , p_correct);
    %display new homographies
    display_new_homographies(new_homographies, images);
    %calculate and display the k matrix
    [k_matrix lambda] = find_new_intrinsic_params(new_homographies, nImages);
    %calculate the extrinsic parameters
    [r_matrices t_matrices] = display_new_extrinsic_parameters(k_matrix, lambda, new_homographies, images);
    %get the corners points after projection using the new homographies
    p_new = project_grid_corners(corners, images, new_homographies , 0);
    %calculate the errors between p_new , p_correct and p_new, p_approx
    calculate_errors(p_new, p_correct, p_approx, images);
end

%% This function displays the projected grid corners on the images
function[all_images_corners] = project_grid_corners(corners, images, homographies, show)
    all_images_corners = [];
    for j = 1 : 4
        image = imread(images{j});
        homography = homographies(:, :, j);
        if show
            heading = ['Figure ', num2str(j), ': Projected grid corners'];
            figure
            imshow(image);
            disp(heading);
            title(heading);
            hold on;
        end
        image_corners = [];
        for i = 1 : 4
            corner = corners(i, :);
            image_corner = homography*transpose(corner);
            x = image_corner(1)/image_corner(3);
            y = image_corner(2)/image_corner(3);
            if show
                plot(x, y, 'r.' , 'MarkerSize', 30);
            end
            corner = [x y];
            image_corners = [image_corners; corner];
        end
        if show
            hold off;
        end
       all_images_corners(:, :, j) = image_corners;
    end
end
        
%% This function displays the Harris corners for all the images
function[all_harris_corners] = display_harris_corners(images)
    nImages = length(images);
    for i = 1 : nImages
        heading = ['Figure ', num2str(i), ': Harris corners'];
        im = imread(images{i});
        figure
        imshow(im);
        title(heading);
        hold on;
        [cim, r, c, rsubp, csubp] = harris(rgb2gray(im), 2, 500, 2, 0);
        nCorners = length(rsubp);
        for j = 1 : nCorners
            x = csubp(j);
            y = rsubp(j);
            plot(x, y, 'r+');
        end
        hold off;
        harris_corners = [csubp rsubp];
        all_harris_corners{i} = harris_corners;
    end    
end

%% This function calculates the closest harris corners
function[all_images_closest_harris_corners] = display_closest_harris_corners(images, all_images_corners, all_harris_corners)
    nImages = length(images);
    nCorners = 4;
    all_images_closest_harris_corners = [];
    for i = 1 : nImages
        heading = ['Figure ', num2str(i), ': grid points'];
        im = imread(images{i});
        figure
        imshow(im);
        title(heading);
        hold on;
        closest_harris_corners = [];
        image_corners = all_images_corners(:, :, i);
        harris_corners = all_harris_corners{i};
        for j = 1 : nCorners
            image_corner = image_corners(j, :);
            dist_matrix = dist2(image_corner, harris_corners);
            [sorted_dist_mat index] = sort(dist_matrix);
            closest_harris_corner = harris_corners(index(1), :);
            x = closest_harris_corner(1);
            y = closest_harris_corner(2);
            plot(x, y, 'r.' , 'MarkerSize', 30);
            closest_harris_corners = [closest_harris_corners; closest_harris_corner];
        end
        hold off;
        all_images_closest_harris_corners(:, :, i) = closest_harris_corners;
    end
end

%% This function calculates the new homographies
function[new_homographies] = calculate_new_homographies(object_corners , all_images_corners)
    nImages = 4;
    nCorners = 4;
    new_homographies = [];
    object_corners = transpose(object_corners);
    ones = [1; 1; 1; 1];
    for i = 1 : nImages
        image_corners = all_images_corners(:, :, i);
        %convert to homogeneous co-ordinates
        image_corners = [image_corners ones];
        image_corners = transpose(image_corners);
        args = [object_corners; image_corners];
        homography = homography2d(args);
        new_homographies(:, :, i) = homography;
    end   
end

%% This function displays new homographies
function[] = display_new_homographies(new_homographies, images)
    nImages = length(images);
    for i = 1 : nImages
        image =  ['image : ' , images{i}];
        display('homography : ')
        disp(new_homographies(:, :, i));
    end
end

%% This function displays the K matrix
function[k_matrix lambda] = find_new_intrinsic_params(homographies, nImages)
    V_matrix = [];
    for i = 1 : nImages
        v_11 = get_v_matrix(homographies(:, :, i)' , 1, 1);
        v_12 = get_v_matrix(homographies(:, :, i)' , 1, 2);
        v_22 = get_v_matrix(homographies(:, :, i)' , 2, 2);
        v_matrix = [v_12 ; v_11 - v_22];
        V_matrix = [V_matrix ; v_matrix];
    end  
    %get the SVD of V matrix
    [U, S, V] = svd(V_matrix);
    v = V(:, end);
    B = [v(1) v(2) v(4); v(2) v(3) v(5); v(4) v(5) v(6)];
    %calculate the intrinsic parameters
    v_0 = (B(1, 2)*B(1, 3) - B(1, 1)*B(2, 3))/(B(1, 1)*B(2, 2)-B(1, 2)^2);
    lambda = B(3, 3) - [B(1, 3)^2 + v_0*(B(1, 2)*B(1, 3)-B(1, 1)*B(2, 3))]/B(1, 1);
    alpha = sqrtm(lambda / B(1, 1));
    beta = sqrtm(lambda*B(1, 1) / (B(1, 1)*B(2, 2) - B(1, 2)^2));
    gamma = -(B(1, 2)*alpha^2*beta/lambda);
    u_0 = gamma*v_0/alpha - B(1, 3)*alpha^2/lambda;
    k_matrix = [alpha gamma u_0; 0 beta v_0; 0 0 1];
    %display the k_matrix
    disp('k matrix : ');
    disp((k_matrix));
end

%% This function calculates and displays the new extrinsic parameters for each image
function[r_matrices t_matrices] = display_new_extrinsic_parameters(k_matrix, lambda, homographies, images)
    k_inverse = inv(k_matrix);
    r_matrices = [];
    t_matrices = [];
    nImages = length(images);
    %calculate the extrinsic parameters for each image
    for i = 1 : nImages
        r_1 = (k_inverse*homographies(:, 1, i)).*lambda;
        r_2 = (k_inverse*homographies(:, 2, i)).*lambda;
        r_3 = cross(r_1, r_2);
        t = (k_inverse*homographies(:, 3, i)).*lambda;
        r_matrix = [r_1 r_2 r_3];
        r_matrices(:, :, i) = r_matrix;
        t_matrices(:, :, i) = t;
    end
    %print the R matrix and the T matrix 
    for i = 1 : nImages              
        image = ['image : ' , images{i}];
        disp(image);
        sprintf('\n');
        R = r_matrices(:, :, i);
        disp('R = ');
        disp((R));
        disp('T = ');
        disp((t_matrices(:, :, i)));
        %%
    end    
end

%% This function calculates the error between p_new , p_correct and p_new , p_approx
function[] = calculate_errors(p_new, p_correct, p_approx, images)
    nImages = 4;
    nCorners = 4;
    for i = 1 : nImages
        image =  ['image : ' , images{i}];
        new_corners = p_new(:, :, i);
        correct_corners = p_correct(:, :, i);
        approx_corners = p_approx(:, :, i);
        for j = 1 : nCorners
            new_corner = new_corners(j, :);
            correct_corner = correct_corners(j, :);
            approx_corner = approx_corners(j, :);
            dist_new_correct = dist2(new_corner, correct_corner);
            dist_new_approx = dist2(new_corner, approx_corner);
            output_corner = ['corner : ' , num2str(j)];
            output = [image , ' , ' , output_corner];
            disp(output);
            image =  ['image : ' , images{i}];
            disp('error squared between the new corner and correct corner = ');
            disp(dist_new_correct);
            disp('error squared between the new corner and approx corner = ');
            disp(dist_new_approx);
            %disp(sprintf('\n'));        
        end
        disp(sprintf('\n'));
    end
end