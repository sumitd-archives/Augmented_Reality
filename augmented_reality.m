%% Part III: Augmented Reality
%% This function calls the augment_image() and augment_object() methods
function[] = part3(homographies, k_matrix, r_matrices, t_matrices)
    base_img1 = 'images2.png';
    base_img2 = 'images9.png';
    base_img3 = 'images12.png';
    base_img4 = 'images20.png';
    images = {base_img1, base_img2, base_img3, base_img4};
    %call the augment_image() method
    augment_image(homographies, images);
    %call the augment_object() method
    augment_object(k_matrix, r_matrices, t_matrices, images);
end

%% This function performs the image augmentation job
function[] = augment_image(homographies, images)

    
    clip_art = imread('2.jpg');
    [nX nY nZ] = size(clip_art);
    %%find the scale factor
    pix_size_x = 210/nX;
    pix_size_y = 270/nY;
    scale_fac = pix_size_x;
    if pix_size_y < pix_size_x
        scale_fac = pix_size_y;
    end
    world_x = [];
    world_y = [];
    %%
    %get world co-ordinates of the image pixels 
    for i = 1 : nX
        for j = 1 : nY
            x = j*scale_fac;
            y = (nX - i)*scale_fac;
            world_x(i, j) = x;
            world_y(i, j) = y;
        end
    end 
   %%
   %find the corresponding image points and copy the clip art pixel at that point.
   for k = 1 : length(images)
       base_img_str = images{k};
       base_img = imread(base_img_str);
        for i = 1 : nX
            for j = 1 : nY
                object_x = world_x(i, j);
                object_y = world_y(i, j);
                P = [object_x object_y 1]';
                %get image points after multiplying by homography
                image_point = homographies(:, :, k)*P;
                image_x = image_point(1) / image_point(3);
                image_y = image_point(2) / image_point(3);
                col1 = clip_art(i, j, 1);
                col2 = clip_art(i, j, 2);
                col3 = clip_art(i, j, 3);
                %Dont include the white pixels 
                if col1 > 200 && col2 > 200 && col3 > 200
                    continue;
                end           
                base_img(round(image_y) , round(image_x) , 1) = clip_art(i, j, 1);
                base_img(round(image_y) , round(image_x) , 2) = clip_art(i, j, 2);
                base_img(round(image_y) , round(image_x) , 3) = clip_art(i, j, 3);
            end
        end
        figure;
        imshow(base_img);
    end 
end

%% This function performs the object augmentation job
function[] = augment_object(k_matrix, r_matrices, t_matrices , images)
    %%
    %Compute the cube word co-ordinates
    corner_1 = [0 0 0 1];
    corner_2 = [90 0 0 1];
    corner_3 = [90 90 0 1];
    corner_4 = [0 90 0 1];
    corner_5 = [0 0 -90 1];
    corner_6 = [90 0 -90 1];
    corner_7 = [90 90 -90 1];
    corner_8 = [0 90 -90 1];
    corners = [corner_1; corner_2; corner_3; corner_4; corner_5; corner_6; corner_7; corner_8];
    nCorners = 8;

    for i = 1 : 4
        image_str = images{i};
        image = imread(image_str);
        figure
        imshow(image);
        hold on;
        disp(size(image));
        extrinsic_matrix = [r_matrices(:, :, i) t_matrices(:, :, i)];
        
        m_matrix = k_matrix*extrinsic_matrix;
        projected_corners = [];
        %%
        %project the cube points on the image
        for j = 1 : nCorners
            projected_corner = m_matrix*corners(j, :)';
            x = projected_corner(1) / projected_corner(3);
            y = projected_corner(2) / projected_corner(3);
            x = round(x);
            y = round(y);
            %plot(round(x), round(y), 'r.' , 'MarkerSize', 30);
            projected_corners = [projected_corners; [round(x) round(y)]];
        end
        %%
        %join the corresponding projected points with lines
        plot(projected_corners(1, :), projected_corners(2, :), 'b');
        plot(projected_corners(1, :), projected_corners(4, :), 'b');
        plot(projected_corners(1, :), projected_corners(5, :), 'b');
        plot(projected_corners(2, :), projected_corners(3, :), 'b');
        plot(projected_corners(2, :), projected_corners(6, :), 'b');
        plot(projected_corners(5, :), projected_corners(6, :), 'b');
        plot(projected_corners(5, :), projected_corners(8, :), 'b');
        plot(projected_corners(6, :), projected_corners(7, :), 'b');
        hold off;
    end
end    