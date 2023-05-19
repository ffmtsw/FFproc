function data = emfieldrot(data,eangle,bangle)

    % Electric field rotation angle
    if eangle ~= 0
        % Electric dipoles
        ex = data(:,1);
        ey = data(:,2);
        % Rotation
        exr = ex.*cosd(eangle) + ey.*sind(eangle);
        eyr = ex.*sind(eangle) + ey.*cosd(eangle);
    else
        exr = data(:,1);
        eyr = data(:,2);
    end

    % Magnetic field rotation angle
    if bangle ~= 0
        % Magnetic field horizontal components
        bx = data(:,3);
        by = data(:,4);
        bzr = data(:,5);
        % Rotation
        bxr = bx.*cosd(bangle) + by.*sind(bangle);
        byr = bx.*sind(bangle) + by.*cosd(bangle);
    else
        bxr = data(:,3);
        byr = data(:,4);
        bzr = data(:,5);
    end 

    data = [exr,eyr,bxr,byr,bzr];

end