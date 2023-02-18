function DY = diff_y(A, dy)

    ny = size(A, 1);
    
    DY = A;
    DY(1,:,:) = (A(2,:,:) - A(1,:,:))/dy;
    DY(ny,:,:) = (A(ny,:,:) - A(ny-1,:,:))/dy;
    DY(2:ny-1,:,:) = (A(3:ny,:,:) - A(1:ny-2,:,:))/(2*dy);

end

