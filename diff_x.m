function DX = diff_x(A, dx)

    nx = size(A, 2);
    
    DX = A;
    DX(:,1,:) = (A(:,2,:) - A(:,1,:))/dx;
    DX(:,nx,:) = (A(:,nx,:) - A(:,nx-1,:))/dx;
    DX(:,2:nx-1,:) = (A(:,3:nx,:) - A(:,1:nx-2,:))/(2*dx);
    
end

