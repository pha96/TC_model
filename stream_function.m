function ST = stream_function(F, G, FG, dx, dy)
ST = FG.*0 + 1e6;
ST_new = ST;
count = 0;
omega = 0.9;
nx = 501;
ny = 501;

while 1 > 0
    count = count +1 ;
    %disp(count)
    %implement BC at the 4 corners
    ST_new(1,1) = ST(1,1) + omega/4*(2*ST(1,2) - 2*dx*F(1,1) + 2*ST(2,1) - 2*dx*G(1,1) - 4*ST(1,1)- dx^2*FG(1,1));
    ST_new(1,nx) = ST(1,nx) + omega/4*(2*ST(1,nx-1) + 2*dx*F(1,nx) + 2*ST(2,nx) - 2*dx*G(1,nx) - 4*ST(1,nx) - dx^2*FG(1,nx));
    ST_new(ny,1) = ST(ny,1) + omega/4*(2*ST(ny,2) - 2*dx*F(ny,1) + 2*ST(ny-1,1) + 2*dx*G(ny,1) - 4*ST(ny,1) - dx^2*FG(ny,1));
    ST_new(ny,nx) = ST(ny,nx) + omega/4*(2*ST(ny,nx-1) + 2*dx*F(ny,nx) + 2*ST(ny-1,nx) + 2*dx*G(ny,nx) - 4*ST(ny,nx) - dx^2*FG(ny,nx));
    %implement BC along the 4 edges
    ST_new(2:ny-1,1) = ST(2:ny-1,1) + omega/4*(2*ST(2:ny-1,2) - 2*dx*F(2:ny-1,1) + ST(1:ny-2,1) + ST(3:ny,1) -4*ST(2:ny-1,1) - dx^2*FG(2:ny-1,1));
    ST_new(2:ny-1,nx) = ST(2:ny-1,nx) + omega/4*(2*ST(2:ny-1,nx-1) + 2*dx*F(2:ny-1,nx) + ST(1:ny-2,nx) + ST(3:ny,nx) - 4*ST(2:ny-1,nx) - dx^2*FG(2:ny-1,nx));
    ST_new(1,2:nx-1) = ST(1,2:nx-1) + omega/4*(ST(1,1:nx-2) + ST(1,3:nx) + 2*ST(2,2:nx-1) - 2*dx*G(1,2:nx-1) - 4*ST(1,2:nx-1) - dx^2*FG(1,2:nx-1));
    ST_new(ny,2:nx-1) = ST(ny,2:nx-1) + omega/4*(ST(ny,1:nx-2) + ST(ny,3:nx) + 2*ST(ny-1,2:nx-1) + 2*dx*G(ny,2:nx-1) - 4*ST(ny,2:nx-1) - dx^2*FG(ny,2:nx-1));
    %inside
    ST_new(2:ny-1,2:nx-1) = ST(2:ny-1,2:nx-1) + omega/4*(ST(1:ny-2,2:nx-1) + ST(3:ny,2:nx-1) + ST(2:ny-1,1:nx-2) + ST(2:ny-1,3:nx) - 4*ST(2:ny-1,2:nx-1) - dx^2*FG(2:ny-1,2:nx-1));
    %compare h_new and h
    residual = abs(ST(:,:) - ST_new(:,:))/mean(ST_new(:));
    ST = ST_new;
    if (mod(count, 200) == 1)
        fprintf('%d  ------ %f \n',count,max(max(residual)));
    end
    if residual <= 0.001 | count > 20000
        break;
    end
end


