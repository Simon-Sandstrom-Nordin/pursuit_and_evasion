function cells = bresenham_cells(x0, y0, x1, y1, fieldmap)

    dx = abs(x1 - x0);
    dy = abs(y1 - y0);
    sx = sign(x1 - x0);
    sy = sign(y1 - y0);
    err = dx - dy;

    cells = [];
    counter = 0;
    while true

        for j = -1:1
            for k = -1:1
                if fieldmap(xy_index_check(x0+j, length(fieldmap)), xy_index_check(y0+k, length(fieldmap))) == 1 
                    kollision = false;
                    return
                end
            end
        end

        if fieldmap(x0, y0) == 1
            collision = true;
            break
        end

        if mod(counter,2) == 0
            cells = [cells; x0, y0];
        end
        counter = counter + 1;

        if any(x0 == x1) && any(y0 == y1)
            collision = true;
            cells = [cells;x1,y1];
            break
        end

        e2 = 2 * err;
        if e2 > -dy
            err = err - dy;
            x0 = x0 + sx;
        end
        if e2 < dx
            err = err + dx;
            y0 = y0 + sy;
        end
    end
end
