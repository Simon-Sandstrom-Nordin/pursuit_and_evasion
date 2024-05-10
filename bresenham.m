function kollision = bresenham(x0, y0, x1, y1, fieldmap, radie)
    dx = abs(x1 - x0);
    dy = abs(y1 - y0);
    sx = sign(x1 - x0);
    sy = sign(y1 - y0);
    err = dx - dy;

    cells = [];
    while true
        cells = [cells; x0, y0];


%         disp("........")
%         disp("x0:" + num2str(x0))
%         disp("x1:" + num2str(x1))
%         disp("y0:" + num2str(y0))
%         disp("y1:" + num2str(y1))
%         disp("x0 == x1:" + num2str(x0 == x1))
%         disp("y0 == y1:" + num2str(y0 == y1))
%         disp("any(x0==x1):" + num2str(any(x0==x1)))
%         disp("any(y0==y1):" + num2str(any(y0==y1)))
%         disp("any(x0 == x1) && any(y0 == y1):" + num2str(any(x0 == x1) && any(y0 == y1)))
%         disp("-----------")
%         
        if isempty(x0)
            break
        end
        if isempty(x1)
            break
        end
        if isempty(y0)
            break
        end
        if isempty(y1)
            break
        end

        if any(x0 == x1) && any(y0 == y1)
            kollision = true;
            break
        end

        for j = -radie:radie
            for k = -radie:radie
                if fieldmap(xy_index_check(x0+j, length(fieldmap)), xy_index_check(y0+k, length(fieldmap))) == 1 
                    kollision = false;
                    return
                end
            end
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
