function  ind = Material_Locate_Handle(point, parameter)


if  isfield(parameter, 'func_type') 
    switch parameter.func_type
         case {'heart'}
             ind = Locate_heart(point, parameter.center);
        case {'circle'}
             ind = Locate_circle(point, parameter.center, parameter.r);
        case {'ellipse'}
             ind = Locate_ellipse(point, parameter.center, parameter.a, parameter.b);
         case {'concentric_circle'}
             ind = Locate_concentric_circle(point, parameter.center, parameter.r, parameter.R);
         case {'concentric_circle_defect'}
             ind = Locate_concentric_circle_defect(point, parameter.center, parameter.r, parameter.R, parameter.theta, parameter.w);
        case {'capsule'}
             ind = Locate_capsule(point, parameter.center, parameter.l, parameter.w);
        case{'random_curv'}
             ind = Locate_random_curv(point, parameter.center);
        case{'four_flower'}
             ind = Locate_four_flower(point, parameter.center);
        case{'Dsquare'}
             ind = Locate_Dsquare(point, parameter.center,parameter.r);
        case{'polygon_random'}
             ind = Locate_polygon_random(point, parameter.center, parameter.x, parameter.y);
    end
else
    switch parameter.type
        case {'polygon'}

            if mod( parameter.num_side, 2 ) == 1 
                L = linspace(parameter.phi +  pi/2, parameter.phi + 2*pi + pi/2, parameter.num_side+1);
            else
                L = linspace(parameter.phi, parameter.phi + 2*pi, parameter.num_side+1);
            end
            parameter.xv = parameter.r * cos(L)';
            parameter.yv = parameter.r * sin(L)'; 

            ind = Locate_polygon(point, parameter.center, parameter.xv, parameter.yv);

        case { 'Darts'}

            % creat nodes 
            L1 = linspace(parameter.phi, 2*pi + parameter.phi , parameter.num + 1);
            L2 = linspace(parameter.phi + pi/parameter.num , 2*pi + parameter.phi + pi/parameter.num, parameter.num + 1);

            L = zeros(1, 2 * parameter.num +1);
            L((1:2:length(L))) = L1; L((2:2:length(L))) = L2(1:end-1);

%             L = [L1(1), L2(1), L1(2), L2(2), L1(3), L2(3), L1(4)];

            for i = 1:length(L)
                if mod(i, 2) == 1
                    parameter.xv(i) = parameter.r1*cos(L(i));
                    parameter.yv(i) = parameter.r1*sin(L(i));
                else
                    parameter.xv(i) = parameter.r2*cos(L(i));
                    parameter.yv(i) = parameter.r2*sin(L(i));
                end

            end 

            ind = Locate_Darts(point, parameter.center, parameter.xv, parameter.yv);
    end
        
      
end

end


function Point_ind = Locate_heart(point, center)
     Point_ind = [];
     for i = 1 : size(center, 1)
         Point_ind = union(Point_ind, Material_Locate_Heart(point,  center(i, :)));
     end
end


function Point_ind = Locate_circle(point, center, r)
     Point_ind = [];
     for i = 1 : size(center, 1)
         Point_ind = union(Point_ind, Material_Locate_circle(point, center(i, :), r));
     end
end

function Point_ind = Locate_ellipse(point, center, a, b)
     Point_ind = [];
     for i = 1 : size(center, 1)
         Point_ind = union(Point_ind, Material_Locate_ellipse(point, center(i, :), a, b));
     end
end

    
function Point_ind = Locate_concentric_circle(point, center, r, R)
     Point_ind = [];
     for i = 1 : size(center, 1)
         Point_ind = union(Point_ind, Material_Locate_concentric_circle(point, center(i, :), r, R));
     end
end

    
function Point_ind = Locate_polygon(point, center, xv, yv)
     Point_ind = [];
     for i = 1 : size(center, 1)
         Point_ind = union(Point_ind, Material_Locate_polygon(point, center(i, :), xv, yv));
     end
end


function Point_ind = Locate_Darts(point, center, xv, yv)
     Point_ind = [];
     for i = 1 : size(center, 1)
         Point_ind = union(Point_ind, Material_Locate_Darts(point, center(i, :), xv, yv));
     end
end

function Point_ind = Locate_concentric_circle_defect(point, center, r, R, theta, w)
     Point_ind = [];
     for i = 1 : size(center, 1)
         Point_ind = union(Point_ind, Material_Locate_concentric_circle_defect(point, center(i, :), r, R, theta, w));
     end
end

function Point_ind = Locate_capsule(point, center, l, w)
     Point_ind = [];
     for i = 1 : size(center, 1)
         Point_ind = union(Point_ind, Material_Locate_capsule(point, center(i, :), l, w));
     end
end

function Point_ind = Locate_random_curv(point, center)
     Point_ind = [];
     for i = 1 : size(center, 1)
         Point_ind = union(Point_ind, Material_Locate_random_curv(point, center(i, :)));
     end
end

function Point_ind = Locate_four_flower(point, center)
     Point_ind = [];
     for i = 1 : size(center, 1)
         Point_ind = union(Point_ind, Material_Locate_four_flower(point, center(i, :)));
     end
end

function Point_ind = Locate_Dsquare(point, center,r)
     Point_ind = [];
     for i = 1 : size(center, 1)
         Point_ind = union(Point_ind, Material_Locate_Dsquare(point, center(i, :),r));
     end
end
    
function Point_ind = Locate_polygon_random(point, center,x,y )
     Point_ind = [];
     for i = 1 : size(center, 1)
         Point_ind = union(Point_ind, Material_polygon_random(point, center(i, :), x, y ));
     end
end




