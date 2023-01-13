function  ind = Material_Locate_circle(point, center, r)


  point_relat = point - center;
  x = point_relat(:, 1); y = point_relat(:, 2);

  ind = find( x.^2 + y.^2 <= r.^2 );