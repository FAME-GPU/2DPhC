function  ind = Material_Locate_polygon(point, center, xv, yv)


  point_relat = point - center;
  x = point_relat(:, 1); y = point_relat(:, 2);

  ind = inpolygon(x, y, xv, yv);

  ind = find ( ind );