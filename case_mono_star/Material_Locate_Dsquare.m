function  ind = Material_Locate_Dsquare(point, center,c)


  point_relat = point - center;
  x = point_relat(:, 1); y = point_relat(:, 2);

  node1 = find(  all([x>0, y>0],2));
  node2 = find(  all([x>0,y<=0],2));
  node3 = find( all([x<=0,y>0],2));
  node4 = find( all([x<=0,y<=0],2));

  t(node1) = atan(y(node1)./x(node1));
  t(node2) = 2*pi - abs(atan(abs(y(node2)) ./x(node2)));
  t(node3) = pi - atan(y(node3)./abs(x(node3)));
  t(node4) = pi+ abs(atan(y(node4)./x(node4))) ;

 
   ind1 = x.^2 + y.^2  <= ( c.*(0.2.*sin(5.*t)+1)).^2' ;
%      ind1 = x.^2 + y.^2  <= ( (0.2.*sin(5.*t)+1)).^2' ;
  ind2 = x.^2 + y.^2 >= (0.15.*(0.2.*sin(5.*t)+1)).^2';

  ind = find(all([ind1],2)) ;


