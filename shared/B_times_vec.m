function  y = B_times_vec(x , muz_array, varargin)
    global  iter_count
      
       y = muz_array .* x;
       if nargin > 2
         tmp = sum(y,1);
         y = y - varargin{1} * tmp;
       end 
    iter_count = iter_count + 1;
end    
 

   