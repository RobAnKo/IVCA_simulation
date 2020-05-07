function c = colony_count_from_area(area)
%c = (area - 76.896)/89.517;
c = (0.0079.*area)+0.2572;
c(c<1) = 1;
end

