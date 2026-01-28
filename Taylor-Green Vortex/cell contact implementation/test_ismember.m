mat = [1, 2; 
    3, 4; 
    5, 6];

seek1 = [3, 4];
seek2 = [3, 5];
seek3 = [2, 1];
seek4 = [7, 8];

m1 = ismember(mat, seek1, "rows");
m2 = ismember(mat, seek2, "rows");
m3 = ismember(mat, seek3, "rows");
m4 = ismember(mat, seek4, "rows");

mat(m1,:) = [];