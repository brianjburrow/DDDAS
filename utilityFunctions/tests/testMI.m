clear all
close all

%% Test function that describes the use of the multiIndex functions
% that are stored in the utility folder

polyOrder = 4;
dim = 1;

a = genTotalOrderMI(polyOrder, dim);

b = genNoMixedMI(polyOrder, dim);

c = genNoCrossMI(polyOrder, dim);
fprintf("TotalOrder: Size %d \n", length(a));

fprintf("No Mixed Terms: Size %d \n", length(b));

fprintf("No Cross Terms: Size %d \n", length(c));
