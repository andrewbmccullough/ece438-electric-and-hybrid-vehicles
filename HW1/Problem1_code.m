us06 = xlsread('us06.xlsx');
plot(us06(:,1), us06(:,2));
title('US06 Drive Cycle');
xlabel('Time (s)');
ylabel('Speed (MPH)');