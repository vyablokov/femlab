set pm3d;
#set palette defined (50 "blue", 125 "green", 200 "red");
set hidden3d;
unset surface;
set view map;
splot 'points.txt';
#pause 0.001
pause mouse "Click to exit.";
