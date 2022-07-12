set terminal gif animate delay 10
set output 'animation.gif'
set size square
set xrange[0:1]
set yrange[0:1]
unset title
unset xtics
unset ytics

set view map
set pm3d map at b
#set contour base
unset colorbox
set cbrange[1:2]

load '/usr/local/gnuplot/jet.pal'

do for [i=0:200:10] {
 splot sprintf("./output/rho_%05d.dat", i) using 1:2:3 notitle w pm3d, \
       sprintf("./output/u_%05d.dat", i) using 1:2:(0):($3/20):($4/20):(0) every 20 notitle with vectors head filled lt 2 lc rgb 'white' front
}
