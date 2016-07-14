reset
f(x,y) = exp(-1/(2*(1-0.25)) \
         *((((x-0.5)**2)/0.25)+(((y-0.5)**2)/0.25)-((2*0.5*(x-0.5)*(y-0.5))/0.25)))
g(x,y) = exp(-1/(2*(1-0.25)) \
         *((((x)**2)/0.25)+(((y)**2)/0.25)-((2*0.5*(x)*(y))/0.25)))
h(x,y) = exp(-((7.05961 * (x-(-6.83351)) * (x-(-6.83351)))) \
         +(7.12101 * (x-(-6.83351)) * (y-(-47.1032))) \
         +(7.12101 * (x-(-6.83351)) * (y-(-47.1032))) \
         +(18.5449 * (y-(-47.1032)) * (y-(-47.1032))) / 2)
i(x,y) = exp(-((7204.8 * (x-(8.01099)) * (x-(8.01099))) \
         +(-173.735 * (x-(8.01099)) * (y-(-26.6039))) \
         +(-173.735 * (x-(8.01099)) * (y-(-26.6039))) \
         +(7.59913 * (y-(-26.6039)) * (y-(-26.6039)))) / 2)
set isosample 500
set contour
set cntrparam levels 50
set yrange[-18.5:-16.5]
set xrange[7:9]
splot i(x,y)
#unset surface
#set view 0,0
#replot