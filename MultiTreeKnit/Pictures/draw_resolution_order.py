## install schemdraw to plot figure
import schemdraw
from schemdraw import dsp

with schemdraw.Drawing() as d:
    d.config(fontsize=14)
    d += (t12 := dsp.Box(w=2, h=2).anchor('W').label('$T_{1,2}$'))
    d += dsp.Arrow().at(t12.E).right(1)
    d += (t13 := dsp.Box(w=2, h=2).anchor('W').label('$T_{1,3}$'))
    d += dsp.Arrow().at(t13.E).right(1)
    d += (t14 := dsp.Box(w=2, h=2).anchor('W').label('$T_{1,4}$'))
    d += dsp.Arrow().at(t14.S).down(1)
    d += dsp.Arrow().at(t13.S).down(1)
    d += dsp.Line().at(t12.S).down(2)
    d += dsp.Arrow().right(2)
    d += (t23 := dsp.Box(w=2, h=2).anchor('W').label('$T_{2,3}$'))
    d += dsp.Arrow().at(t23.E).right(1)
    d += (t24 := dsp.Box(w=2, h=2).anchor('W').label('$T_{2,4}$'))
    d += dsp.Arrow().at(t24.S).down(1)
    d += dsp.Line().at(t23.S).down(2)
    d += dsp.Arrow().right(2)
    d += (t34 := dsp.Box(w=2, h=2).anchor('W').label('$T_{3,4}$'))

d.save("resolution_order.png")