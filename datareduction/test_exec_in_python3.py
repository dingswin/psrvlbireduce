#!/usr/bin/env python
import numpy as np
def a_func():
    a = 'row_no'
    b = {'row_no':2}
    #exec(('%ss = np.array([])' % a), globals())
    exec(("%ss = b['%s']" % (a,a)), {'a':a, 'b':b, 'row_nos':None})
    return row_nos
a_func()
