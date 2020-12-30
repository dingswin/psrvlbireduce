#!/usr/bin/env python
import multiprocessing, mspsrpifun
a = mspsrpifun.single_out_quasars_with__cone_searched_Gaia_sources__and__AgnCrossId('')
p1 = multiprocessing.Process(target=a.make_AgnCrossID_RA_Dec_table_step1, args=(0))
p2 = multiprocessing.Process(target=a.make_AgnCrossID_RA_Dec_table_step1, args=(1100))
p3 = multiprocessing.Process(target=a.make_AgnCrossID_RA_Dec_table_step1, args=(2200))

p1.start()
p2.start()
p3.start()

p1.join()
p2.join()
p3.join()
