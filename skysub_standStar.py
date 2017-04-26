# -*- coding: utf-8 -*-
from numpy import * 
from string import *
import glob
import os

#20160404
# sci = ['0219','0220','0221','0222','0223','0224']
# tim = 'morn'

#20160406
# sci = ['0449','0450','0451','0452','0453','0454']
# tim = 'morn'

#20160407
sci = ['0542','0543','0544','0545','0546','0547']
tim = 'morn'

#20160403
# sci = ['0102','0103','0104','0105','0106','0107']
# tim = 'morn'

#20160405
# sci = ['0334','0335','0336','0337','0338','0339']
# tim = 'morn'

sci_list = ''
for s in sci:
	sci_list += 'pesvp'+s+' '

print "Reducing Sci Frame: "+sci_list

os.system("echo "+sci_list+" | tr ' ' '\n' | $CUREVP/parallel \"$CUREVP/subtractsky -f ../flat/"+tim+"/{}.fmod -d ../flat/"+tim+"/{}.dist "+tim+"/{}.fits\"")
