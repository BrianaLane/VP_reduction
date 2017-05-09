# -*- coding: utf-8 -*-
from numpy import * 
from string import *
import os

#run this inside of the object folder (within outerdir/redux/subdir/data_folder/object)

#20170330
#sci1 = ['0233','0234','0235']
sci1 = ['0233']
sky1 = '0232'
tim1 = 'eve'

# sci2 = ['0236','0237','0238','0240','0241','0242']
# sky2 = '0239'
# tim2 = 'eve'

# sci3 = ['0243','0244','0245']
# sky3 = '0246'
# tim3 = 'morn'

#sci4 = ['0256','0257','0258']
sci4 = ['0257']
sky4 = '0255'
tim4 = 'morn'

# sci5 = ['0259','0260','0261']
# sky5 = '0262'
# tim5 = 'morn'

# sci6 = ['0264','0265','0266']
# sky6 = '0263'
# tim6 = 'morn'

# sci7 = ['0267','0268','0269']
# sky7 = '0270'
# tim7 = 'morn'

# #20170331
# sci1 = ['0359','0360','0361']
# sky1 = '0358'
# tim1 = 'eve'

# sci2 = ['0362','0363','0364','0366','0367','0368']
# sky2 = '0365'
# tim2 = 'eve'

# sci3 = ['0369','0370','0371','0373','0374','0375']
# sky3 = '0372'
# tim3 = 'morn'

# sci4 = ['0376','0377','0378']
# sky4 = '0379'
# tim4 = 'morn'

# sci5 = ['0381','0382','0383']
# sky5 = '0380'
# tim5 = 'morn'

# sci6 = ['0384','0385','0386']
# sky6 = '0387'
# tim6 = 'morn'

#20160403
# sci1 = ['0082','0083','0084']
# sky1 = '0081'
# tim1 = 'eve'

# sci2 = ['0085','0086','0087']
# sky2 = '0088'
# tim2 = 'eve'

# sci3 = ['0089','0090','0091']
# sky3 = '0088'
# tim3 = 'eve'

# sci4 = ['0092','0093','0094','0096','0097','0098']
# sky4 = '0095'
# tim4 = 'morn'

# sci5 = ['0099','0100','0101','0102','0103','0104','0105','0106','0107','0108','0109','0110','0111','0112','0113','0114']
# sky5 = '0108'
# tim5 = 'morn'

#20160404
# sci1 = ['0205','0206','0207']
# sca1 = [1.0,1.0,0.9]
# sky1 = '0204'
# tim1 = 'eve'

# sci2 = ['0208','0209','0210','0212','0213','0214']
# sca2 = [1.1,1.1,1.1,1.0,1.0,1.0]
# sky2 = '0211'
# tim2 = 'eve'

# sci3 = ['0215','0216','0217']
# sca3 = [1.1,1.0,1.0]
# sky3 = '0218'
# tim3 = 'morn'

# sci4 = ['0226','0227','0228']
# sca4 = [0.9,1.0,1.3]
# sky4 = '0225'
# tim4 = 'morn'

# sci5 = ['0229','0230','0231','0233','0234','0235']
# sca5 = [1.0,0.9,1.15,1.0,1.0,1.0]
# sky5 = '0232'
# tim5 = 'morn'

# sci6 = ['0236','0237','0238']
# sca6 = [1.0,1.0,1.0]
# sky6 = '0239'
# tim6 = 'morn'

#20160405
# sci1 = ['0320','0321','0322']
# sky1 = '0319'
# tim1 = 'eve'

# sci2 = ['0323','0324','0325','0327','0328','0329']
# sky2 = '0326'
# tim2 = 'eve'

# sci3 = ['0330','0331','0332']
# sky3 = '0333'
# tim3 = 'morn'

# sci4 = ['0341','0342','0343']
# sky4 = '0340'
# tim4 = 'morn'

# sci5 = ['0344','0345','0346']
# sky5 = '0347'
# tim5 = 'morn'

# sci6 = ['0349','0350','0351']
# sky6 = '0348'
# tim6 = 'morn'

# sci7 = ['0352','0353','0354']
# sky7 = '0355'
# tim7 = 'morn'

#20160406
# sci1 = ['0435','0436','0437']
# sky1 = '0434'
# tim1 = 'eve'

# sci2 = ['0438','0439','0440','0442','0443','0444']
# sky2 = '0441'
# tim2 = 'eve'

# sci3 = ['0445','0446','0447']
# sky3 = '0448'
# tim3 = 'morn'

#20160407
# sci1 = ['0528','0529','0530']
# sky1 = '0527'
# tim1 = 'eve'

# sci2 = ['0531','0532','0533','0535','0536','0537']
# sky2 = '0534'
# tim2 = 'eve'

# sci3 = ['0538','0539','0540']
# sky3 = '0541'
# tim3 = 'morn'

# sci4 = ['0549','0550','0551']
# sky4 = '0548'
# tim4 = 'morn'

# sci5 = ['0552','0553','0554']
# sky5 = '0555'
# tim5 = 'morn'

#20170228
# sci1 = ['0117','0118','0119']
# sky1 = '0116'
# tim1 = 'morn'

# sci2 = ['0120','0121','0122']
# sky2 = '0125'
# tim2 = 'morn'

# sci3 = ['0127','0128','0129']
# sky3 = '0126'
# tim3 = 'morn'

# sci4 = ['0130','0131','0132']
# sky4 = '0133'
# tim4 = 'morn'

#sci_frame = [sci1,sci2,sci3,sci4,sci5,sci6,sci7]
#sky_frame = [sky1,sky2,sky3,sky4,sky5,sky6,sky7]
#tim_frame = [tim1,tim2,tim3,tim4,tim5,tim6,tim7]
#sca_frame = [sca1,sca2,sca3,sca4,sca5,sca6,sca7]

sci_frame = [sci1,sci4]
sky_frame = [sky1,sky4]
tim_frame = [tim1,tim4]
#sca_frame = [sca1,sca2,sca3,sca4,sca5,sca6]

#If this equals True it runs CURE's parallel routine to run skysubtract in parallel. Else runs them one by one
parallel = True

for i in range(len(sci_frame)):

	if parallel:
		print "Running sky subtraction in parallel!"

		sci_list = ''
		for s in sci_frame[i]:
			sci_list += 'pesvp'+s+' '

		print "Reducing Sci Frame: "+sci_list
		print "Using Sky Frame: pesvp"+sky_frame[i]

		os.system("echo "+sci_list+" | tr ' ' '\n' | $CUREVP/parallel \"$CUREVP/subtractsky -f ../flat/"+tim_frame[i]+"/{}.fmod -d ../flat/"+tim_frame[i]+"/{}.dist -X "+tim_frame[i]+"/pesvp"+sky_frame[i]+".fits -F ../flat/"+tim_frame[i]+"/pesvp"+sky_frame[i]+".fmod -D ../flat/"+tim_frame[i]+"/pesvp"+sky_frame[i]+".dist "+tim_frame[i]+"/{}.fits\"")

	else: 
		print "Running sky subtraction individually"

		for s in range(len(sci_frame[i])):
			sci_im = 'pesvp'+sci_frame[i][s]

			print "Reducing Sci Frame: "+sci_im
			print "Using Sky Frame: pesvp"+sky_frame[i]

			print "$CUREVP/subtractsky -T 30 -w 250 --x-sky-scaling 1.0 -f ../flat/"+tim_frame[i]+"/"+sci_im+".fmod -d ../flat/"+tim_frame[i]+"/"+sci_im+".dist -X "+tim_frame[i]+"/pesvp"+sky_frame[i]+".fits -F ../flat/"+tim_frame[i]+"/pesvp"+sky_frame[i]+".fmod -D ../flat/"+tim_frame[i]+"/pesvp"+sky_frame[i]+".dist "+tim_frame[i]+"/"+sci_im+".fits"
			#os.system("$CUREVP/subtractsky -T 30 -w 250 --x-sky-scaling "+str(sca_frame[i][s])+" -f ../flat/"+tim_frame[i]+"/"+sci_im+".fmod -d ../flat/"+tim_frame[i]+"/"+sci_im+".dist -X "+tim_frame[i]+"/pesvp"+sky_frame[i]+".fits -F ../flat/"+tim_frame[i]+"/pesvp"+sky_frame[i]+".fmod -D ../flat/"+tim_frame[i]+"/pesvp"+sky_frame[i]+".dist "+tim_frame[i]+"/"+sci_im+".fits")









