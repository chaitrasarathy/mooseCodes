#------------------------------------------------------------------------------#
#This is version 1 of moose code for a detailed model of rat CA3 pyramidal neuron
#Channel morphologies are from Barrionuevo G. et al., 1996
# http://www.ncbi.nlm.nih.gov/pubmed/8743416
#
#Channel distributions from Lazarewicz 2002 et al., 
#   http://www.ncbi.nlm.nih.gov/pubmed/12459292
#------------------------------------------------------------------------------#
import moose
#import moogli
import pylab
import numpy
import time
from moose import neuroml
from moose.neuroml.ChannelML import ChannelML
import sys

#this path must be changed acc to user settings
sys.path.append('/home/sarathy/Work/moose/moose-core/python/rdesigneur') 
import rdesigneur as rd

useGssa = True
combineSegments = True

#specifying the .swc file names
elecFileNames = ("cell4zr.CNG.swc", "cell5zr.CNG.swc") #only 4 and 5 in v1


PI = 3.14159265359

frameRunTime = 1e-2
baselineTime = 0.05
currPulseTime = 0.05
postPulseTime = 0.1
runtime = baselineTime + currPulseTime + postPulseTime

injectionCurrent = 0.8e-9
deltaCurrent = 0.1e-9
maxCurrent = 1.2e-9

somaVm = []
spineVm = []
iList = []

def buildRdesigneur():
	cellProto = [ [ "./cells/" + elecFileNames[0], "elec" ] ]
	chanProto = [
        ['./chans/caL.xml'], \
        #['./chans/caN.xml'], \
        #['./chans/caT.xml'], \
        ['./chans/hd.xml'], \
        ['./chans/kad.xml'], \
        #['./chans/kahp.xml'], \
        ['./chans/kap.xml'], \
        #['./chans/kc.xml'], \
#        ['./chans/kd.xml'], \ #slow subthreshold potassium current
        ['./chans/kdr.xml'], \
        ['./chans/km.xml'], \
        ['./chans/na3.xml'], \
    ]
	spineProto = [\
		['makeSpineProto()', 'spine']
	]
	chemProto = []
	passiveDistrib = [ 
	        [ ".", "#", "RM", "3", "CM", "0.02", "RA", "2",  \
	            "Em", "-65e-3", "initVm", "-65e-3" ], \
            [ ".", "#soma#", "RM", "6", "CM", "0.01", "RA", "2" ] \
        ]
	chanDistrib = [ \
	        ["L", "#dend#,#soma#", "Gbar", "p < 50 ? 13 : 0" ], \
            #["Cav2.2", "#soma#,#dend#", "Gbar", "15" ], \
            #["CaT", "#soma#,#dend#", "Gbar", "10" ], \
            ["hd", "#dend#,#soma#", "Gbar", "1*(1+3/100*p)" ], \
            ["kad", "#soma#,#dend#", "Gbar", "p > 100e-6 ? 110*(1+(11/(7*100))*p) : 0" ], \
            #["kahp", "#apical#", "Gbar", "10" ], \
            #["kahp", "#soma#", "Gbar", "0" ], \
            #["kahp", "#basal", "Gbar", "5" ], \
            ["kap", "#soma#,#dend#", "Gbar", "p <= 100e-6 ? 110*(1+(11/(7*100))*p) : 0" ], \
            #["KC", "#soma#,#dend#", "Gbar", "p < 150e-6 ? 40*(150-p)/150 : 0" ], \
           #["kd", "#", "Gbar", "p < 50e-6 ? 500 : 100" ], \
            ["kdr", "#soma#,#dend#", "Gbar", "20"], \
            ["M", "#soma#,#dend#", "Gbar", "p <= 100e-6 ? 1 : 0"  ], \
            ["na3", "#soma#,#dend#", "Gbar", "350" ], \
        ]
	spineDistrib = [ \
		["spine", '#apical#', "spineSpacing", "20e-6", \
                "spineSpacingDistrib", "2e-6", \
                "angle", "0", \
                "angleDistrib", str( 2*PI ), \
                "size", "1", \
                "sizeDistrib", "0.5" ] \
        ]
	chemDistrib = []
	adaptorList = []
	rd.addSpineProto()
	rdes = rd.rdesigneur(
        useGssa = useGssa, \
        combineSegments = combineSegments, \
        stealCellFromLibrary = True, \
        passiveDistrib = passiveDistrib, \
        spineDistrib = spineDistrib, \
        chanDistrib = chanDistrib, \
        chemDistrib = chemDistrib, \
        cellProto = cellProto, \
        chanProto = chanProto, \
        chemProto = chemProto, \
        adaptorList = adaptorList
    )
	return rdes

def deliverStim(currTime):
	global injectionCurrent	
	global spineVm
	global somaVm
	if numpy.fabs( currTime - baselineTime ) < frameRunTime/2.0 :
		#start
		eList = moose.wildcardFind( '/model/elec/#soma#' )
		assert( len(eList) > 0 )
		eList[0].inject = injectionCurrent
		injectionCurrent += deltaCurrent
	if numpy.fabs( currTime - baselineTime - currPulseTime) < frameRunTime/2.0 :
		#end
		eList = moose.wildcardFind( '/model/elec/#soma#' )
		assert( len(eList) > 0 )
		eList[0].inject = 0.0
	if runtime - currTime < frameRunTime * 2.0 :
		somaVm.append( moose.element( '/graphs/VmTab' ).vector )
		spineVm.append( moose.element( '/graphs/eSpineVmTab' ).vector )
		iList.append(injectionCurrent)
		if injectionCurrent < maxCurrent :
			moose.reinit()	

def displayPlots():
    pylab.subplot( 2,1,1 )
    t = numpy.arange( 0, len( somaVm[0] ), 1 ) * 50e-6
    for i in somaVm:
        pylab.plot( t, i[:len(t)] )
    pylab.legend()
    pylab.title( 'somaVm' )
    pylab.subplot( 2,1,2 )
    t = numpy.arange( 0, len( spineVm[0] ), 1 ) * 50e-6
    for i in spineVm:
        pylab.plot( t, i[:len(t)] )
    pylab.legend()
    pylab.title( 'spineVm' )
    pylab.show()
    
def dontBuild3dDisplay( rdes ):
    currTime = moose.element( '/clock' ).currentTime
    while currTime < runtime:
        moose.start( frameRunTime )
        currTime = moose.element( '/clock' ).currentTime
        if ( currTime < runtime ):
            deliverStim( currTime )
        else:
			displayPlots()
			break
    
def buildPlots( rdes ):
    if not moose.exists( '/graphs' ):
        moose.Neutral( '/graphs' )
    numPlots = 10
    vtab = moose.Table( '/graphs/VmTab' )
    moose.connect( vtab, 'requestOut', rdes.soma, 'getVm' )
    eSpineCaTab = moose.Table( '/graphs/eSpineCaTab' )
    #moose.le( '/model/elec/head50' )
    elist = moose.wildcardFind( '/model/elec/head#' )
    numSpines = len( elist )
    assert( numSpines > 0 )
    path = elist[ numSpines / 2 ].path + "/Ca_conc"
    moose.connect( eSpineCaTab, 'requestOut', path, 'getCa' )
    eSpineVmTab = moose.Table( '/graphs/eSpineVmTab' )
    moose.connect( eSpineVmTab, 'requestOut', elist[ numSpines / 2], 'getVm' )
    eSpineGkTab = moose.Table( '/graphs/eSpineGkTab' )
    path = elist[ numSpines / 2 ].path + "/NMDA"
    moose.connect( eSpineGkTab, 'requestOut', path, 'getGk' )
    
def main():
	global synSpineList
	global synDendList
	numpy.random.seed( 1234 )
	rdes = buildRdesigneur()
	rdes.buildModel( '/model' )
	assert( moose.exists( '/model' ) )
	synSpineList = moose.wildcardFind( "/model/elec/#head#/glu,/model/elec/#head#/NMDA" )
	temp = set( moose.wildcardFind( "/model/elec/#/glu,/model/elec/#/NMDA" ) )
	synDendList = list( temp - set( synSpineList ) )
	print "num spine, dend syns = ", len( synSpineList ), len( synDendList )
	moose.reinit()
	buildPlots( rdes )
	t1 = time.time()
	dontBuild3dDisplay(rdes)
	print 'real time = ', time.time() - t1

if __name__ == '__main__':
	main()
	
	
	
	
	
	
	
	
