#------------------------------------------------------------------------------#
#Description: model of rat prelimbic Layer 5 pyramidal neuron with distance-based channel conductances
#Reference:
#	1. Cell morphologies are from  archive: Amaral from neuromorpho.org
#	2. Channel distributions are from Narayanan et al., 2010 http://www.ncbi.nlm.nih.gov/pubmed/20457854
#
#Version:
#	v1.0 dummy channels for calcium- L-type T-type and N-type and km; 
#	
#------------------------------------------------------------------------------#
import moose
#import moogli
import pylab
import numpy
import time
from moose import neuroml
from moose.neuroml.ChannelML import ChannelML
import sys
import os

#this path must be changed acc to user settings
sys.path.append('/home/sarathy/Work/moose/moose-core/python/rdesigneur') 
import rdesigneur as rd

useGssa = True
combineSegments = True

#specifying the .swc file names
elecFileNames = ("c10861.CNG.swc", "cell5zr.CNG.swc") #only 4 and 5 in v1


PI = 3.14159265359

frameRunTime = 0.1e-2#delta time
baselineTime = 0.01#delay
currPulseTime = 2#duration
postPulseTime = 0.1
runtime = baselineTime + currPulseTime + postPulseTime

injectionCurrent = 1.5e-9
deltaCurrent = 0.1e-9
maxCurrent = 0.7e-9

somaVm = []
spineVm = []
iList = []

def buildRdesigneur():
	path = os.path.abspath(os.path.dirname('currentStep_CA3PC-Narayanan2010_v1.0.py'))
	cellProto = [ [ "./cells/" + elecFileNames[0], "elec" ] ]
	chanProto = [
	    [path+'/proto20.make_Ca_conc()', 'Ca_conc'], \
        ['./chans/CaL.xml'], \
        ['./chans/CaN.xml'], \
        ['./chans/CaT.xml'], \
        ['./chans/hd.xml'], \
        ['./chans/kad.xml'], \
        [path+'/proto20.make_K_AHP()', 'kahp'], \
        ['./chans/kap.xml'], \
        [path+'/proto20.make_K_C()', 'kca'], \
        ['./chans/kdr.xml'], \
        ['./chans/km.xml'], \
        ['./chans/na3.xml'] \
    ]
	spineProto = [\
		['makeSpineProto()', 'spine']
	]
	chemProto = []
	passiveDistrib = [ 
	        [ ".", "#", "RM", "6", "CM", "0.0075", "RA", "2",  \
	            "Em", "-65e-3", "initVm", "-65e-3" ] \
    ]
        chanDistrib = [ \
	        ["Ca_conc", "#soma#,#dend#", "tau", "0.0133" ], \
	        ["CaL", "#dend#,#soma#", "Gbar", "p < 50e-6 ? 25 : 0" ], \
            ["CaN", "#soma#,#dend#", "Gbar", "28" ], \
            ["CaT", "#soma#,#dend#", "Gbar", "2.5" ], \
            ["hd", "#dend#,#soma#", "Gbar", "0.1*(1+1.75*p/100)" ], \
            ["kad", "#dend#", "Gbar", "p >= 100e-6 ? 1e-3*(7+(11*p)/100) : 0" ], \
            ["kad", "#soma#", "Gbar", "70" ], \
            ["kahp", "#", "Gbar", "3" ], \
            ["kap", "#dend#", "Gbar", "p < 100e-6 ? 1e-3*(7+(11*p)/100 : 0" ], \
            ["kap", "#soma#", "Gbar", "70" ], \
            ["kca", "#soma#,#dend#", "Gbar", "8" ], \
            #["kdr", "#soma#,#dend#", "Gbar", "300"], \
            ["kdr", "#soma#,#dend#", "Gbar", "88"], \
            #["M", "#soma#,#dend#", "Gbar", "p < 100e-6 ? 1 : 0"  ], \
            ["M", "#soma#,#dend#", "Gbar", "p < 100e-6 ? 0.1 : 0"  ], \
            #["na3", "#soma#,#dend#", "Gbar", "250" ], \
            ["na3", "#soma#,#dend#", "Gbar", "180" ] \
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
		#print "1. injected current = ", injectionCurrent
		injectionCurrent += deltaCurrent
		#print "del stim first ", moose.element('/clock').currentTime
	if numpy.fabs( currTime - baselineTime - currPulseTime) < frameRunTime/2.0 :
		#end
		eList = moose.wildcardFind( '/model/elec/#soma#' )
		assert( len(eList) > 0 )
		eList[0].inject = 0.0
		#print "2. injected current = ", injectionCurrent
		#print "del stim second ", moose.element('/clock').currentTime
	if runtime - currTime < frameRunTime * 2.0 :
		#print "3. reinit-ing"
		somaVm.append( moose.element( '/graphs/VmTab' ).vector )
		spineVm.append( moose.element( '/graphs/eSpineVmTab' ).vector )
		iList.append(injectionCurrent)
		if injectionCurrent < maxCurrent :
			moose.reinit()	

def displayPlots():
    pylab.figure(1, figsize= (8,10))
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
    
    
def chanDisPlot():
    pylab.figure(1, figsize= (8,10))
    ax = pylab.subplot( 1,1,1 )
    neuron = moose.element( '/model/elec' )
    comptDistance = dict( zip( neuron.compartments, neuron.pathDistanceFromSoma ) )
    for i in moose.wildcardFind( '/library/#[ISA=ChanBase]' ):
        chans = moose.wildcardFind( '/model/elec/#/' + i.name )
        print i.name, len( chans )
        p = [ 1e6*comptDistance.get( j.parent, 0) for j in chans ]
        Gbar = [ j.Gbar/(j.parent.length * j.parent.diameter * PI) for j in chans ]
        if len( p ) > 2:
            pylab.plot( p, Gbar, linestyle = 'None', marker = ".", label = i.name )
    ax.set_yscale( 'log' )
    pylab.xlabel( "Distance from soma (microns)" )
    pylab.ylabel( "Channel density (Seimens/sq mtr)" )
    pylab.legend()
    pylab.title( 'Channel distribution' )
    pylab.show()
    
def dontBuild3dDisplay( rdes ):
    currTime = moose.element( '/clock' ).currentTime
    while currTime < runtime:
        moose.start( frameRunTime )
        #print 'in dont build while:', moose.element('/clock').currentTime
        currTime = moose.element( '/clock' ).currentTime
        if ( currTime < runtime ):
			#print "in dont build if "
			deliverStim( currTime )
        else:
			displayPlots()
			#chanDisPlot()
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
	
	
	
	
	
	
	
	
