#------------------------------------------------------------------------------#
#This is version 1.1 of moose code for a detailed model of rat CA3 pyramidal neuron
#Channel morphologies are from Barrionuevo G. et al., 1996
# http://www.ncbi.nlm.nih.gov/pubmed/8743416
#
#Channel distributions for soma and densrites are from Lazarewicz 2002 et al., 
#   http://www.ncbi.nlm.nih.gov/pubmed/12459292
#Channel distributions for axon are from Hemmond 2008 et al., 
#   http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4339291/
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
elecFileNames = ("cell6zr.CNG.swc", "cell5zr.CNG.swc") #only 4 and 5 in v1


PI = 3.14159265359

frameRunTime = 1e-2
baselineTime = 0.05
currPulseTime = 0.05
postPulseTime = 0.1
runtime = baselineTime + currPulseTime + postPulseTime

injectionCurrent = 0.1e-9
deltaCurrent = 0.1e-9
maxCurrent = 0.6e-9

somaVm = []
spineVm = []
iList = []

#dummy channels for calcium- L-type T-type and N-type, kahp, kc; no kd; 


def buildRdesigneur():
	cellProto = [ [ "./cells/" + elecFileNames[0], "elec" ] ]
	#chanProto = [
	    #['/home/sarathy/Work/templates/neurons/cA3_PC_rat/chans/proto20.make_Ca_conc()', 'Ca_conc'], \
        #['./chans/CaL.xml'], \
        #['./chans/CaN.xml'], \
        #['./chans/CaT.xml'], \
        #['./chans/hd.xml'], \
        #['./chans/kad.xml'], \
        #['/home/sarathy/Work/templates/neurons/cA3_PC_rat/chans/proto20.make_K_AHP()', 'kahp'], \
        #['./chans/kap.xml'], \
        #['/home/sarathy/Work/templates/neurons/cA3_PC_rat/chans/proto20.make_K_C()', 'kca'], \
##        ['./chans/kd.xml'], \ #slow subthreshold potassium current
        #['./chans/kdr.xml'], \
        #['./chans/km.xml'], \
        #['./chans/na3.xml'], \
    #]
        chanProto = [
        ['/home/sarathy/Work/templates/neurons/ca3_PC_rat/chans/proto20.make_Ca_conc()', 'Ca_conc'], \
        ['./chans/CaL.xml'], \
        ['./chans/CaN.xml'], \
        ['./chans/CaT.xml'], \
        ['./chans/hd-90.xml'], \
        ['./chans/hd-82.xml'], \
        ['./chans/kad.xml'], \
        ['/home/sarathy/Work/templates/neurons/ca3_PC_rat/chans/proto20.make_K_AHP()', 'kahp'], \
        ['./chans/kap.xml'], \
        ['/home/sarathy/Work/templates/neurons/ca3_PC_rat/chans/proto20.make_K_C()', 'kca'], \
#        ['./chans/kd.xml'], \ #slow subthreshold potassium current
        ['./chans/kdr.xml'], \
        ['./chans/km.xml'], \
        ['./chans/na3.xml'], \
    ]
	spineProto = [\
		['makeSpineProto()', 'spine']
	]
	chemProto = []
	#------------------------Lazarewicz et al 2002----------------------------------------#
	passiveDistrib = [ 
	        [ ".", "#", "RM", "3", "CM", "0.02", "RA", "2",  \
	            "Em", "-65e-3", "initVm", "-65e-3" ], \
            [ ".", "#soma#", "RM", "6", "CM", "0.01" ] \
        ]
	chanDistrib = [ \
	        ["Ca_conc", "#soma#,#dend#", "tau", "0.0133" ], \
	        ["CaL", "#soma#", "Gbar", "13" ], \
            ["CaL", "#dend#,#apical#,#basal#", "Gbar", "p < 50e-6 ? 13 : 0" ], \
            ["CaN", "#soma#,#dend#", "Gbar", "15" ], \
            ["CaT", "#soma#,#dend#", "Gbar", "10" ], \
            ["hd-82", "#dend#,#apical#", "Gbar", "p > 100 ? 1*(1+3/100*p) : 0" ], \
            ["hd-90", "#dend#,#basal#", "Gbar", "p <= 100 ? 1*(1+3/100*p) : 0" ], \
            ["hd-90", "#soma#", "Gbar", "1" ], \
            ["kad", "#dend#,#apical#,#basal#", "Gbar", "p > 100e-6 ? 110*(1+(11/(7*100))*p) : 0" ], \
            ["kad", "#soma#", "Gbar", "110" ], \
            ["kahp", "#soma#", "Gbar", "0" ], \
            ["kahp", "#dend#,#apical#", "Gbar", "5" ], \
            ["kahp", "#dend#,#basal#", "Gbar", "10" ], \
            ["kap", "#dend#,#apical#,#basal#", "Gbar", "p > 100e-6 ? 110*(1+(11/(7*100))*p) : 0" ], \
            ["kap", "#soma#", "Gbar", "110" ], \
            ["kap", "#axon#", "Gbar", "200"], \
            ["kca", "#dend#,#apical#,#basal#", "Gbar", "p < 150e-6 ? 40*(150-p)/150 : 0" ], \
			["kca", "#soma#", "Gbar", "40"], \
			##["kd", "#", "Gbar", "p < 50e-6 ? 500 : 100" ], \
            ["kdr", "#soma#,#dend#", "Gbar", "100"], \
            ["kdr", "#axon#", "Gbar", "100"], \
            ["M", "#soma#,#dend#,#apical#,#basal#", "Gbar", "p <= 100e-6 ? 1 : 0"  ], \
            ["na3", "#soma#,#dend#", "Gbar", "350" ], \
            ["na3", "#axon#", "Gbar", "1100" ], \
            #["na3", "#axon#", "Gbar", "1100" ], \
        ]
        #------------------------Narayanan et al 2010----------------------------------------#
    	#passiveDistrib = [ 
	        #[ ".", "#", "RM", "6", "CM", "0.0075", "RA", "2",  \
	            #"Em", "-65e-3", "initVm", "-65e-3" ], \
            ##[ ".", "#soma#", "RM", "6", "CM", "0.01", "RA", "2" ] \
        #]
        #chanDistrib = [ \
	        #["Ca_conc", "#soma#,#dend#", "tau", "0.0133" ], \
	        #["CaL", "#dend#,#soma#", "Gbar", "p < 50e-6 ? 25 : 0" ], \
            #["CaN", "#soma#,#dend#", "Gbar", "28" ], \
            #["CaT", "#soma#,#dend#", "Gbar", "2.5" ], \
            #["hd", "#dend#,#soma#", "Gbar", "0.1*(1+1.75*p/100)" ], \
            #["kad", "#dend#", "Gbar", "p >= 100e-6 ? 1e-3*(7+(11*p)/100) : 0" ], \
            #["kad", "#soma#", "Gbar", "70" ], \
            #["kahp", "#", "Gbar", "3" ], \
            #["kap", "#dend#", "Gbar", "p < 100e-6 ? 1e-3*(7+(11*p)/100 : 0" ], \
            #["kap", "#soma#", "Gbar", "70" ], \
            #["kca", "#soma#,#dend#", "Gbar", "8" ], \
           ##["kd", "#", "Gbar", "p < 50e-6 ? 500 : 100" ], \
            #["kdr", "#soma#,#dend#", "Gbar", "300"], \
            #["M", "#soma#,#dend#", "Gbar", "p < 100e-6 ? 1 : 0"  ], \
            #["na3", "#soma#,#dend#", "Gbar", "250" ], \
            ##["na3", "#axon#", "Gbar", "p < 80e-6 ? 5000 : 0" ], \
        #]
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
		print "1. injected current = ", injectionCurrent
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
        currTime = moose.element( '/clock' ).currentTime
        if ( currTime < runtime ):
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
	
	
	
	
	
	
	
	
