'''
--------------------------------------------------------------------------------
monte carlo simulation processing functions
--------------------------------------------------------------------------------
1. read_sim_header
2. read_sim_data
3. combine_sim_data_s_p_ratio
--------------------------------------------------------------------------------
'''



################################################################################
# 1. read_sim_header
################################################################################
def read_sim_header(simfilename):
	'''
	----------------------------------------------------------------------------
	read Monte Carlo simulation header file
	----------------------------------------------------------------------------
	function read_sim_header reads the header output of the Monte Carlo 
	simulations and	stores the simulation/source/structural parameters in a 
	python dictionary called simpa.
	----------------------------------------------------------------------------
	:param simfilename: 	name of the simulation header
	----------------------------------------------------------------------------
	:return simpa: 			dictionary containing simulation, source and 
							structural parameters
	----------------------------------------------------------------------------
	'''
	import glob																	# import glob module
	#---------------------------------------------------------------------------
	file = glob.glob(simfilename)												# glob search for filename
	simheader   	 = open(file[0]) 											# open the simulation header
	lines         	 = simheader.readlines()									# read the lines of the header file
	simpa       	 = {}														# create empty dictionary for storing the parameters
	#---------------------------------------------------------------------------
	# Simulation parameters
	#---------------------------------------------------------------------------
	simpa['tag'] 		      =    str(lines[1].split()[-1])					# Simulation Tag
	simpa['npart'] 	  	      =    int(lines[2].split()[-1])					# Number of Particles						 
	simpa['q'] 		 	      =    int(lines[3].split()[-1])					# Number of Intervals in probability table			
	simpa['Nangle'] 	      =    int(lines[4].split()[-1])					# Number of Intervals in (0 PI/2) for reflection/conversion tables	
	simpa['GridSize'] 	      = [int(lines[5].split()[-5]),	
                                 int(lines[5].split()[-4]),					# Grid size (X, y, z, theta, t)		
								 int(lines[5].split()[-3]),					# Grid size (x, Y, z, theta,t)
                                 int(lines[5].split()[-2]),					# Grid size (x, y, Z, theta,t)
                                 int(lines[5].split()[-1])]					# Grid size (x, y, z, theta,T)
	simpa['GridSpacing']      = [float(lines[6].split()[-5]),
                                 float(lines[6].split()[-4]),					# Grid spacing (X, y, z, theta,t) in km or s
								 float(lines[6].split()[-3]),					# Grid spacing (x, Y, z, theta,t) in km or s
								 float(lines[6].split()[-2]),					# Grid spacing (x, y, Z, theta,t) in km or s
								 float(lines[6].split()[-1])]					# Grid spacing (x, y, z, theta,T) in km or s
	simpa['TimeStep'] 	      =  float(lines[7].split()[-1])					# Simulation time interval in s
	simpa['SimLen'] 	      =  float(lines[8].split()[-1])					# Simulation length in s								
	simpa['Qref'] 	          =  float(lines[9].split()[-1])					# Reference 1/Q (internal amplification)					
	simpa['emodulo'] 	  	  =  float(lines[10].split()[-1])					# Energy modulo		
	simpa['SimulationTime']   =  float(lines[11].split()[-1])					# Execution time in seconds									
	simpa['looptime'] 	      =  float(lines[12].split()[-1])					# Simulation time in seconds											
	simpa['outputtype']       =    int(lines[13].split()[-1])					# Output type (0 for 4D, 1 for section)'	
	#---------------------------------------------------------------------------
	# Source Parameters
	#---------------------------------------------------------------------------
	simpa['SourcePos']        = [float(lines[16].split()[-3]),					# Source position (X, y, z) in km from grid origin
								 float(lines[16].split()[-2]),					# Source position (x, Y, z) in km from grid origin
								 float(lines[16].split()[-1])]					# Source position (x, y, Z) in km from grid origin
	simpa['SourceType']       =    int(lines[17].split()[-1])					# Source type (0 for P, 1 for S-Source)       			
	simpa['SourceOri']        = [float(lines[18].split()[-3]),					# Source orientation (STRIKE, dip, rake) in degrees											
								 float(lines[18].split()[-2]),					# Source orientation (strike, DIP, rake) in degrees											
								 float(lines[18].split()[-1])]					# Source orientation (strike, dip, RAKE) in degrees											
	simpa['frc'] 		      =  float(lines[19].split()[-1])					# Angular frequency			
	simpa['sourcebody']       =    int(lines[20].split()[-1])					# Source body	
	#---------------------------------------------------------------------------
	# Structural settings			
	#---------------------------------------------------------------------------
	simpa['Nbody']            =    int(lines[23].split()[-1])					# Number of bodies
	simpa['body'] 			  =  {}												# create empty body dictionary
	#---------------------------------------------------------------------------
	for i in range(simpa['Nbody']):												# loop over all the bodies
		simpa['body'][i] = {}													# create dictionary for body[i]
		simpa['body'][i]['x'] 	 		=  [float(lines[25+i*32].split()[-2]),	# X coordinate 1
										    float(lines[25+i*32].split()[-1])]	# X coordinate 2
		simpa['body'][i]['y']      		 = [float(lines[26+i*32].split()[-2]),	# Y coordinate 1	
											float(lines[26+i*32].split()[-1])]	# Y coordinate 2	
		simpa['body'][i]['z']        	=  [float(lines[27+i*32].split()[-2]),	# Z coordinate 1
											float(lines[27+i*32].split()[-1])]	# Z coordinate 2	
		simpa['body'][i]['TimStepFrac'] =   float(lines[28+i*32].split()[-1])	# Time Step fraction
		simpa['body'][i]['v']         	=[0,float(lines[29+i*32].split()[-2])]	# S Velocity						
		simpa['body'][i]['gam']        	=   float(lines[29+i*32].split()[-1])	# gamma				
		simpa['body'][i]['rho']        	=   float(lines[30+i*32].split()[-1])	# Density						
		simpa['body'][i]['type']      	=     str(lines[31+i*32].split()[-1])	# Type of ACF (gauss, expo, karman)						
		simpa['body'][i]['kap']         =   float(lines[32+i*32].split()[-1])	# Kappa (only evaluated for karman ACF) 						
		simpa['body'][i]['a']       	=   float(lines[33+i*32].split()[-1])	# Correlation distance a						
		simpa['body'][i]['e']         	=   float(lines[34+i*32].split()[-1])	# Fractional fluctutation e 						
		simpa['body'][i]['ny']         	=   float(lines[35+i*32].split()[-1])	# Birch ny 							
		simpa['body'][i]['Qi'] 			=  [float(lines[36+i*32].split()[-2]),	# Intrinsic Attenuation (1/Q_P 1/Q_s)
											float(lines[36+i*32].split()[-1])]	# Intrinsic Attenuation (1/Q_p 1/Q_S)					
		simpa['body'][i]['Qex'] 	    =  [float(lines[37+i*32].split()[-2]),	# Exponent of absobtion factor of single timestep	
											float(lines[37+i*32].split()[-1])]  # Exponent of absobtion factor of single timestep		
		simpa['body'][i]['pg'] 	 		=  [float(lines[38+i*32].split()[-2]),	# Inverse of free path of P- and S-waves
											float(lines[38+i*32].split()[-1])]  # Inverse of free path of P- and S-waves		
		simpa['body'][i]['tau'] 	    =  [float(lines[39+i*32].split()[-2]),	# Mean free times  of P- and S-waves 
											float(lines[39+i*32].split()[-1])]  # Mean free times  of P- and S-waves 		
		simpa['body'][i]['path'] 	  	=  [float(lines[40+i*32].split()[-2]),	# Length of pathes in one time step for P- and S-waves 
											float(lines[40+i*32].split()[-1])]  # Length of pathes in one time step for P- and S-waves 
		#------------------------------------------------------------------------
		simpa['body'][i]['scattab']     	 =   {}								# empty dictionary for scattering parameters
		simpa['body'][i]['scattab']['gpp0']  = float(lines[42+i*32].split()[-1])# 'gpp0'
		simpa['body'][i]['scattab']['gps0']  = float(lines[43+i*32].split()[-1])# 'gps0'			
		simpa['body'][i]['scattab']['gsp0']  = float(lines[44+i*32].split()[-1])# 'gsp0'		
		simpa['body'][i]['scattab']['gssl0'] = float(lines[45+i*32].split()[-1])# 'gssl0'						
		simpa['body'][i]['scattab']['gssr0'] = float(lines[46+i*32].split()[-1])# 'gssr0'			
		simpa['body'][i]['scattab']['t_pp']  = float(lines[48+i*32].split()[-1])# '<cos(theta)>_pp'			
		simpa['body'][i]['scattab']['t_ps']  = float(lines[49+i*32].split()[-1])# '<cos(theta)>_ps'			
		simpa['body'][i]['scattab']['t_sp']  = float(lines[50+i*32].split()[-1])# '<cos(theta)>_sp'					
		simpa['body'][i]['scattab']['t_ssl'] = float(lines[51+i*32].split()[-1])# '<cos(theta)>_ssl'						
		simpa['body'][i]['scattab']['t_ssr'] = float(lines[52+i*32].split()[-1])# '<cos(theta)>_ssr'
		simpa['body'][i]['scattab']['l^*_p'] = float(lines[54+i*32].split()[-1])# 'l^*_p'					
		simpa['body'][i]['scattab']['l^*_s'] = float(lines[55+i*32].split()[-1])# 'l^*_s'					
		simpa['body'][i]['v'][0] 			 = (simpa['body'][i]['v'][1] * 		# 'P velocity'	
												simpa['body'][i]['gam'])		# 'P velocity'		
		simpa['body'][i]['l']   	= simpa['frc'] / simpa['body'][i]['v'][1]	# 'l'						
		simpa['body'][i]['q']       = simpa['q']								# 'body q'							
		simpa['body'][i]['path'] = [simpa['body'][i]['v'][0]*simpa['TimeStep'],	# 'body path'
									simpa['body'][i]['v'][1]*simpa['TimeStep']]	# 'body path'											
	#---------------------------------------------------------------------------
	return simpa																# return simulations parameter dictionary
################################################################################


import numpy as np

def read_sim(filename,simpa):
    data = np.fromfile(filename, dtype='f')
    data = np.reshape(data,simpa['GridSize'],order='F')
    axes = {}
    for ind, label in enumerate(['x','y','z','theta','t']):
        axes.update({label:np.arange(simpa['GridSize'][ind])*
                simpa['GridSpacing'][ind]})
    return data, axes
"""
def read_sec(filename,simpa):
    data = np.fromfile(filename, dtype='f')
    #data = np.reshape(data,simpa['GridSize'])
    return data

def read_trc(filename,simpa):
    data = np.fromfile(filename, dtype='f')
    #data = np.reshape(data,simpa['GridSize'])
    return data
"""


################################################################################
# 2. read_sim_data
################################################################################
def read_sim_data(simfilename,simpa,distances,components):
	'''
	----------------------------------------------------------------------------
	read Monte Carlo simulation data files
	----------------------------------------------------------------------------
	function read_sim_data reads the data ouput of the Monte Carlo 
	simulations and stores the data in a obspy.stream object for given 
	components and distances.
	----------------------------------------------------------------------------
	:param simfilename: 		name of the simulation 
	:param simpa:				name of the simpa dictionary
	:distances:					distances to create traces for
	:components:				components to create traces for
	----------------------------------------------------------------------------
	:return stream: 			obspy.stream object containing 
								simulation data for given components 
								and distances
	----------------------------------------------------------------------------
	'''
	import obspy.core															# import obspy.core package
	import array																# import array package
	import numpy as np															# import numpy package
	import glob																	# import glob package
	#---------------------------------------------------------------------------
	dmin_1 = max(int(simpa['GridSize'][0]*simpa['GridSpacing'][0]-				# calculate d_min part 1
					simpa['SourcePos'][0]),int(simpa['SourcePos'][0]))			# calculate d_min part 1
	dmin_2 = max(int(simpa['GridSize'][1]*simpa['GridSpacing'][1]-				# calculate d_min part 2
					simpa['SourcePos'][1]),int(simpa['SourcePos'][1]))			# calculate d_min part 2
	dmin   = max(dmin_1,dmin_2)													# calculate d_min part 3
	GridSpacing = (simpa['GridSpacing'][0]+simpa['GridSpacing'][1])/2			# calculate GridSpacing
	GridSize = dmin/GridSpacing													# calculate GridSize
	stream = obspy.core.Stream()												# create empty stream object
	count = 0																	# set counter
	for i in range(len(distances)):												# loop over all the distances
		for j in range(len(components)):										# loop over all the components
			file = glob.glob(simfilename+'*'+components[j]+'.sec')				# find according file
			data = np.fromfile(file[0], dtype='f')								# load data from file
			stream.append(obspy.core.Trace())									# append trace to stream object
			stream[count].stats.network		  = 'MC'							# set network name
			stream[count].stats.station		  = str(int(distances[i]))				# set station name to distance
			stream[count].stats.sampling_rate = 1/float(simpa['GridSpacing'][3])# set sampling rate
			stream[count].stats.channel		  = components[j][0:2]	 			# set channel entry
			stream[count].data =  data[int((distances[i]/GridSpacing)*			# load data for distance
											simpa['GridSize'][3]):
											int((distances[i]/GridSpacing)*
											simpa['GridSize'][3]+
											simpa['GridSize'][3])]
			count = count + 1													# raise count by 1
	#---------------------------------------------------------------------------
	return stream																# return stream
################################################################################



################################################################################
# 3. combine_sim_data_s_p_ratio
################################################################################
def combine_sim_data_s_p_ratio(stream_p,stream_s,s_p_ratio):
	'''
	----------------------------------------------------------------------------
	combines stream_p and stream_s to stream_ratio_s_p
	----------------------------------------------------------------------------
	function stream_ratio_ps combines the streams for P-Source and S-Source
	Monte Carlo Simulations in one new stream stream_ratio_s_p. S-energy to 
	P-energy ratio is given by Sato and Fehler, 1998
	------------------------------------------------------------------------
	:param stream_p: 			stream containing data for the p-pource simulations for 
								given components and distances
	:param stream_s: 			stream containing data for the s-source simulations for 
								given components and distances
	:param s_p_ratio:			S to P energy ratio
	------------------------------------------------------------------------
	:return stream_ratio_s_p:	stream containing combined P and S data
	------------------------------------------------------------------------
	'''	
	import obspy.core															# import obspy.core module
	#---------------------------------------------------------------------------
	stream_ratio_s_p = obspy.core.Stream()										# create empty stream object
	for i in range(len(stream_p)):												# loop over all traces in stream_p
		stream_ratio_s_p.append(obspy.core.Trace())								# append empty trace 
		stream_ratio_s_p[i].stats.network	   = stream_p[i].stats.network 		# set network name
		stream_ratio_s_p[i].stats.station	   = stream_p[i].stats.station      # set station name
		stream_ratio_s_p[i].stats.sampling_rate= stream_p[i].stats.sampling_rate# set sampling rate
		stream_ratio_s_p[i].stats.channel	   = stream_p[i].stats.channel      # set channel name
		stream_ratio_s_p[i].data 			   = (stream_p[i].data*(1-s_p_ratio)# create new data entry 
												 + stream_s[i].data*(s_p_ratio))# with P to S energy ratio
	#---------------------------------------------------------------------------
	return stream_ratio_s_p														# return the stream
################################################################################
