import sys, os, string
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table

class synClean:

	def __init__(self):

			self.rootdir = os.getcwd()+'/'


	def writeSynageTable(self,freq,fluxIntegrated,fluxErr,title):

			synTable = self.rootdir+'synageFl'+title+'.tbl'

			title1='point 1\nlabel: '+title

			f = open(synTable, 'w')
			f.write('freq_units: MHz\nflux_units: Jy\n\n')

			f.write(title1+'\n')

			for i in xrange(0,len(freq)):
					line = 'n '+str(round(freq[i]/1e6,1))+' '+str(fluxIntegrated[i])+' '+str(fluxErr[i])+'\n'        
		
					f.write(line)

			f.close()
			

			return synTable

	def readCIStats(self,CI_table):
			
			table = ascii.read(CI_table)

			mod_chi = float(table['col3'])
			mod_ndf = float(table['col4'])
			mod_chired = float(table['col5'])
			br_f = float(table['col6'])
			nuinf_err = float(table['col7'])
			nusup_err = float(table['col8'])
			norm = float(table['col9'])
			norminf_err = float(table['col10'])
			normsup_err = float(table['col11'])
			alpha_inj = float(table['col12'])
			alphainf_err = float(table['col13'])
			alphasup_err = float(table['col14'])

			stats={'alpha':alpha_inj,
						 'alpha_errinf':alphainf_err,
						 'alpha_errsup':alphasup_err,
						 'break':br_f,
						 'break_sup': nusup_err,
						 'break_inf': nuinf_err,
						 'norm':norm,
						 'chisq':mod_chi,
						 'ndf':mod_ndf,
						 'chired':mod_chired}         		
			
			return stats

	def readCIMod(self, CI_table):
		
		f= open(CI_table)

		freq_mod = []
		flux_mod = []
		i=0
		for line in f:

				if not line.strip():    
						continue
				line = string.strip(line)
				if line[0] == '#':
						i+=1
						continue

				line = string.split(line,'\t')
				freq_mod.append(line[0])
				flux_mod.append(line[1])

				
		freq_mod = np.array(freq_mod[:-12],dtype=float)
		flux_mod = np.array(flux_mod[:-12],dtype=float)
		f.close()
		array=np.array([freq_mod,flux_mod])
		
		return array

	def readCIOFFStats(self,CIOFF_table):
			
			table = ascii.read(CIOFF_table)

			mod_chi = float(table['col3'])
			mod_ndf = float(table['col4'])
			mod_chired = float(table['col5'])
			br_f = float(table['col6'])
			nuinf_err = float(table['col7'])
			nusup_err = float(table['col8'])
			norm = float(table['col9'])
			norminf_err = float(table['col10'])
			normsup_err = float(table['col11'])
			alpha_inj = float(table['col12'])
			alphainf_err = float(table['col13'])
			alphasup_err = float(table['col14'])
			t_off = float(table['col15'])
			tinf_err = float(table['col16'])
			tsup_err = float(table['col17'])

			stats={'alpha':alpha_inj,
						 'alpha_errinf':alphainf_err,
						 'alpha_errsup':alphasup_err,
						 'break':br_f,
						 'break_sup': nusup_err,
						 'break_inf': nuinf_err,
						 'norm':norm,
						 'chisq':mod_chi,
						 'ndf':mod_ndf,
						 'chired':mod_chired,
						 'tratio':t_off,          
						 't_errinf':tinf_err,
						 't_errsup':tsup_err}         		
			
			return stats

	def readCIOFFMod(self, CIOFF_table):

		f= open(CIOFF_table)

		freq_mod = []
		flux_mod = []
		i=0
		for line in f:

				if not line.strip():    
						i+=1
						continue
				line = string.strip(line)
				if line[0] == '#':
						i+=1
						continue

				line = string.split(line,'\t')
				freq_mod.append(line[0])
				flux_mod.append(line[1])
				i+=1

		freq_mod = np.array(freq_mod[:-12],dtype=float)
		flux_mod = np.array(flux_mod[:-12],dtype=float)
		f.close()
		array=np.array([freq_mod,flux_mod])
		
		return array
