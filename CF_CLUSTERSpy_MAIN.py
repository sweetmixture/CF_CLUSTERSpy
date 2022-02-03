#!/usr/bin/python

from CF_CLUSTERSpy_MOD import * 

class cf_cluster(cf_cluster_mod):

	def scf_support_get_cost(self,etha):			# INPUT PARAMETERS ... CONFIG 1 / 2 , COM 1 / 2
		Return = 0.					# EXPANTION COEFFICIENT 'etha'
		r1 = 0.; r2 = 0.
		for i in range(self.c1_n):
			r1 = math.sqrt((self.c1_config[i][0]-self.c1_com[0])**2. + (self.c1_config[i][1]-self.c1_com[1])**2. + (self.c1_config[i][2]-self.c1_com[2])**2.)
			r2 = math.sqrt((self.c2_config[i][0]-self.c2_com[0])**2. + (self.c2_config[i][1]-self.c2_com[1])**2. + (self.c2_config[i][2]-self.c2_com[2])**2.)
			Return = Return + (r1-etha*r2)**2.
		return Return

	def scf_support_get_cost_deriv(self,etha):		# INPUT PARAMETERS ... CONFIG 1 / 2 , COM 1 / 2
		Return = 0.					# EXPANSION COEFFICIENT 'etha'
		for i in range(self.c1_n):
			r1 = math.sqrt((self.c1_config[i][0]-self.c1_com[0])**2. + (self.c1_config[i][1]-self.c1_com[1])**2. + (self.c1_config[i][2]-self.c1_com[2])**2.)
			r2 = math.sqrt((self.c2_config[i][0]-self.c2_com[0])**2. + (self.c2_config[i][1]-self.c2_com[1])**2. + (self.c2_config[i][2]-self.c2_com[2])**2.)
			Return = Return + (2.*etha*r2*r2 - 2.*r1*r2)
		return Return

	def get_scf_rms(self,etha):
		Return = 0.
		new_config_1 = self.c1_config.copy()
		new_config_2 = self.c2_config.copy()

		for i in range(self.c1_n):
			dx = new_config_1[i][0] - etha*new_config_2[i][0]
			dy = new_config_1[i][1] - etha*new_config_2[i][1]
			dz = new_config_1[i][2] - etha*new_config_2[i][2]
			Return = Return + dx**2. + dy**2. + dz**2.
		return math.sqrt(Return)/3./float(self.c1_n)


	def get_scf(self):					# HAS TO BE CALLED AFTER STRUCTURE MATCHING ( RECENTRE / ROTATION ) DONE
		self.etha = 1.					# INITIAL GUESSED EXPANSION COEFFICIENT
		alpha = 0.5;     beta = 0.96;	 t = 1.		# BACKTRACKING LINE SERACH WORKING VARIABLES
		trial_cost = 0.; dummy_step = 0.		#

		self.__SCF_ITER_MAX = 500				# scf max parameter
		self.__SCF_LIMIT_STEP_SIZE = 0.50			# scf optimisation (opti eta w.r.t cost function) step size limit

		self.__scf_out = []

		for n in range(self.__SCF_ITER_MAX):
			cost       = self.scf_support_get_cost(self.etha)		# CALCULATE CURRENT COST
			cost_deriv = self.scf_support_get_cost_deriv(self.etha)	# CALCULATE DERIVATIVE d(COST)/d(self.etha)

			scf_rms = self.get_scf_rms(self.etha)

			if n%4 == 0:						# LOG RECORDING
				self.__scf_out.append(" Cyc: %3d  Gnorm: %12.6f  RMS: %12.6f  Etha: %12.6f" % (n+1,math.fabs(cost_deriv),scf_rms,self.etha))

			if math.fabs(cost_deriv) < 10E-7 :			# DERIVATIVE TOLERANCE
				break
			t = 1.		
			if math.fabs(cost_deriv*t) > self.__SCF_LIMIT_STEP_SIZE :	# APPLY LIMITED STEP SIZE
				while True:
					t = t*beta
					if math.fabs(cost_deriv)*t < self.__SCF_LIMIT_STEP_SIZE :
						break
			while True:							# MAIN LOOP FOR OPTIMIZATION
				trial_cost = self.scf_support_get_cost( (self.etha - alpha*t*cost_deriv) )
				if trial_cost >= (cost - alpha*t*math.fabs(cost_deriv)) :
					t = t*beta
					if trial_cost < cost :
						self.etha = self.etha - alpha*t*cost_deriv
						break
					if t < 10E-10 :
						break
				else:
					self.etha = self.etha - alpha*t*cost_deriv
					break
		return self.etha	# THIS IS EXPANSION COEFFICIENT IN .... ( r1 - self.etha*r2 )

		'''
			for i in range(self.c1_n):
				out.write(" %2s%12.6f%12.6f%12.6f\n" % (self.c1_config[i][3],self.c1_config[i][0],self.c1_config[i][1],self.c1_config[i][2])) 
			for i in range(self.c2_n):
				out.write(" %2s%12.6f%12.6f%12.6f\n" % (self.c2_config[i][3],self.c2_config[i][0],self.c2_config[i][1],self.c2_config[i][2])) 
		'''


	def wrt_out(self):

		config_rms = self.get_rms()
		c1_pmoi    = [self.c1_eval_f.item(0),self.c1_eval_f.item(1),self.c1_eval_f.item(2)]
		c2_pmoi    = [self.c2_eval_f.item(0),self.c2_eval_f.item(1),self.c2_eval_f.item(2)]
		c1_pmoi.sort();	c2_pmoi.sort()

		pmoi_rms   = math.sqrt( (c1_pmoi[0]-c2_pmoi[0])**2. + (c1_pmoi[1]-c2_pmoi[1])**2. + (c1_pmoi[2]-c2_pmoi[2])**2. )
		pmoi_rms   = pmoi_rms/float(self.c1_n)/3.	

		with open("CF.out","w") as out:

			# INTRO COMMENT
			out.write("#"*60); out.write("\n")
			out.write("#"); out.write(" "*58); out.write("#"); out.write("\n")
			out.write("#"); out.write(" "*21); out.write("COMPARE CLUSTERS"); out.write(" "*21); out.write("#\n"); 
			out.write("#"); out.write(" "*58); out.write("#"); out.write("\n")
			out.write("#"); out.write(" "*21); out.write("Python Version 1"); out.write(" "*21); out.write("#\n"); 
			out.write("#"); out.write(" "*21); out.write("   23.01.2020   "); out.write(" "*21); out.write("#\n"); 
			out.write("#"); out.write(" "*58); out.write("#"); out.write("\n")
			out.write("#"); out.write(" "*21); out.write("Scott M Woodley "); out.write(" "*21); out.write("#\n"); 
			out.write("#"); out.write(" "*21); out.write("  Woongkyu Jee  "); out.write(" "*21); out.write("#\n"); 
			out.write("#"); out.write(" "*58); out.write("#"); out.write("\n")
			out.write("#"*60); out.write("\n")

			# OUTPUT COMMENT
			out.write("="*60); out.write("\n")
			out.write(" "*60); out.write("\n")
			out.write(" RESULTS FROM COMPARING CLUSTERS                            \n")
			out.write(" "*60); out.write("\n")
			out.write(" Using recentred / rotated structural data                  \n")
			out.write(" "*60); out.write("\n")
			out.write(" CLUSTER FILE INPUTS :\n")
			out.write(" "*60); out.write("\n")
			out.write(" CLUSTER 1 :  %s\n" % (self.xyz_file_1))
			out.write(" CLUSTER 2 :  %s\n" % (self.xyz_file_2))
			out.write(" "*60); out.write("\n")
			out.write("="*60); out.write("\n")
			out.write(" "*60); out.write("\n")

			out.write(" Degeneracy in the pricipal moment of inertia :  %d         \n" % (self.pmoi_deg))
			out.write(" "*60); out.write("\n")

			out.write(" Principal Moment of Inertia\n")
			out.write(" "*60); out.write("\n")

			out.write(" CLUSTER 1 : %12.6lf%12.6lf%12.6lf\n" % (c1_pmoi[0],c1_pmoi[1],c1_pmoi[2]))
			out.write(" CLUSTER 2 : %12.6lf%12.6lf%12.6lf\n" % (c2_pmoi[0],c2_pmoi[1],c2_pmoi[2]))
			out.write(" "*60); out.write("\n")
			out.write(" RMS PMOI  : %12.6f\n" % (pmoi_rms))			
			out.write(" "*60); out.write("\n")
			out.write("="*60); out.write("\n")
			out.write(" "*60); out.write("\n")

			out.write(" Final Cluster Configurations ( Recentred / Rotated ) \n")
			out.write(" "*60); out.write("\n")

			out.write(" CLUSTER 1 :\n")
			out.write(" "*60); out.write("\n")
			for i in range(self.c1_n):
				out.write(" %2s%12.6f%12.6f%12.6f\n" % (self.c1_config[i][3],self.c1_config[i][0],self.c1_config[i][1],self.c1_config[i][2])) 
			out.write(" "*60); out.write("\n")
			out.write(" CLUSTER 2 :\n")
			out.write(" "*60); out.write("\n")
			for i in range(self.c2_n):
				out.write(" %2s%12.6f%12.6f%12.6f\n" % (self.c2_config[i][3],self.c2_config[i][0],self.c2_config[i][1],self.c2_config[i][2])) 
			out.write(" "*60); out.write("\n")
			'''
			DEBUG
			'''
			if self.rms_config_org < config_rms:
				out.write(" O_RMS_CONFIG: %12.6f\t%12.10s\n" % (self.rms_config_org, "RAW"))
				out.write(" RMS CONFIG: %12.6f\t%12.10s\n" % (config_rms, "ROTATED"))
			else:
				out.write(" RMS CONFIG: %12.6f\t%12.10s\n" % (config_rms, "ROTATED"))
			#out.write(" RMS CONFIG: %12.6f\n" % (config_rms))			
			#out.write(" ORIGINAL_RMS_REFERENCE_ONLY: %12.6f\n" %(self.rms_config_org))
			out.write(" "*60); out.write("\n")
			out.write("="*60); out.write("\n")

			# RESCALING 
			out.write(" "*60); out.write("\n")
			out.write(" Scaling Factor\n")
			out.write(" "*60); out.write("\n")
			out.write(" Calculating optimised expansion coefficient (Etha)\n")
			out.write(" on 'cluster 2' with respect to 'cluster 1'\n")
			out.write(" "*60); out.write("\n")
			out.write(" ( Etha          )(x2)   (x1)\n")
			out.write(" (      Etha     )(y2) = (y1)  Find optimised 'Etha' value\n")
			out.write(" (          Etha )(z2)   (z1)\n")
			out.write(" "*60); out.write("\n")
			out.write(" Scaling Factor(Etha) Optimisation Log \n")
			out.write(" "*60); out.write("\n")
			for n in range(len(self.__scf_out)):
				out.write(" %s\n" % self.__scf_out[n])
			out.write(" "*60); out.write("\n")
			out.write(" Scaling Factor : %12.6f\n" % (self.etha))
			#out.write("                         COMING SOON !                      \n")
			out.write(" "*60); out.write("\n")
			out.write("="*60); out.write("\n")



			# END COMMENT
			out.write("#"*60); out.write("\n")
			out.write("#               COMPLETE SUCCESS - ALL DONE                #\n")
			out.write("#"*60); out.write("\n")
			out.write("#             THANK YOU FOR USING CF-CLUSTERS              #\n")
			out.write("#"*60); out.write("\n")

	def wrt_xyz(self):
		
		with open("CF_1.xyz","w") as f:
			f.write("\t%d\n" % (self.c1_n))
			f.write("%s\n" % (self.xyz_file_1))
			for i in range(self.c1_n):
				f.write("%2s%12.6f%12.6f%12.6f\n" % (self.c1_config[i][3],self.c1_config[i][0],self.c1_config[i][1],self.c1_config[i][2])) 

		with open("CF_2.xyz","w") as f:
			f.write("\t%d\n" % (self.c2_n))
			f.write("%s\n" % (self.xyz_file_2))
			for i in range(self.c2_n):
				f.write("%2s%12.6f%12.6f%12.6f\n" % (self.c2_config[i][3],self.c2_config[i][0],self.c2_config[i][1],self.c2_config[i][2])) 

		with open("CF_3.xyz","w") as f:
			f.write("\t%d\n" % (self.c1_n + self.c2_n))
			f.write("\n")
			for i in range(self.c1_n):
				f.write("%2s%12.6f%12.6f%12.6f\n" % (self.c1_config[i][3],self.c1_config[i][0],self.c1_config[i][1],self.c1_config[i][2])) 
			for i in range(self.c2_n):
				f.write("%2s%12.6f%12.6f%12.6f\n" % (self.c2_config[i][3],self.c2_config[i][0],self.c2_config[i][1],self.c2_config[i][2])) 
'''
	PLEASE READ COMMENTS BELOW
'''


if __name__ == '__main__':


	main = cf_cluster(sys.argv[1],sys.argv[2])		# INITIATE 'cf_cluster' INSTANCE INTO 'main'

	main.load_xyz()						# READ THE INPUT CLUSTER XYZ FILES
	main.get_com()						# CALCULATE CENTRE OF MASS
	main.get_shift()					# RECENTRE THE CLUSTERS WITH RESPECT TO CALCULATED CENTRE OF MASS
	main.get_moi()						# CALCULATE INERTIA TENSOR
	deg = main.get_deg()					# CHECK IF THERE IS DEGENERACY IN THE PRICIPAL MOMENTS

								# IF NO DEGENERACY - 'deg' = 1
								# TWO  - FOLD        'deg' = 2
								# THREE- FOLD        'deg' = 2

	main.get_eig(deg)					# BASED ON THE FOUND DEGENERACY, FIND THE MOST OPTIMAL TRANSFORMATION MATRICES 
	main.get_tra(deg)					# DO TRANSFORMATION
	main.get_com()						# RE-CALCULATE CENTRE OF MASS	( ONLY EFFECTIVE WHEN THERE IS DEGENERACY, SINCE IN THIS CASE THE ROTATION IS 
	main.get_shift()					# RE-CENTRE CLUSTERS		  CARRIED OUT WITH DUMMY ATOMS, WHICH ARE FOR REMOVING THE DEGENERACY )

	etha = main.get_scf()

	main.wrt_out()						# WRITE ANALYSIS OUTPUT (WILL BE WRITTEN IN 'CF.out')
	main.wrt_xyz()						# WRITE RECENTRED / ROTATED CLUSTERS (WILL BE WRITTEN IN 'CF_1.xyz' and 'CF_2.xyz' )


	# IF YOU WANT TO RUN MANY FILES AT THE SAME TIME, LOOk .cmp.sh shell script in this directory
