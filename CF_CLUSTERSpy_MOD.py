#!/usr/bin/bash

import sys,math
import numpy as np
from random import seed, randint

kd = lambda i,j : 1. if (i==j) else 0.

class cf_cluster_mod:

	def __init__(self,*f):

		self.xyz_file_1 = f[0]
		self.xyz_file_2 = f[1]

		self.__DEGEN_TOL     = 0.025
		self.__DEGEN_DL      = 7.5

	def load_xyz(self):

		try:
			with open(self.xyz_file_1,"r") as f1:

				self.c1_n = f1.readline()
				self.c1_n = int(self.c1_n)
				rws       = f1.readline()
				self.c1_config = [[0 for i in range(4)] for j in range(self.c1_n+2)]					# HERE +2 IS FOR DUMMY ATOM WS

				for i in range(self.c1_n):
					rws = f1.readline()
					spl = rws.split()
					self.c1_config[i][0] = float(spl[1])
					self.c1_config[i][1] = float(spl[2])
					self.c1_config[i][2] = float(spl[3])
					self.c1_config[i][3] = spl[0]
		
		except FileNotFoundError:
			print("Error, input file 1 not found ...")
			sys.exit()

	
		try: 
			with open(self.xyz_file_2,"r") as f2:

				self.c2_n = f2.readline()
				self.c2_n = int(self.c2_n)
				rws       = f2.readline()
				self.c2_config = [[0 for i in range(4)] for j in range(self.c2_n+2)]					# HERE +2 IS FOR DUMMY ATOM WS

				for i in range(self.c1_n):
					rws = f2.readline()
					spl = rws.split()
					self.c2_config[i][0] = float(spl[1])
					self.c2_config[i][1] = float(spl[2])
					self.c2_config[i][2] = float(spl[3])
					self.c2_config[i][3] = spl[0]
		except FileNotFoundError:
			print("Error, input file 2 not found ...")
			sys.exit()



		if self.c1_n != self.c2_n:
			print("ERROR: number of atoms doesn't match")
			sys.exit()

		return None


	def get_com(self,mode=1):
		'''
			DESCRIPTION:
			Calculate the centre of mass of cluster 1 and 2,
			The results will be saved into the internal variables 'c1_com / c2_com'
			If there is no degeneracy in the pricipal moments of clusters, use the method with mode 1,
			or if there is 2/3 fold degeneracy then use with mode 2/3
		'''

		self.c1_com = [ 0. for i in range(3) ]
		self.c2_com = [ 0. for j in range(3) ]

		if   mode == 1:			# NO DEG
			an = self.c1_n + 0
		elif mode == 2:			# 2-FOLD DEG
			an = self.c1_n + 1
		elif mode == 3:			# 3-FOLD DEG
			an = self.c1_n + 2	# THE FINAL 'an' VARAIBLE SAVES THE NUMBER OF ATOMS IN THE SYSTEM (+ DUMMY ATOMS FOR GETTING RID OF DEGENERACY)

		for i in range(an):
			for j in range(3):
				self.c1_com[j] = self.c1_com[j] + self.c1_config[i][j]
				self.c2_com[j] = self.c2_com[j] + self.c2_config[i][j]

		self.c1_com[0] = self.c1_com[0]/float(an);	self.c2_com[0] = self.c2_com[0]/float(an)
		self.c1_com[1] = self.c1_com[1]/float(an);	self.c2_com[1] = self.c2_com[1]/float(an)
		self.c1_com[2] = self.c1_com[2]/float(an);	self.c2_com[2] = self.c2_com[2]/float(an)

		return None

	def get_shift(self,mode=1):
		'''
			DESCRIPTION:
			This method has to be called after 'get_com()' achieved.
			Based on the calculated centre of mass, the method recentre (translate) the com to the origin.

			After the shift, you might need to recalculate the centre of mass again (see 'get_com()').
		'''
		if   mode == 1:			# NO DEG
			an = self.c1_n + 0
		elif mode == 2:			# 2-FOLD DEG
			an = self.c1_n + 1
		elif mode == 3:			# 3-FOLD DEG
			an = self.c1_n + 2

		for i in range(an):
			for j in range(3):
				self.c1_config[i][j] = self.c1_config[i][j] - self.c1_com[j]
				self.c2_config[i][j] = self.c2_config[i][j] - self.c2_com[j]
		return None


	def get_rms(self):
		'''
			DESCRIPTION:
			This method calculates RMS of the coordinates of the two different clusters.
			It returns the RMS value directly.

			Note that this method does not reflect the dummy atoms when it calculates RMS.
		'''
		Return = 0.
		for i in range(self.c1_n):
			dx = self.c1_config[i][0] - self.c2_config[i][0]
			dy = self.c1_config[i][1] - self.c2_config[i][1]
			dz = self.c1_config[i][2] - self.c2_config[i][2]
			Return = Return + dx**2. + dy**2. + dz**2.
		return math.sqrt(Return)/3./float(self.c1_n)


	def get_moi(self,mode=1):
		'''
			DESCRIPTION:
			This method calculates inertia tensor w.r.t the given cluster configuration.
			You might need to carry out 'get_com() -> get_shift()' on the cluster to find uniquely defined
			principal moment of inertia.
		'''

		if   mode == 1:			# NO DEG
			an = self.c1_n + 0
		elif mode == 2:			# 2-FOLD DEG
			an = self.c1_n + 1
		elif mode == 3:			# 3-FOLD DEG
			an = self.c1_n + 2

		# DO NOT TRUST THE MOI TENSOR WHEN THIS METHOD IS USED WITH MODE EITHER 2 OR 3 ... WHEN IT INCLUDES DUMMY ATOMS

		self.c1_moi_np = np.matrix([[0. for i in range(3)] for j in range(3)])
		self.c2_moi_np = np.matrix([[0. for i in range(3)] for j in range(3)])

		for i in range(an):
			c1_r2 = self.c1_config[i][0]**2. + self.c1_config[i][1]**2. + self.c1_config[i][2]**2.
			c2_r2 = self.c2_config[i][0]**2. + self.c2_config[i][1]**2. + self.c2_config[i][2]**2.

			for j in range(3):
				for k in range(3):
					self.c1_moi_np.itemset((j,k),self.c1_moi_np.item((j,k)) \
						+ (c1_r2*kd(j,k) - self.c1_config[i][j]*self.c1_config[i][k]))
					self.c2_moi_np.itemset((j,k),self.c2_moi_np.item((j,k)) \
						+ (c2_r2*kd(j,k) - self.c2_config[i][j]*self.c2_config[i][k]))
		return None

	def get_deg(self):
		'''
			DESCROIPTION:
		'''

		Return = 0
		self.__DEGEN_TOL_MIN = 1. - self.__DEGEN_TOL					# IF A PAIR OF EVALS IN 0.95 SIMILARITY -> DEGENERACY RANK GOES UP
		self.__DEGEN_TOL_MAX = 1. + self.__DEGEN_TOL
		c2_deg_eval = None;	c2_deg_evec = None;					# THIS METHOD TAKES ONLY THE SECOND CLUSTER TO SEE IF THERE IS DEGENERACY IN IT'S PRICIPAL MOMENTS
		flag = False
		c2_deg_eval, c2_deg_evec = np.linalg.eig(self.c2_moi_np)

		if self.__DEGEN_TOL_MIN < c2_deg_eval.item((1))/c2_deg_eval.item((0)) and c2_deg_eval.item((1))/c2_deg_eval.item((0)) < self.__DEGEN_TOL_MAX:
			Return = Return + 1
			flag = True
		if self.__DEGEN_TOL_MIN < c2_deg_eval.item((2))/c2_deg_eval.item((0)) and c2_deg_eval.item((2))/c2_deg_eval.item((0)) < self.__DEGEN_TOL_MAX:
			Return = Return + 1
			if flag == True:
				self.pmoi_deg = Return + 1
				return Return + 1
		if self.__DEGEN_TOL_MIN < c2_deg_eval.item((2))/c2_deg_eval.item((1)) and c2_deg_eval.item((2))/c2_deg_eval.item((1)) < self.__DEGEN_TOL_MAX:
			Return = Return + 1

		self.pmoi_deg = Return + 1
		return Return + 1

	def set_dummy_atom(self,mode):
		'''
			DESCROIPTION:
			Cluster recentred is necessary before taking this step
		'''
		seed(1)
		if   mode == 2:											# 2-FOLD DEG
			an = self.c1_n + 1									# NEED EXTRA 1 DUMMY ATOM TO GET RID OF DEGENERACY
			rnd_id = randint(0,self.c1_n)								# PICK ONE
			c1_v_mag = math.sqrt( self.c1_config[rnd_id][0]**2. + self.c1_config[rnd_id][1]**2. + self.c1_config[rnd_id][2]**2. )
			c2_v_mag = math.sqrt( self.c2_config[rnd_id][0]**2. + self.c2_config[rnd_id][1]**2. + self.c2_config[rnd_id][2]**2. )
			c1_dv = [ self.c1_config[rnd_id][0]/c1_v_mag, self.c1_config[rnd_id][1]/c1_v_mag, self.c1_config[rnd_id][2]/c1_v_mag ]
			c2_dv = [ self.c2_config[rnd_id][0]/c2_v_mag, self.c2_config[rnd_id][1]/c2_v_mag, self.c2_config[rnd_id][2]/c2_v_mag ]
			for i in range(3):
				c1_dv[i] = self.__DEGEN_DL*c1_dv[i];	c2_dv[i] = self.__DEGEN_DL*c2_dv[i]	# GET DUMMY ATOM POS
				self.c1_config[an-1][i] = c1_dv[i];	self.c2_config[an-1][i] = c2_dv[i]	# SET DUMMY ATOM IN CONFIG

		elif mode == 3:											# 3-FOLD DEG
			an = self.c1_n + 2									# NEED EXTRA 2 DUMMY ATOMS TO GET RID OF DEGENERACY

			rnd_id_a = randint(0,self.c1_n)
			while True:
				rnd_id_b = randint(0,self.c1_n)
				if rnd_id_a != rnd_id_b:							# PICK TWO DIFFERENT ATOMS
					break

			c1_v_mag_a = math.sqrt( self.c1_config[rnd_id_a][0]**2. + self.c1_config[rnd_id_a][1]**2. + self.c1_config[rnd_id_a][2]**2. )
			c1_v_mag_b = math.sqrt( self.c1_config[rnd_id_b][0]**2. + self.c1_config[rnd_id_b][1]**2. + self.c1_config[rnd_id_b][2]**2. )
			c2_v_mag_a = math.sqrt( self.c2_config[rnd_id_a][0]**2. + self.c2_config[rnd_id_a][1]**2. + self.c2_config[rnd_id_a][2]**2. )
			c2_v_mag_b = math.sqrt( self.c2_config[rnd_id_b][0]**2. + self.c2_config[rnd_id_b][1]**2. + self.c2_config[rnd_id_b][2]**2. )
			c1_dv_a = [ self.c1_config[rnd_id_a][0]/c1_v_mag_a, self.c1_config[rnd_id_a][1]/c1_v_mag_a, self.c1_config[rnd_id_a][2]/c1_v_mag_a ]
			c1_dv_b = [ self.c1_config[rnd_id_b][0]/c1_v_mag_b, self.c1_config[rnd_id_b][1]/c1_v_mag_b, self.c1_config[rnd_id_b][2]/c1_v_mag_b ]
			c2_dv_a = [ self.c2_config[rnd_id_a][0]/c2_v_mag_a, self.c2_config[rnd_id_a][1]/c2_v_mag_a, self.c2_config[rnd_id_a][2]/c2_v_mag_a ]
			c2_dv_b = [ self.c2_config[rnd_id_b][0]/c2_v_mag_b, self.c2_config[rnd_id_b][1]/c2_v_mag_b, self.c2_config[rnd_id_b][2]/c2_v_mag_b ]

			for i in range(3):
				c1_dv_a[i] = self.__DEGEN_DL*c1_dv_a[i];	c1_dv_b[i] = self.__DEGEN_DL*c1_dv_b[i]
				c2_dv_a[i] = self.__DEGEN_DL*c2_dv_a[i];	c2_dv_b[i] = self.__DEGEN_DL*c2_dv_b[i]
				self.c1_config[an-2][i] = c1_dv_a[i];	self.c1_config[an-1][i] = c1_dv_b[i]
				self.c2_config[an-2][i] = c2_dv_a[i];	self.c2_config[an-1][i] = c2_dv_b[i]


	def eig_support_sort(self,val,vec):
		'''
			DESCROIPTION:
		'''
		vec_tmp = [[0.for i in range(3)] for j in range(3)]
		for i in range(3):
			for j in range(3):
				vec_tmp[i][j] = vec.item((i,j))
		val_dic = {1:val.item(0),2:val.item(1),3:val.item(2)}
		val_dic = {k: v for k, v in sorted(val_dic.items(), key=lambda item: item[1])}
		key_idx = list(val_dic.keys())
		val = np.sort(val)
		for i in range(3):
			for j in range(3):
				vec.itemset((j,i),vec_tmp[j][key_idx[i]-1])
		return val, vec

	def eig_support_rms(self,config1,config2):
		'''
			DESCROIPTION:
		'''
		Return = 0.
		for i in range(self.c1_n):
			dx = config1[i][0] - config2[i][0]
			dy = config1[i][1] - config2[i][1]
			dz = config1[i][2] - config2[i][2]
			Return = Return + dx**2. + dy**2. + dz**2.
		return math.sqrt(Return)/3./float(self.c1_n)

	def eig_support_match(self, config1, vec1, config2, vec2, noa):
		'''
			DESCROIPTION:
		'''

		cf1 = [[0.for i in range(3)] for j in range(noa)]		# DUMMY CONFIG WORKSPACE
		cf2 = [[0.for i in range(3)] for j in range(noa)]
		for i in range(noa):
			for j in range(3):
				cf1[i][j] = config1[i][j]; cf2[i][j] = config2[i][j]
		# 8 DIFFERENT WAYS
		rms = {0:0.,1:0.,2:0.,3:0.,4:0.,5:0.,6:0.,7:0.}
		# TRANSFORM REFERENCE CLUSTER - 2
		c2_trans_mat = vec2.transpose()
		for i in range(noa):
			vec_tmp = np.array([config2[i][0],config2[i][1],config2[i][2]])
			vec_res = np.dot(c2_trans_mat,vec_tmp)
			cf2[i][0] = vec_res.item(0);	cf2[i][1] = vec_res.item(1);	cf2[i][2] = vec_res.item(2)
		# TRANSFORM REFERENCE CLUSTER - 2 DONE

		for i in range(8):

			if i == 0 :	# DEFAULT
				c1_trans_mat = np.copy(vec1)
				c1_trans_mat = c1_trans_mat.transpose()
				for n in range(noa):
					vec_tmp = np.array([config1[n][0],config1[n][1],config1[n][2]])
					vec_res = np.dot(c1_trans_mat,vec_tmp)
					cf1[n][0] = vec_res.item(0);	cf1[n][1] = vec_res.item(1);	cf1[n][2] = vec_res.item(2)
				rms[i] = self.eig_support_rms(cf1,cf2)
			if i == 1 :	# 1
				c1_trans_mat = np.copy(vec1)
				c1_trans_mat.itemset((0,0),-1.*c1_trans_mat.item(0,0))
				c1_trans_mat.itemset((1,0),-1.*c1_trans_mat.item(1,0))
				c1_trans_mat.itemset((2,0),-1.*c1_trans_mat.item(2,0))

				c1_trans_mat = c1_trans_mat.transpose()
				for n in range(noa):
					vec_tmp = np.array([config1[n][0],config1[n][1],config1[n][2]])
					vec_res = np.dot(c1_trans_mat,vec_tmp)
					cf1[n][0] = vec_res.item(0);	cf1[n][1] = vec_res.item(1);	cf1[n][2] = vec_res.item(2)
				rms[i] = self.eig_support_rms(cf1,cf2)
			if i == 2 :	# 2
				c1_trans_mat = np.copy(vec1)
				c1_trans_mat.itemset((0,1),-1.*c1_trans_mat.item(0,1))
				c1_trans_mat.itemset((1,1),-1.*c1_trans_mat.item(1,1))
				c1_trans_mat.itemset((2,1),-1.*c1_trans_mat.item(2,1))

				c1_trans_mat = c1_trans_mat.transpose()
				for n in range(noa):
					vec_tmp = np.array([config1[n][0],config1[n][1],config1[n][2]])
					vec_res = np.dot(c1_trans_mat,vec_tmp)
					cf1[n][0] = vec_res.item(0);	cf1[n][1] = vec_res.item(1);	cf1[n][2] = vec_res.item(2)
				rms[i] = self.eig_support_rms(cf1,cf2)
			if i == 3 :	# 3
				c1_trans_mat = np.copy(vec1)
				c1_trans_mat.itemset((0,2),-1.*c1_trans_mat.item(0,2))
				c1_trans_mat.itemset((1,2),-1.*c1_trans_mat.item(1,2))
				c1_trans_mat.itemset((2,2),-1.*c1_trans_mat.item(2,2))

				c1_trans_mat = c1_trans_mat.transpose()
				for n in range(noa):
					vec_tmp = np.array([config1[n][0],config1[n][1],config1[n][2]])
					vec_res = np.dot(c1_trans_mat,vec_tmp)
					cf1[n][0] = vec_res.item(0);	cf1[n][1] = vec_res.item(1);	cf1[n][2] = vec_res.item(2)
				rms[i] = self.eig_support_rms(cf1,cf2)
			if i == 4 :	# 1/2
				c1_trans_mat = np.copy(vec1)
				c1_trans_mat.itemset((0,0),-1.*c1_trans_mat.item(0,0))
				c1_trans_mat.itemset((1,0),-1.*c1_trans_mat.item(1,0))
				c1_trans_mat.itemset((2,0),-1.*c1_trans_mat.item(2,0))
				c1_trans_mat.itemset((0,1),-1.*c1_trans_mat.item(0,1))
				c1_trans_mat.itemset((1,1),-1.*c1_trans_mat.item(1,1))
				c1_trans_mat.itemset((2,1),-1.*c1_trans_mat.item(2,1))

				c1_trans_mat = c1_trans_mat.transpose()
				for n in range(noa):
					vec_tmp = np.array([config1[n][0],config1[n][1],config1[n][2]])
					vec_res = np.dot(c1_trans_mat,vec_tmp)
					cf1[n][0] = vec_res.item(0);	cf1[n][1] = vec_res.item(1);	cf1[n][2] = vec_res.item(2)
				rms[i] = self.eig_support_rms(cf1,cf2)
			if i == 5 :	# 1/3
				c1_trans_mat = np.copy(vec1)
				c1_trans_mat.itemset((0,0),-1.*c1_trans_mat.item(0,0))
				c1_trans_mat.itemset((1,0),-1.*c1_trans_mat.item(1,0))
				c1_trans_mat.itemset((2,0),-1.*c1_trans_mat.item(2,0))
				c1_trans_mat.itemset((0,2),-1.*c1_trans_mat.item(0,2))
				c1_trans_mat.itemset((1,2),-1.*c1_trans_mat.item(1,2))
				c1_trans_mat.itemset((2,2),-1.*c1_trans_mat.item(2,2))

				c1_trans_mat = c1_trans_mat.transpose()
				for n in range(noa):
					vec_tmp = np.array([config1[n][0],config1[n][1],config1[n][2]])
					vec_res = np.dot(c1_trans_mat,vec_tmp)
					cf1[n][0] = vec_res.item(0);	cf1[n][1] = vec_res.item(1);	cf1[n][2] = vec_res.item(2)
				rms[i] = self.eig_support_rms(cf1,cf2)
			if i == 6 :	# 2/3
				c1_trans_mat = np.copy(vec1)
				c1_trans_mat.itemset((0,1),-1.*c1_trans_mat.item(0,1))
				c1_trans_mat.itemset((1,1),-1.*c1_trans_mat.item(1,1))
				c1_trans_mat.itemset((2,1),-1.*c1_trans_mat.item(2,1))
				c1_trans_mat.itemset((0,2),-1.*c1_trans_mat.item(0,2))
				c1_trans_mat.itemset((1,2),-1.*c1_trans_mat.item(1,2))
				c1_trans_mat.itemset((2,2),-1.*c1_trans_mat.item(2,2))

				c1_trans_mat = c1_trans_mat.transpose()
				for n in range(noa):
					vec_tmp = np.array([config1[n][0],config1[n][1],config1[n][2]])
					vec_res = np.dot(c1_trans_mat,vec_tmp)
					cf1[n][0] = vec_res.item(0);	cf1[n][1] = vec_res.item(1);	cf1[n][2] = vec_res.item(2)
				rms[i] = self.eig_support_rms(cf1,cf2)
			if i == 7 :	# 2/3
				c1_trans_mat = np.copy(vec1)
				c1_trans_mat.itemset((0,0),-1.*c1_trans_mat.item(0,0))
				c1_trans_mat.itemset((1,0),-1.*c1_trans_mat.item(1,0))
				c1_trans_mat.itemset((2,0),-1.*c1_trans_mat.item(2,0))
				c1_trans_mat.itemset((0,1),-1.*c1_trans_mat.item(0,1))
				c1_trans_mat.itemset((1,1),-1.*c1_trans_mat.item(1,1))
				c1_trans_mat.itemset((2,1),-1.*c1_trans_mat.item(2,1))
				c1_trans_mat.itemset((0,2),-1.*c1_trans_mat.item(0,2))
				c1_trans_mat.itemset((1,2),-1.*c1_trans_mat.item(1,2))
				c1_trans_mat.itemset((2,2),-1.*c1_trans_mat.item(2,2))

				c1_trans_mat = c1_trans_mat.transpose()
				for n in range(noa):
					vec_tmp = np.array([config1[n][0],config1[n][1],config1[n][2]])
					vec_res = np.dot(c1_trans_mat,vec_tmp)
					cf1[n][0] = vec_res.item(0);	cf1[n][1] = vec_res.item(1);	cf1[n][2] = vec_res.item(2)
				rms[i] = self.eig_support_rms(cf1,cf2)

		rms = {k: v for k, v in sorted(rms.items(), key=lambda item: item[1])}
		key_idx = list(rms.keys())

		if key_idx[0] == 0:
			a = +1.;	b = +1.;	c = +1.;	# Default
		if key_idx[0] == 1:
			a = -1.;	b = +1.;	c = +1.;	# 1
		if key_idx[0] == 2:
			a = +1.;	b = -1.;	c = +1.;	# 2
		if key_idx[0] == 3:
			a = +1.;	b = +1.;	c = -1.;	# 3
		if key_idx[0] == 4:
			a = -1.;	b = -1.;	c = +1.;	# 1/2
		if key_idx[0] == 5:
			a = -1.;	b = +1.;	c = -1.;	# 1/3
		if key_idx[0] == 6:
			a = +1.;	b = -1.;	c = -1.;	# 2/3
		if key_idx[0] == 7:
			a = -1.;	b = -1.;	c = -1.;	# 1/2/3

		Return = np.copy(vec1)
		Return.itemset((0,0),a*Return.item(0,0))
		Return.itemset((1,0),a*Return.item(1,0))
		Return.itemset((2,0),a*Return.item(2,0))
		Return.itemset((0,1),b*Return.item(0,1))
		Return.itemset((1,1),b*Return.item(1,1))
		Return.itemset((2,1),b*Return.item(2,1))
		Return.itemset((0,2),c*Return.item(0,2))
		Return.itemset((1,2),c*Return.item(1,2))
		Return.itemset((2,2),c*Return.item(2,2))

		return Return

	def get_eig(self,mode=1):
		'''
			DESCROIPTION:
		'''
		# IN THE FIRST STEP, THE METHOD GETS EIG SYS OF THE ORIGINAL INPUTS
		self.c1_eval, self.c1_evec = np.linalg.eig(self.c1_moi_np)
		self.c2_eval, self.c2_evec = np.linalg.eig(self.c2_moi_np)

		self.c1_eval_f = np.copy(self.c1_eval)										# FOR LATTER USE
		self.c2_eval_f = np.copy(self.c2_eval)										# SAVE ORIGINAL EVALS

		if   mode == 1:			# NO DEG
			an = self.c1_n + 0
		elif mode == 2:			# 2-FOLD DEG
			an = self.c1_n + 1
		elif mode == 3:			# 3-FOLD DEG
			an = self.c1_n + 2

		# CASE 1 - IF THERE IS NO DEGENERACY IN CLUSTER 2 PRINCIPAL MOMENT OF INERTIA
		if mode == 1:
			# EIG SYS ALREADY DONE	# SORT EIGSYS
			self.c1_eval, self.c1_evec = self.eig_support_sort(self.c1_eval,self.c1_evec)
			self.c2_eval, self.c2_evec = self.eig_support_sort(self.c2_eval,self.c2_evec)
			self.c1_evec = self.eig_support_match(self.c1_config,self.c1_evec,self.c2_config,self.c2_evec,an)	# AFTER THIS POINT, CORRECT EIGENVECTOR FOR TRANSFORMATION ACHIEVED

		if mode == 2:	# FOR 2 FOLD DEG
			self.c1_eval_d = None;	self.c2_eval_d = None;								# EMPTY EIGSYS WORKSPACE FOR WORKING WITH DUMMY ATOMS
			self.c1_evec_d = None;	self.c2_evec_d = None;

			self.set_dummy_atom(mode)										# SET DUMMY ATOM
			self.get_com(mode)											# RECALCULATE COM WITH DUMMY ATOM
			self.get_shift(mode)											# RECENTRE COM WITH DUMMY ATOM
			self.get_moi(mode)											# RECALCULATE MOI WITH DUMMY ATOM / IN SHIFTED CENTRE

			self.c1_eval_d, self.c1_evec_d = np.linalg.eig(self.c1_moi_np)						# GET NEW EIGSYS
			self.c2_eval_d, self.c2_evec_d = np.linalg.eig(self.c2_moi_np)

			self.c1_eval_d, self_c1_evec_d = self.eig_support_sort(self.c1_eval_d,self.c1_evec_d)			# SORT THE EIGENVALUES BY MAGNITUDES
			self.c2_eval_d, self_c2_evec_d = self.eig_support_sort(self.c2_eval_d,self.c2_evec_d)

			self.c1_evec_d = self.eig_support_match(self.c1_config,self.c1_evec_d,self.c2_config,self.c2_evec_d,an)	# AFTER THIS POINT, CORRECT EIGENVECTOR FOR THE DUMMY TRANSFORMATION ACHIEVED

		if mode == 3:	# FOR 3 FOLD DEG
			self.c1_eval_d = None;	self.c2_eval_d = None;								# EMPTY EIGSYS WORKSPACE FOR WORKING WITH DUMMY ATOMS
			self.c1_evec_d = None;	self.c2_evec_d = None;

			self.set_dummy_atom(mode)										# SET DUMMY ATOMS
			self.get_com(mode)											# RECALCULATE COM WITH DUMMY ATOMS
			self.get_shift(mode)											# RECENTRE COM WITH DUMMY ATOMS
			self.get_moi(mode)											# RECALCULATE MOI WITH DUMMY ATOMS / IN SHIFTED CENTRE

			self.c1_eval_d, self.c1_evec_d = np.linalg.eig(self.c1_moi_np)						# GET NEW EIGSYS
			self.c2_eval_d, self.c2_evec_d = np.linalg.eig(self.c2_moi_np)

			self.c1_eval_d, self_c1_evec_d = self.eig_support_sort(self.c1_eval_d,self.c1_evec_d)			# SORT THE EIGENVALUES BY MAGNITUDES
			self.c2_eval_d, self_c2_evec_d = self.eig_support_sort(self.c2_eval_d,self.c2_evec_d)

			self.c1_evec_d = self.eig_support_match(self.c1_config,self.c1_evec_d,self.c2_config,self.c2_evec_d,an)	# AFTER THIS POINT, CORRECT EIGENVECTOR FOR THE DUMMY TRANSFORMATION ACHIEVED

	def get_tra(self,mode=1):
		'''
			DESCROIPTION:
		'''
		if   mode == 1:			# NO DEG
			an = self.c1_n + 0

			c1_trans = np.copy(self.c1_evec);	c2_trans = np.copy(self.c2_evec)
			c1_trans = np.transpose(c1_trans);	c2_trans = np.transpose(c2_trans)

			for i in range(an):
				v_1 = np.array([self.c1_config[i][0],self.c1_config[i][1],self.c1_config[i][2]])
				v_2 = np.array([self.c2_config[i][0],self.c2_config[i][1],self.c2_config[i][2]])
				v_res_1 = np.dot(c1_trans,v_1)
				v_res_2 = np.dot(c2_trans,v_2)

				self.c1_config[i][0] = v_res_1.item(0); self.c1_config[i][1] = v_res_1.item(1); self.c1_config[i][2] = v_res_1.item(2)
				self.c2_config[i][0] = v_res_2.item(0); self.c2_config[i][1] = v_res_2.item(1); self.c2_config[i][2] = v_res_2.item(2)

		elif mode == 2:			# 2-FOLD DEG
			an = self.c1_n + 1

			c1_trans = np.copy(self.c1_evec_d);	c2_trans = np.copy(self.c2_evec_d)
			c1_trans = np.transpose(c1_trans);	c2_trans = np.transpose(c2_trans)

			for i in range(an):
				v_1 = np.array([self.c1_config[i][0],self.c1_config[i][1],self.c1_config[i][2]])
				v_2 = np.array([self.c2_config[i][0],self.c2_config[i][1],self.c2_config[i][2]])
				v_res_1 = np.dot(c1_trans,v_1)
				v_res_2 = np.dot(c2_trans,v_2)

				self.c1_config[i][0] = v_res_1.item(0); self.c1_config[i][1] = v_res_1.item(1); self.c1_config[i][2] = v_res_1.item(2)
				self.c2_config[i][0] = v_res_2.item(0); self.c2_config[i][1] = v_res_2.item(1); self.c2_config[i][2] = v_res_2.item(2)

		elif mode == 3:			# 3-FOLD DEG
			an = self.c1_n + 2

			c1_trans = np.copy(self.c1_evec_d);	c2_trans = np.copy(self.c2_evec_d)
			c1_trans = np.transpose(c1_trans);	c2_trans = np.transpose(c2_trans)

			for i in range(an):
				v_1 = np.array([self.c1_config[i][0],self.c1_config[i][1],self.c1_config[i][2]])
				v_2 = np.array([self.c2_config[i][0],self.c2_config[i][1],self.c2_config[i][2]])
				v_res_1 = np.dot(c1_trans,v_1)
				v_res_2 = np.dot(c2_trans,v_2)

				self.c1_config[i][0] = v_res_1.item(0); self.c1_config[i][1] = v_res_1.item(1); self.c1_config[i][2] = v_res_1.item(2)
				self.c2_config[i][0] = v_res_2.item(0); self.c2_config[i][1] = v_res_2.item(1); self.c2_config[i][2] = v_res_2.item(2)

	def show_info(self,*args):
		'''
			DESCROIPTION:
		'''

		if 'config' in args:
			 print("#"*80); print("CONFIG CLUSTER 1"); print("#"*80)
			 for i in range(self.c1_n):
				 print("%2s%12.6lf%12.6lf%12.6lf" % (self.c1_config[i][3], \
					 self.c1_config[i][0],self.c1_config[i][1], \
					 self.c1_config[i][2]))
			 print("#"*80); print("CONFIG CLUSTER 2"); print("#"*80)
			 for i in range(self.c2_n):
				 print("%2s%12.6lf%12.6lf%12.6lf" % (self.c2_config[i][3], \
					 self.c2_config[i][0],self.c2_config[i][1], \
					 self.c2_config[i][2]))
		if 'com' in args:
			print("#"*80); print("CENTRE OF MASS CLUSTER 1 : ",end="")
			print("%.6lf%12.6lf%12.6lf" % ( self.c1_com[0], self.c1_com[1], self.c1_com[2] ))
			print("#"*80)
			print("#"*80); print("CENTRE OF MASS CLUSTER 2 : ",end="")
			print("%.6lf%12.6lf%12.6lf" % ( self.c2_com[0], self.c2_com[1], self.c2_com[2] ))
			print("#"*80)
		if 'moi' in args:
			print("#"*80); print("INERTIA TENSOR CLUSTER 1 : ")
			for i in range(3):
				print("%12.6f%12.6f%12.6f" % (self.c1_moi_np.item((i,0)),\
					self.c1_moi_np.item((i,1)),self.c1_moi_np.item((i,2))))
			print("#"*80); print("INERTIA TENSOR CLUSTER 2 : ")
			for i in range(3):
				print("%12.6f%12.6f%12.6f" % (self.c2_moi_np.item((i,0)),\
					self.c2_moi_np.item((i,1)),self.c2_moi_np.item((i,2))))
			print("#"*80)





if __name__ == '__main__':
	'''
		DESCROIPTION:
	'''

	test = cf_cluster_mod(sys.argv[1],sys.argv[2])
	test.load_xyz()
	#test.show_info('config')
	test.get_com()
	#test.show_info('com')
	test.get_shift()
	test.get_com()
	#test.show_info('com')
	#test.show_info('config')
	test.get_moi()
	#test.show_info('moi')

	deg = test.get_deg()
	print("cluster 2 degeneracy : %d" % (deg))

	test.get_eig(deg)
	test.get_tra(deg)
	test.get_com()
	test.get_shift()
	test.get_com()
	test.show_info('config')
	test.get_moi()
	test.show_info('moi')








