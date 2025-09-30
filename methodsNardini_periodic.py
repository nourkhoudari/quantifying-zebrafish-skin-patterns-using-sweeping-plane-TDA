import numpy as np
import pdb, glob, time, os, imageio, copy

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from ripser import ripser
from persim import plot_diagrams
import gudhi as gd
from scipy.stats import multivariate_normal
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('png', 'pdf') 


def param_sweep(N,X,filename=None,iter_num = 50,plane_dir='less',periodic_y = False):
	
	'''
	param_sweep

	inputs:
	N : Binary image
	X : independent variable vector meshgrid (nparray)
	iter_num : number of flooding events to compute
        filename: where to save the output (None to avoid saving)
        plane dir: 'less'/'greater' to look at pixels behind/ahead of plane location (string)
        periodic_y: whether or not to use periodic boundary condition for the second axis

	output
	diag : Birth and death times for all topological features
	'''

	if N.shape != X.shape:
		raise Exception("Shape of N and X must the equal")

	xm,ym = N.shape

	r_range = np.linspace(0,np.max(X),iter_num)

	if plane_dir == 'greater':
                r_range = r_range[::-1]

	#to help with indexing
	if periodic_y:
		y_int = np.arange(ym)
		y_int_p1 = (np.arange(ym)+1)%ym
	else:
		y_int = np.arange(ym-1)
		y_int_p1 = np.arange(1,ym)


        	
	st = gd.SimplexTree()

	for k,rr in enumerate(r_range):

		#find nonzero pixels in N that are to the correct direction of the plane
		if plane_dir == 'less':
			N_update = np.logical_and(N,X<=rr)
		elif plane_dir == 'greater':
			N_update = np.logical_and(N,X>=rr)
		else:
			raise Exception("N_update must be \'less\' or \'greater\'")

		## look for vertical neighbors -- create 0-simplices
		cell_loc = N_update==1
		a = np.where(cell_loc)
		a = np.hstack((a[0][:,np.newaxis],a[1][:,np.newaxis]))
		locs = a[:,0] + xm*a[:,1]
		#locs = xm*a[:,0] + a[:,1]
		for j in locs:
			st.insert([j],filtration = k)

		## look for vertical neighbors -- create 1-simplices
		vert_neighbors = np.logical_and(N_update[:-1,:]==1,N_update[1:,:]==1)
		a = np.where(vert_neighbors)
		a = np.hstack((a[0][:,np.newaxis],a[1][:,np.newaxis]))
		locs = a[:,0] + xm*a[:,1]
		#locs = xm*a[:,0] + a[:,1]
		for j in locs:
			st.insert([j,j+1],filtration = k)

		## look for horizontal neighbors
		horiz_neighbors = np.logical_and(N_update[:,y_int]==1,N_update[:,y_int_p1]==1)
		a = np.where(horiz_neighbors)
		a = np.hstack((a[0][:,np.newaxis],a[1][:,np.newaxis]))
		locs = a[:,0] + xm*a[:,1]
		
		for j in locs:
			st.insert([j,(j+xm)%(xm*ym)],filtration = k)


		#look for diagonal neighbors (top left to bottom right)
		diag_neighbors = np.logical_and(N_update[:-1,y_int]==1,N_update[1:,y_int_p1]==1)
		a = np.where(diag_neighbors)
		a = np.hstack((a[0][:,np.newaxis],a[1][:,np.newaxis]))
		locs = a[:,0] + xm*a[:,1]
		#locs = xm*a[:,0] + a[:,1]
		for j in locs:
			st.insert([j,(j+xm+1)%(xm*ym)],filtration = k)
		

		#look for diagonal neighbors (bottom left to top right)
		diag_neighbors = np.logical_and(N_update[1:,y_int]==1,N_update[:-1,y_int_p1]==1)
		a = np.where(diag_neighbors)
		a = np.hstack((a[0][:,np.newaxis],a[1][:,np.newaxis]))
		locs = a[:,0] + xm*a[:,1]
		
		for j in locs:
			st.insert([j+1,(j+xm)%(xm*ym)],filtration = k)

		st.set_dimension(2)

		###include 2-simplices (looking for four different types of corners)
	
		if periodic_y:
			y_iter_end = ym
		else:
			y_iter_end = ym-1
        

		for j in np.arange(y_iter_end):
			for i in np.arange(xm-1):

				
				#top left corner:
				if N_update[i,j]==1 and N_update[i+1,j]==1 and N_update[i,(j+1)%ym]==1:
					st.insert([i + xm*j,(i+1) + xm*j , (i + xm*(j+1))%(xm*ym)],filtration = k)
				
				#top right corner
				if N_update[i,j]==1 and N_update[i+1,j]==1 and N_update[i+1,(j+1)%ym]==1:
					st.insert([i + j*xm, (i+1)+j*xm, ((i+1)  + (j+1)*xm)%(xm*ym)],filtration = k)

				#bottom left corner
				if N_update[i,j]==1 and N_update[i,(j+1)%ym]==1 and N_update[i+1,(j+1)%ym]==1:
					st.insert([i + j*xm, (i + (j+1)*xm)%(xm*ym), ((i+1) + (j+1)*xm)%(xm*ym)],filtration = k)

				#bottom right corner
				if N_update[i+1,(j+1)%ym]==1 and N_update[i+1,j]==1 and N_update[i,(j+1)%ym]==1:
					st.insert([((i+1) + (j + 1)*xm)%(xm*ym), (i+1) + j*xm, (i + (j + 1)*xm)%(xm*ym)],filtration = k)

	
	diag = st.persistence()

	if filename is not None:

		data = {}
		data['BD'] = diag
		np.save(filename,data)


	return diag


def weight_fun_ramp(x,**options):

	'''
	Weight function for persistence images

	inputs 

	x : function input
	b : max x value

	outputs 

	y: function output
	'''

	b = options.get("b")

	y = np.zeros(x.shape)

	samp = np.where(x<=0)[0]
	y[samp] = np.zeros(samp.shape)

	samp = np.where(np.logical_and(x>0,x<b))[0]
	y[samp] = x[samp]/b

	samp = np.where(x>=b)[0]
	y[samp] = np.ones(samp.shape)

	return y


def weight_fun_1(x,**options):

	'''
	Weight function of 1's for persistence images

	inputs 

	x: function input

	outputs 

	y: function output
	'''

	y = np.ones(x.shape)

	return y



def betti_curve(diag=None,filename=None,filename_save=None,r0=0,r1=1,rN=40):

	'''
	betti_curve construction

	inputs

	diag :          Input Birth-death interval list. If none, then this will be loaded in
	filename: 		Where Birth-death interval list is stored
	filename_save	Where to save persistence image

	output

	IP : 			Persistence Image
	'''


	if diag is None:
		
		if filename is None:
			raise Exception("Either interval data or filename for one must be provided")
		mat = np.load(filename + '.npy',allow_pickle=True, encoding='latin1').item()

		diag = mat['BD']

	r_range = np.linspace(r0,r1,rN)

	b0 = np.zeros(r_range.shape)
	b1 = np.zeros(r_range.shape)

	for i,r in enumerate(r_range):
		for dd in diag:
			if r >= dd[1][0] and r < dd[1][1]:
				if dd[0] == 0:
					b0[i] += 1
				elif dd[0] == 1:
					b1[i] += 1

	if filename_save is not None:
		data = {}
		data['b0'] = b0
		data['b1'] = b1
		data['r'] = r_range
		np.savetxt(filename_save + "_b0", b0)
		np.savetxt(filename_save + "_b1", b1)
		np.savetxt(filename_save + "_r", r_range)
		#np.save(filename_save,data)

	return b0,b1,r_range



def Persist_im(diag=None,filename=None,filename_save=None,inf_val=25,sigma=1e-1,weight_fun=weight_fun_ramp):

	'''
	create persistence image

	inputs

	diag :          Input Birth-death interval list. If none, then this will be loaded in
	filename: 		Where Birth-death interval list is stored
	filename_save	Where to save persistence image

	output

	IP : 			Persistence Image
	'''
    
	if diag is None:
		
		if filename is None:
			raise Exception("Either interval data or filename for one must be provided")
		mat = np.load(filename + '.npy',allow_pickle=True, encoding='latin1').item()

		diag = mat['BD']

	#resolution of final persistance image will be res**2
	res = 50	
	
	### Convert to non-diagonal form
	BD_list = [np.zeros((1,2)),np.zeros((1,2))]

	b0 = 0
	b1 = 0
	for dd in diag:
		if dd[0] == 0:

			if b0 == 0:
				BD_list[0][0,:] = dd[1]
			else:
				BD_list[0] = np.vstack((BD_list[0],dd[1]))
			
			b0 += 1

		elif dd[0] == 1:

			if b1 == 0:
				BD_list[1][0,:] = dd[1]
			else:
				BD_list[1] = np.vstack((BD_list[1],dd[1]))

			b1 += 1

	Ip_ones = [np.zeros((res,res)),np.zeros((res,res))]
	Ip_ramp = [np.zeros((res,res)),np.zeros((res,res))]

	for i,BD in enumerate(BD_list):

		BD[np.isinf(BD)] = inf_val
		BD_adjust = np.hstack([BD[:,0][:,np.newaxis],(BD[:,1] - BD[:,0])[:,np.newaxis]])

		width,height = np.max(BD_adjust,axis=0)
		length = inf_val#np.max((width,height))
		U = BD_adjust.shape[0]

		x = np.linspace(0,length,res+1)
		y = np.linspace(0,length,res+1)

		X,Y = np.meshgrid(x,y)

		shape = X.shape

		weights_ones = weight_fun_1(BD_adjust[:,1],b=height)
		weights_ramp = weight_fun_ramp(BD_adjust[:,1],b=height)

		for j,bd in enumerate(BD_adjust):

			Ip_tmp = np.zeros((res+1,res+1))
			for k,xx in enumerate(x):
				for l,yy in enumerate(y):
					Ip_tmp[k,l] = multivariate_normal.cdf(np.hstack((xx,yy)),
															mean=bd,
															cov=sigma)
			
			#Use summed area table (coordinates reverse of those described in wikipedia)
			Ip_ones[i] +=  weights_ones[j]*(Ip_tmp[1:,1:] + Ip_tmp[:-1,:-1] - Ip_tmp[1:,:-1] - Ip_tmp[:-1,1:])
			Ip_ramp[i] +=  weights_ramp[j]*(Ip_tmp[1:,1:] + Ip_tmp[:-1,:-1] - Ip_tmp[1:,:-1] - Ip_tmp[:-1,1:])


	if filename_save is not None:
		data = {}
		data['Ip'] = Ip_ones
		np.save(filename_save[0],data)

		data = {}
		data['Ip'] = Ip_ramp
		np.save(filename_save[1],data)

	return Ip_ones,Ip_ramp
