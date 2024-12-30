from brian2 import *
import matplotlib.pyplot as plt
import numpy as np
from math import cos, tan, sin, exp, floor, sqrt


#Izhikevich model with exponentially decaying current inputs:
class izzy:
    eqns = '''
    dv/dt = (0.04*v**2 + 5*v + 140 - u + I)/ms : 1
    du/dt = (a*(b*v - u))/ms : 1
    dI/dt = (0 - I)/tau : 1
    a : 1
    b : 1
    c : 1
    d : 1
    x : 1
    y : 1
    z : 1
    ecn : 1
    '''

    #neuron reset behavior
    res = 'v = c ; u = u + d'

    #neuron threshold
    thresh = 'v>30'
    
    #random indices for excitatory and inhibitory neurons
    def ecn_inhib_init(N, Percent_ecn):
        excitatory = []
        inhibitory = []
        
        for ii in range (0,N):
            if rand()>Percent_ecn:
                inhibitory.append(ii)
            else:
                excitatory.append(ii)
                
        return excitatory, inhibitory
    
    #assigning of behavior variables for excitatory (ecn) and inhibitory (inhib) neurons
    a_inhib = '0.02 + 0.08 * rand()'
    b_inhib = '0.25 - 0.05 * rand()'
    c_inhib = '(-65)'
    d_inhib = '2.0'
    a_ecn = '0.02'
    b_ecn = '0.2'
    c_ecn = '(-65 + 15 * rand()**2)'
    d_ecn = '8 - 6 * rand()**2'
    
    #calculation of the resting state values for 'v' and 'u'
    rest = '((b-5) - sqrt((5-b)**2 - 4*0.04*140))/0.08'
    restu = 'v * b'
    
    ###lower stimulation currents for inhibit'ory neurons
    #regular random stimulation w/random magnitude
    stim_type = 'I_post += rand()*STIMRANGE*(ecn_post + (2/5)*(1-ecn_post))'
    #time-specific stimulation
    stim_type3 = 'I_post += ((t/ms)%(1000) < 30 and (t/ms)%(1000) > 0) * SpecStim*(ecn_post + (2/5)*(1-ecn_post))'

    #2 alternating sites
    stim_alt21 = 'I_post += ((t/ms)%(2000) < 30 and (t/ms)%(2000) > 0) * SpecStim*(ecn_post + (2/5)*(1-ecn_post))'
    stim_alt22 = 'I_post += ((t/ms + 1000)%(2000) < 30 and (t/ms + 1000)%(2000) > 0) * SpecStim*(ecn_post + (2/5)*(1-ecn_post))'
    
    #distance dependent connection probability 
    SE_prob = "ecn_pre*C*exp(-((x_post - x_pre)**2 + (y_post - y_pre)**2 + (z_post - z_pre)**2)/lamb**2)"
    SI_prob = "(1-ecn_pre)*C*exp(-((x_post - x_pre)**2 + (y_post - y_pre)**2 + (z_post - z_pre)**2)/lamb**2)"
    
    #distance-dependant probability if distinction between excitatory and inhibitory is not necessary
    allprob = 'connprob'
    gen_prob = "C*exp(-((x_post - x_pre)**2 + (y_post - y_pre)**2 + (z_post - z_pre)**2)/lamb**2)"
    
    #random range of initial weights
    SE_w = '0.5*K*rand()' 
    SI_w = '-K*rand()'
    
    def speed(start, SPIKES_i, SPIKES_t, G_x, G_y, G_z):
        speeds = []
        induse = 0
        time = SPIKES_t[induse]

        Search = True
        while Search == True:
            induse += 1
            time = SPIKES_t[induse]
            if time >= start:
                Search = False
        
        avedist = []
        for i in range(16):
            dists = []
            current = time
            while time < current + 5:
                neuron = SPIKES_i[induse]
                xval = abs(G_x[neuron] - 49.5)
                yval = abs(G_y[neuron] - 49.5)
                zval = abs(G_z[neuron] - 1)
                dist = sqrt(xval**2 + yval**2 + zval**2)
                dists.append(dist)
                induse += 1
                time = SPIKES_t[induse]
            if i > 0:
                speeds.append((np.mean(dists) - preval)/5)
            preval = np.mean(dists)
    
        return(np.mean(speeds))
                
    def props_div(SE,K):
        top = 0
        bottom = 0
        mid = 0
        total = len(SE)
        for i in range(total):
            if (SE[i] >= (0.5 * K)):
                top += 1
            if (SE[i] <= 0.0):
                bottom += 1
            if (SE[i] > 0.0) & (SE[i] < (0.5*K)): 
                mid += 1
        top = top/total
        bottom = bottom/total
        mid = mid/total
        return([mid,top,bottom])
    
    #distance dependent signal time-delay
    delay = 'DELAYMULT*sqrt((x_post - x_pre)**2 + (y_post - y_pre)**2 + (z_post - z_pre)**2) * ms'
        
#position classification for random quasi-2D (or 3D) lattice   
class pos_3D:
    x = '(i%Width) - (rand_pos/2) + rand_pos*rand()'
    y = '((floor(i/Width))%Height) - (rand_pos/2) + rand_pos*rand()'
    z = '(floor(i/(Height*Width))) - (rand_pos/2) + rand_pos*rand()'

    #plot frames of activity within the network
    def plot_frames_flat(start, start_time, timestep, Nx, Ny, Nz, v, colorchoice, max, name):
        sum = 0.0
        initial = int(start/0.1)
        xi = 0
        yi = 0
        step = int(timestep/0.1)
        neurons = Nx * Ny * Nz
        plt.figure(figsize=(27, 3), dpi=80)
        for jj in range (0,9):
            time = initial + step * jj
            current_time = start_time + start + jj * timestep
            VoltageArray = np.zeros((Nx,Ny))
            for ii in range (0,neurons-1):
                xi = int(ii%Nx)
                yi = int(floor(ii/Nx) - Ny*floor(ii/(Nx*Ny)))
                if ii<(Ny*Nx):
                    VoltageArray[xi][yi] = v[ii][time]
                if ii>= (Ny*Nx):
                    VoltageArray[xi][yi] = (VoltageArray[xi][yi] + v[ii][time])/2                    
            title = 'T = ' + str(current_time) + ' ms'        
            plt.subplot(1,9,jj+1)
            plt.pcolormesh(VoltageArray,cmap=colorchoice,vmin=-80,vmax=max)
            plt.title(title, fontsize = '14')
            if jj + 1 == 9:
                plt.colorbar()
        plt.savefig(name)
        plt.show()
    
#function for creating vectormaps of weight change:
class vectormap:
    # function for a directional divergence parameter
    def dir_div_3D(SE_i, SE_j, SE_w, G_x, G_y, G_z, Nx , Ny , Nz):
        W = np.zeros((Nx,Ny,Nz))
        Wy = np.zeros((Nx,Ny,Nz))
        Wx = np.zeros((Nx,Ny,Nz))
        Wz = np.zeros((Nx,Ny,Nz))
        Nvals = np.zeros((Nx,Ny,Nz))
       
        for i in range(len(SE_i)):
            ival = SE_i[i]
            jval = SE_j[i]
            ipos = np.array([int(G_x[ival]),int(G_y[ival]),int(G_z[ival])])
            jpos = np.array([int(G_x[jval]),int(G_y[jval]),int(G_z[jval])])
            dpos = jpos - ipos
            dphat = dpos/(sqrt(dpos[0]**2 + dpos[1]**2 + dpos[2]**2))
            wvec = dphat * SE_w[i]
            W[ipos[0],ipos[1],ipos[2]] += SE_w[i]
            Wx[ipos[0],ipos[1],ipos[2]] += wvec[0]
            Wy[ipos[0],ipos[1],ipos[2]] += wvec[1]
            Wz[ipos[0],ipos[1],ipos[2]] += wvec[2]
            Nvals[ipos[0],ipos[1],ipos[2]] += 1
        
        for j in range(Nx):
            for k in range(Ny):
                for b in range(Nz):
                    if W[j,k,b] != 0:
                        if Wx[j,k,b] != 0:
                            #Wx[j,k,b] = Wx[j,k,b]/(sqrt(Wx[j,k,b]**2 + Wy[j,k,b]**2 + Wz[j,k,b]**2))
                            Wx[j,k,b] = Wx[j,k,b]/(sqrt(Wx[j,k,b]**2 + Wy[j,k,b]**2))
                        if Wy[j,k,b] != 0:
                            #Wy[j,k,b] = Wy[j,k,b]/(sqrt(Wx[j,k,b]**2 + Wy[j,k,b]**2 + Wz[j,k,b]**2))
                            Wy[j,k,b] = Wy[j,k,b]/(sqrt(Wx[j,k,b]**2 + Wy[j,k,b]**2))
                        #if Wz[j,k,b] != 0:
                        #    Wz[j,k,b] = Wz[j,k,b]/(sqrt(Wx[j,k,b]**2 + Wy[j,k,b]**2 + Wz[j,k,b]**2))
        
        
        dir_div = []
        
        for x in range(Nx):
            for y in range(Ny):
                for z in range(Nz):
                    if W[x,y,z]!=0:
                        if y > 1 and x > 1 and y < Ny-2 and x < Nx-2:
                            if z == 0:
                                for ll in range(3*3*2):
                                    if W[x,y,z] != 0:
                                        xp = ll%3
                                        yp = floor(ll/3)%3
                                        zp = floor(ll/9)
                                        if ll != 4 and W[x-1+xp, y-1+yp, z+zp] != 0:
                                            xdot = Wx[x-1+xp, y-1+yp, z+zp] * Wx[x,y,z]
                                            ydot = Wy[x-1+xp, y-1+yp, z+zp] * Wy[x,y,z]
                                            #zdot = Wz[x-1+xp, y-1+yp, z+zp] * Wz[x,y,z]
                                            #dir_div.append(xdot+ydot+zdot)
                                            dir_div.append(xdot+ydot)
                            if z == 1:
                                for ll in range(3*3*3):
                                    if Nvals[x,y,z] != 0:
                                        xp = ll%3
                                        yp = floor(ll/3)%3
                                        zp = floor(ll/9)
                                        if ll != 13 and W[x-1+xp, y-1+yp, z-1+zp] != 0:
                                            xdot = Wx[x-1+xp, y-1+yp, z-1+zp] * Wx[x,y,z]
                                            ydot = Wy[x-1+xp, y-1+yp, z-1+zp] * Wy[x,y,z]
                                            #zdot = Wz[x-1+xp, y-1+yp, z-1+zp] * Wz[x,y,z]
                                            #dir_div.append(xdot+ydot+zdot)
                                            dir_div.append(xdot+ydot)
                            if z == 2:
                                for ll in range(3*3*2):
                                    if Nvals[x,y,z] != 0:
                                        xp = ll%3
                                        yp = floor(ll/3)%3
                                        zp = floor(ll/9)
                                        if ll != 13 and W[x-1+xp, y-1+yp, z-1+zp] != 0:
                                            xdot = Wx[x-1+xp, y-1+yp, z-1+zp] * Wx[x,y,z]
                                            ydot = Wy[x-1+xp, y-1+yp, z-1+zp] * Wy[x,y,z]
                                            #zdot = Wz[x-1+xp, y-1+yp, z-1+zp] * Wz[x,y,z]
                                            #dir_div.append(xdot+ydot+zdot)
                                            dir_div.append(xdot+ydot)
                                    
        dir_corrected = []
        for i in range(len(dir_div)):
            if dir_div[i] != float('nan'):
                dir_corrected.append(dir_div[i])
        dir_return = np.mean(dir_corrected)
        return(dir_return)
    
    #vectormap for quasi-2D example
    def vector_3D(Delt_SE, i, j, Nx, Ny, Nz, name):
        size = len(Delt_SE)
        length = int(Ny/5)
        width = int(Nx/5)
        x_pos = np.zeros(length * width)
        y_pos = np.zeros(length * width)
        x_mag = np.zeros(length * width)
        y_mag = np.zeros(length * width)
        N_array = np.zeros(length * width)
        for bb in range (length*width):
            x = bb%(width)
            y = floor(bb/width)
            x_pos[bb] = x
            y_pos[bb] = y
        for kk in range(int(size)):
            deltaW = Delt_SE[kk]
            xpre = i[kk]%Nx            
            xpost = j[kk]%Nx
            deltx = xpost - xpre
            ypre = floor(i[kk]/(Nx))%Ny
            ypost = floor(j[kk]/(Nx))%Ny
            delty = ypost - ypre
            zpre = floor(i[kk]/(Nx * Ny))
            zpost = floor(j[kk]/(Nx * Ny))
            deltz = zpost - zpre
            distxy = sqrt(deltx**2 + delty**2)
            dist = sqrt(deltx**2 + delty**2 + deltz**2)
            XW = 0.0
            YW = 0.0
            index = 0
            
            if delty == 0 and deltx == 0:
                continue
            
            XW = deltaW * (distxy/dist) * (deltx/distxy)
            YW = deltaW * (distxy/dist) * (delty/distxy)
            
            yspace = floor(i[kk]/(Nx * 5)) - floor(Ny/5) * floor(i[kk]/(Ny * Nx))   
            xspace = floor((i[kk]%Nx)/5)
            index = int(xspace + yspace * ((Nx)/5))
            
            x_mag[index] += XW
            y_mag[index] += YW
            N_array[index] += 1.0
            
        for ll in range (len(x_mag)):
            mag = sqrt(x_mag[ll]**2 + y_mag[ll]**2)
            x_mag[ll] = x_mag[ll] / N_array[ll]
            y_mag[ll] = y_mag[ll] / N_array[ll]
            
        for i in range(len(y_pos)):
            y_pos[i] = y_pos[i] * 5
            x_pos[i] = x_pos[i] * 5
            
        xtick = np.arange(0,110,10)
        ytick = np.arange(0,110,10)
            
        plt.figure(figsize=(10, 10), dpi=80)
        plt.quiver(y_pos, x_pos, y_mag, x_mag)
        plt.title(r'$\vec{\Delta W}$', fontsize = '26')
        plt.xlabel('X', fontsize = '20')
        plt.ylabel('Y', fontsize = '20')
        plt.xticks(xtick, fontsize = '18')
        plt.yticks(ytick, fontsize = '18')
        plt.show()
        
        return([y_pos, x_pos, y_mag, x_mag])
      
#regular current-input synapse
class reg_synapse:
    eqns = '''
    w : 1
    '''
    on_pre='''
    I_post += w
    ''' 
#Regular STDP, no Inhibitory Rules, No Metaplasticity
class STDP1:
    eqns = '''
    w : 1
    dapre/dt = -apre/tau_pre : 1 (event-driven)
    dapost/dt = -apost/tau_post : 1 (event-driven)
    xdist : 1
    ydist : 1
    zdist : 1
    dist : 1
    '''
    on_post='''
    apost += Apost
    w = clip(w+apre, 0, wmax)
    '''
    on_pre='''
    I_post += w
    apre += Apre 
    w = clip(w+apost, 0, wmax)
    '''
    inhib = '''
    w : 1
    '''
    pre_inhib = ''' 
    I_post += w 
    '''
   