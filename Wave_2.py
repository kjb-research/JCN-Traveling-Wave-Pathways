from brian2 import *
from math import exp, floor, sqrt
import matplotlib.pyplot as plt
import numpy as np
from setup import izzy, pos_3D, STDP1, write_video, vectormap

def Wave(seednum, batch, RATE):
    seed(seednum)
    vidName = str(batch) +  '/' + str(batch) + str(seednum)
    fileName1 = str(batch) + '/' + str(batch) + str(seednum) + 'rate.npy'
    fileName2 = str(batch) + '/' + str(batch) + str(seednum) + 'div.npy'
    fileName3 = str(batch) + '/' + str(batch) + str(seednum) + 'times.npy'
    fileName4 = str(batch) + '/' + str(batch) + str(seednum) + 'vector.npy'
    fileName6 = str(batch) + '/' + str(batch) + str(seednum) + 'order.npy'

    #Define dimensionality of Network
    Height = 100
    Width = 100
    Layers = 3
    N = Height * Width * Layers

    #Define length of simulation
    sim_t = 30000*ms

    #connectivity constants
    #K - controls range of initial weights
    K = 11.0
    #C - controls maximum connection probability
    C = 0.6
    #lambda - controls the characteristic connectivity length
    lamb = 2.5

    #Proportion excitatory neurons
    Percent_ecn = 0.8

    #Time-decay of current inputs
    tau = 4*ms

    #Range of input current and frequency from Poisson Point Processes (random stimulation)
    STIMRANGE = 1.8
    stim_freq = 180.0*Hz

    #Multiplier for scale of distance dependent time-delay
    DELAYMULT = 0.5

    #Proportion of space neurons can exist between one another (random positioning)
    rand_pos = 0.0

    #STDP Constants for Excitatory Neurons (ignore if not using an STDP synapse model)
    tau_pre = 16.8*ms
    tau_post = 33.7*ms
    rate = float(RATE)
    Apre = 0.0016*rate
    Apost = -Apre
    wmax = 0.5*K

    #define izhikevich model (current modeled as a differential equation with uniform time-decay constant tau)
    eqns = izzy.eqns
    res = izzy.res
    thresh = izzy.thresh

    a_inhib = izzy.a_inhib
    b_inhib = izzy.b_inhib
    c_inhib = izzy.c_inhib
    d_inhib = izzy.d_inhib
    a_ecn = izzy.a_ecn
    b_ecn = izzy.b_ecn
    c_ecn = izzy.c_ecn
    d_ecn = izzy.d_ecn

    #determining which neurons are excitatory
    excitatory, inhibitory = izzy.ecn_inhib_init(N, Percent_ecn)

    #Create array of neurons with pre-defined threshold/spiking behavior
    G = NeuronGroup(N, eqns, threshold = thresh, reset = res, method = 'euler')

    #if statements to set up the parameters of excitatory and inhibitory neurons based on the original Izhikevich paper
    #'simple model of spiking neurons'
    if len(inhibitory) > 0:
        G.ecn[inhibitory] = '0'
        G.a[inhibitory] = izzy.a_inhib
        G.b[inhibitory] = izzy.b_inhib
        G.c[inhibitory] = izzy.c_inhib
        G.d[inhibitory] = izzy.d_inhib

    if len(excitatory) > 0:
        G.ecn[excitatory] = '1'
        G.a[excitatory] = izzy.a_ecn
        G.b[excitatory] = izzy.b_ecn
        G.c[excitatory] = izzy.c_ecn
        G.d[excitatory] = izzy.d_ecn


    #Initialize voltage 'v' and recovery variable 'u' values
    G.v = izzy.rest
    G.u = izzy.restu

    #Create positional x,y,z coordinates w/random displacement from lattice positions
    G.x = pos_3D.x
    G.y = pos_3D.y
    G.z = pos_3D.z

    #Create a single poisson point process per neuron for random stimulation
    P = PoissonGroup(N, stim_freq)

    #Connect 1 point process to each neuron. Different range of stimulation currents for inhibitory and excitatory
    PS = Synapses(P,G,on_pre = izzy.stim_type)
    PS.connect(j='i')

    #create excitatory and inhibitory synapses between neurons
    SE = Synapses(G, G, STDP1.eqns, on_pre = STDP1.on_pre, on_post = STDP1.on_post)
    SI = Synapses(G,G, STDP1.inhib, on_pre = STDP1.pre_inhib)

    #distance-dependent connection probability
    SE.connect(condition='i!=j', p = izzy.SE_prob)
    SI.connect(condition='i!=j', p = izzy.SI_prob)

    #different range of random weights for excitatory and inhibitory
    SE.w = izzy.SE_w
    SI.w = izzy.SI_w

    #distance dependent time-delay
    SI.delay = izzy.delay
    SE.delay = izzy.delay

    #excitatory weight arrays
    SE_pre = []
    SE_i = []
    SE_j = []
    deltaSE = []

    #for saving of initial excitatory weights
    for i in range (len(SE.w)):
        SE_pre.append(SE.w[i])
        SE_i.append(SE.i[i])
        SE_j.append(SE.j[i])

    #various empty lists for saved quantities
    #the network population firing rate for each 100 ms period
    rate = []
    #starting time for each 100 ms datapoint
    times = []
    #divergence of weights
    div = []
    #local directional order of synaptic weights
    dir_save = []

    Nsteps = int((sim_t/ms)/100)
    time = 0.0 

    div.append(izzy.props_div(SE.w,K))
    
    Rates = PopulationRateMonitor(G)
    for i in range(Nsteps):
        if i ==5 or i == 100 or i == 200 or i == 250:
            #saving network activity for plots
            Voltage = StateMonitor(G,'v',record=True)
        run(100*ms)
        time += 100
        if i ==7 or i == 102 or i == 202 or i == 252:
            #plotting network activity
            name = str(batch) + '/' + str(seednum) + '_' + str(i) + '.png'
            pos_3D.plot_frames_flat(0, i*100, 15, Width, Height, Layers, Voltage.v, 'afmhot', -40, name)
            del(Voltage)
        rate.append(np.mean(Rates.rate/Hz))
        del(Rates)
        times.append(time)
        dir_save.append(vectormap.dir_div_3D(SE.i,SE.j,SE.w,G.x,G.y,G.z,Width,Height,Layers))
        div.append(izzy.props_div(SE.w,K))
        Rates = PopulationRateMonitor(G)
        print('percent complete = ', 100.0*((i+1)/Nsteps))
        

    filename11 = str(batch) +  '/' + str(batch) + str(seednum) + 'wdist.npy'
    np.save(filename11, [w1vals, w2vals, w3vals])
        
    #saving relevant quantities
    np.save(fileName1, rate)
    np.save(fileName2, div)
    np.save(fileName3, times)
    np.save(fileName6, dir_save)

    #saving change in weights for vectorplots
    for j in range(len(SE.w)):
        deltaSE.append(SE.w[j] - SE_pre[j])
    np.save(fileName4, [deltaSE, SE_i, SE_j])



