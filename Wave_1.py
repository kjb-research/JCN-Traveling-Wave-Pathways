from brian2 import *
from math import exp, floor, sqrt
from matplotlib.pyplot import *   
from setup import izzy, pos_3D, STDP1, reg_synapse, write_video, vectormap

def Wave(seednum, batch, RATE):

    #setting the random seed based on input
    seed(seednum)

    #names for numpy '.npy' datafiles (missing numbers are due to the removal of a few originally measured quantities)
    vidName = str(batch) +  '/' + str(batch) + str(seednum) + '.avi'
    fileName1 = str(batch) +  '/' + str(batch) + str(seednum) + 'rate.npy'
    fileName2 = str(batch) +  '/' + str(batch) + str(seednum) + 'trials.npy'
    fileName3 = str(batch) +  '/' + str(batch) + str(seednum) + 'div.npy'
    fileName4 = str(batch) +  '/' + str(batch) + str(seednum) + 'vector.npy'
    fileName5 = str(batch) + '/' + str(batch) + str(seednum) + 'speed.npy'
    fileName6 = str(batch) + '/' + str(batch) + str(seednum) + 'order.npy'
    fileName9 = str(batch) + '/' + str(batch) + str(seednum) + 'dw_dist.npy'
    fileName10 = str(batch) + '/' + str(batch) + str(seednum) + 'dw_std.npy'

    #Define dimensionality of Network
    Height = 100
    Width = 100
    Layers = 3
    N = Height * Width * Layers

    #Define length of simulation
    sim_t = 100000*ms

    #connectivity constants
    #K - range of possible weights
    K = 11.0
    #C - controls the maximum connection probability
    C = 0.6
    #lambda - controls the characteristic length of neuronal connectivity
    lamb = 2.5

    #Proportion excitatory neurons
    Percent_ecn = 0.8

    #Time-decay of current inputs
    tau = 4*ms

    #Sub-Threshold Background (random stimulation)
    STIMRANGE = 0.5
    stim_freq = 100.0*Hz

    #Center Stimulation
    SpecStim = 4.0
    SpecFreq = 500.0*Hz

    #Multiplier for scale of distance dependent time-delay
    DELAYMULT = 0.5

    #Proportion of space neurons can move between one another (random positioning)
    rand_pos = 0.0

    #STDP Constants for Excitatory Neurons
    tau_pre = 16.8*ms
    tau_post = 33.7*ms 
    rate = RATE
    Apre = 0.0016*rate
    Apost = -Apre
    wmax = 0.5*K

    print('Izhikevich Model...')
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

    #setting excitatory and inhibitory neuron parameters
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

    #Initialize voltage and recovery variable values
    G.v = izzy.rest
    G.u = izzy.restu

    #Create positional x,y,z coordinates w/random displacement from lattice positions
    G.x = pos_3D.x
    G.y = pos_3D.y
    G.z = pos_3D.z

    #Create a single poisson point process per neuron for random stimulation (background)
    P = PoissonGroup(N, stim_freq)

    #create poisson inputs for the bursts of central stimulation
    P1_Spec = PoissonGroup(192, SpecFreq)

    #Connect 1 point process to each neuron. Different range of stimulation currents for inhibitory and excitatory
    PS = Synapses(P,G,on_pre = izzy.stim_type)
    PS.connect(j='i')

    #connecting the central stimulation spike trains to the correct neurons (only works for 100X100)
    P1S = Synapses(P1_Spec, G, on_pre = izzy.stim_type3)
    P1S.connect(j = 'int((Width*Height)*(floor(i/64))+ (46 + floor(i/8)%8)*Width + (46 + i%8))')

    #creating the excitatory (SE) and inhibitory (SI) synapses
    SE = Synapses(G, G, STDP1.eqns, on_pre = STDP1.on_pre, on_post = STDP1.on_post)
    SI = Synapses(G,G, reg_synapse.eqns, on_pre = reg_synapse.on_pre)

    #distance-dependent connection probability
    SE.connect(condition='i!=j', p = izzy.SE_prob)
    SI.connect(condition='i!=j', p = izzy.SI_prob)

    #different range of randomly initialized weights for excitatory and inhibitory
    SE.w = izzy.SE_w
    SI.w = izzy.SI_w

    #distance dependent time-delay
    SI.delay = izzy.delay
    SE.delay = izzy.delay

    #saving initial weights and indices for vectorfield plot
    SEpre = []
    SEi = []
    SEj = []
    for i in range (len(SE.w)):
        SEpre.append(SE.w[i])
        SEi.append(SE.i[i])
        SEj.append(SE.j[i])

    #creating empty lists for various saved quantities
    trials = []
    div = []
    rate = []
    dir_save = []
    speeds = []
    dw_dist = []
    dw_std = []
    
    #directional local order calculation
    div.append(izzy.props_div(SE.w,K))
    nonzeros = (SE.w != 0.0)
    nonzero = np.sum(nonzeros)
    dir_save.append(vectormap.dir_div_3D(SE.i,SE.j,SE.w,G.x,G.y,G.z,Width,Height,Layers))

    Nsteps = int((sim_t/ms)/1000) + 1
    Rates = PopulationRateMonitor(G)
    time = 0.0
    for i in range(Nsteps):
        trial = i+1
        
        Wpre = []
        for ll in range(len(SE.w)):
            Wpre.append(SE.w[ll])
        Wpre = np.array(Wpre)
        

        SPIKES = SpikeMonitor(G)
        run (30 * ms)
        Rates = PopulationRateMonitor(G)
        if trial%10 == 0:
            #voltage saving for periodic plotting of behavior
            Voltage = StateMonitor(G, 'v', record=True)
        run(100*ms)
        rate.append(mean(Rates.rate/Hz))
        speeds.append(izzy.speed(time+30,SPIKES.i,SPIKES.t/ms,G.x,G.y,G.z))
        del(Rates)
        del(SPIKES)
        run(50*ms)
        if trial%10 == 0:
            #plotting behavior
            title = str(batch) + '/trial' + str(trial) + '_seed' + str(seednum) + '.png'
            pos_3D.plot_frames_flat(0, trial*1000 + 30, 15, Width, Height, Layers, Voltage.v, 'afmhot', -40, title)
            del(Voltage)
        run(820*ms)
        div.append(izzy.props_div(SE.w,K))
        dir_save.append(vectormap.dir_div_3D(SE.i,SE.j,SE.w,G.x,G.y,G.z,Width,Height,Layers))
        trials.append(trial)
        
        Wpost = []
        for ll in range(len(SE.w)):
            Wpost.append(SE.w[ll])
        Wpost = np.array(Wpost)

            
        dwvals = Wpost - Wpre
        dw_dist.append(np.mean(dwvals))
        dw_std.append(np.std(dwvals))

        time += 1000.0
        
    #saving various measured quantities
    np.save(fileName1, rate)
    np.save(fileName2, trials)
    np.save(fileName3, div)
    np.save(fileName5, speeds)
    np.save(fileName6, dir_save)
    np.save(fileName9, dw_dist)
    np.save(fileName10, dw_std)

    deltaSE = []              
    for j in range(len(SE.w)):
        deltaSE.append(SE.w[j] - SEpre[j])
    
    np.save(fileName4, [deltaSE,SEi,SEj])







