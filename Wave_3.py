from brian2 import *
from math import exp, floor, sqrt
from matplotlib.pyplot import *   
from setup import izzy, pos_3D, STDP1, reg_synapse, write_video, vectormap, wave_propagation

def Wave(seednum,batch, RATE):
    seed(seednum)
    
    vidName = str(batch) +  '/' + str(batch) + str(seednum) + '.avi'
    fileName1 = str(batch) +  '/' + str(batch) + str(seednum) + 'rate.npy'
    fileName2 = str(batch) +  '/' + str(batch) + str(seednum) + 'trials.npy'
    fileName3 = str(batch) +  '/' + str(batch) + str(seednum) + 'div.npy'
    fileName4 = str(batch) +  '/' + str(batch) + str(seednum) + 'vector.npy'
    fileName6 = str(batch) + '/' + str(batch) + str(seednum) + 'order.npy'
    
    #Define dimensionality of Network
    Height = 100
    Width = 100
    Layers = 3
    N = Height * Width * Layers

    #Define length of simulation
    sim_t = 100000*ms

    #connectivity constants
    #K - controls the initial range of weights
    K = 11.0
    #C - controls the maximum connection probability
    C = 0.6
    #lambda - controls the characteristic connectivity length
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
    print('Setting model parameters...')
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


    #Initialize voltage values
    G.v = izzy.rest
    G.u = izzy.restu

    #Create positional x,y,z coordinates w/random displacement from lattice positions
    G.x = pos_3D.x
    G.y = pos_3D.y
    G.z = pos_3D.z

    #Create a single poisson point process per neuron for random stimulation
    P = PoissonGroup(N, stim_freq)
    
    #creating poisson point process spike trains for the top and bottom stimulation locations
    P1_Spec = PoissonGroup(192, SpecFreq)
    P2_Spec = PoissonGroup(192, SpecFreq)


    #Connect 1 point process to each neuron. Different range of stimulation currents for inhibitory and excitatory
    PS = Synapses(P,G,on_pre = izzy.stim_type)
    PS.connect(j='i')

    #connecting point stimulation to the correct areas
    P1S = Synapses(P1_Spec, G, on_pre = izzy.stim_alt21)
    P1S.connect(j = 'int((Width*Height)*(floor(i/64)) + (46+floor(i/8)%8)*Width + (71 + i%8))')
    P2S = Synapses(P2_Spec, G, on_pre = izzy.stim_alt22)
    P2S.connect(j = 'int((Width*Height)*(floor(i/64)) + (46+floor(i/8)%8)*Width + (21 + i%8))')

    #create excitatory (SE) and inhibitory (SI) synapses between neurons
    SE = Synapses(G, G, STDP1.eqns, on_pre = STDP1.on_pre, on_post = STDP1.on_post)
    SI = Synapses(G,G, reg_synapse.eqns, on_pre = reg_synapse.on_pre)

    #distance-dependent connection probability
    SE.connect(condition='i!=j', p = izzy.SE_prob)
    SI.connect(condition='i!=j', p = izzy.SI_prob)

    #different range of random weights for excitatory and inhibitory
    SE.w = izzy.SE_w
    SI.w = izzy.SI_w

    #distance dependent time-delay
    SI.delay = izzy.delay
    SE.delay = izzy.delay

    #saving the initial synaptic weights
    SEpre = []
    SEi = []
    SEj = []

    for i in range (len(SE.w)):
        SEpre.append(SE.w[i])
        SEi.append(SE.i[i])
        SEj.append(SE.j[i])

    #empty lists for saving of various quantities
    #stimulation index
    trials = []
    #divergence of weights
    div = []
    #network population firing rates 100ms after each stimulation
    rate = []
    #local directional order of synaptic weights
    dir_save = []

    div.append(izzy.props_div(SE.w,K))
    dir_save.append(vectormap.dir_div_3D(SE.i,SE.j,SE.w,G.x,G.y,G.z,Width,Height,Layers))
    
    Nsteps = int((sim_t/ms)/1000)
    time = 0.0
    for i in range(Nsteps):
        trial = i+1
        run(30*ms)
        if trial%10 == 0 or trial%10 == 1:
            Voltage = StateMonitor(G, 'v', record=True)
        Rates = PopulationRateMonitor(G)
        run(100*ms)
        rate.append(mean(Rates.rate/Hz))
        del(Rates)
        run(50*ms)
        if trial%10 == 0 or trial%10 == 1:
            title = str(batch) + '/trial' + str(trial) + '_seed' + str(seednum) + '.png'
            pos_3D.plot_frames_flat(0, trial*1000 + 30, 15, Width, Height, Layers, Voltage.v, 'afmhot', -40, title)
            del(Voltage)
        run(820*ms)
        div.append(izzy.props_div(SE.w,K))
        dir_save.append(vectormap.dir_div_3D(SE.i,SE.j,SE.w,G.x,G.y,G.z,Width,Height,Layers))
        trials.append(trial)
        print('stim trial: ', trial)
        time += 1000.0

    np.save(fileName1, rate)
    np.save(fileName2, trials)
    np.save(fileName3, div)
    np.save(fileName6, dir_save)
    
    deltaSE = []              
    for j in range(len(SE.w)):
        deltaSE.append(SE.w[j] - SEpre[j])
    
    np.save(fileName4, [deltaSE,SEi,SEj])







