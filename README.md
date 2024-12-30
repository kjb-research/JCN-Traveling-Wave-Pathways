# Neuronal traveling waves form preferred pathways using synaptic plasticity - Reproducibility Code
## What is this for?:
This code is meant to recreate results from the paper "Neuronal Traveling Waves Form Preferred Pathways Using Synaptic Plasticity" by Kendall Butler and Luis Cruz. It implements the BRIAN neural simulator to simulate traveling waves through a model neuronal tissue with Spike-Timing Dependent Plasticity (STDP). This research was conducted through the Physics Department at Drexel University. If you can't find your way back to the original work, it is open access at https://link.springer.com/article/10.1007/s10827-024-00890-2. 

For reference and debugging, the versions of various python modules used while putting together this research are below. If any issues or questions arise, feel free to email me directly at (kendall.butler.research@gmail.com), or email the corresponding author, Dr. Luis Cruz.

### Module Versions (as used at the time of publication - newer versions may require debugging the posted code):
- python = 3.11.9
- brian2 = 2.6.0
- numpy = 1.26.4
- matplotlib = 3.8.4

### What does the code do?:
The codes posted recreates the experiments in our paper. It creates a lattice of simulated neurons using the izhikevich model (https://ieeexplore.ieee.org/document/1257420), connected locally via simple synapses. When a neuron spikes, a delayed signal is sent to any neurons it is connected to via a synapse. Each of these has a 'weight' that controls how much the neuron receiving the signal has it's input current incremented. This locally connected lattice allows for the propagation of travelling waves. We also include STDP, which controls positive and negative changes in excitatory synaptic weights based on the timing of spikes. A presynaptic neuron is the cell sending the signal along a synapse after a spike, and the postsynaptic neuron is the cell recieving that signal. The equations describing this are below:

$$\Delta t = t_i - t_j$$ (difference in spike times between presynaptic neuron 'i' and postsynaptic neuron 'j')

$$\Delta W_{ij} = A_+ e^{-|\Delta t|/\tau_+}$$ (if $$\Delta t$$ is negative)

$$\Delta W_{ij} = A_- e^{-|\Delta t|/\tau_-}$$ (if $$\Delta t$$ is positive)

Put simply, the goal of simulating a biologically plausible network, with plasticity, that propagates travelling waves is to better understand the effect that traveling waves can have on real networks with plasticity as a mechanism. As discussed in our paper, we have found that they are capable of forming preferred synaptic pathways over which traveling waves can travel along more consistently. We also investigated the robustness of this phenomenon, and whether this can be interpreted as a mechanism of large-scale competition between initially available pathways.

There are three experiments:
1. Central Stimulation (bursts of stimulation at the center of the network every 1000 ms, meant to investigate the formation of pathways with repetetive wave propagation)
2. Stochastic (or random) Stimulation (random input spikes spread continuously over the entire network, creating less predictable wave formation and subsequent formed pathways)
3. Alternating Stimulation (alternating bursts of stimulation on opposite sides of the network, meant to investigate the competition between initially available pathways)


## How to run the code:
As you can see, there are four python files included on this repository. Three represent the code to run our three distinct experiments. The fourth is a shared file, that has code used to create the network and run various analysis. In order to run the code, you will need to download at least one of the three 'experiment' files, along with the shared 'setup' file. Here is how you can tell which is which:

- 'Wave_1.py' (python file corresponding to the central stimulation experiment)
- 'Wave_2.py' (python file corresponding to the stochastic stimulation experiment)
- 'Wave_3.py' (python file corresponding to the alternating stimulation experiment)
- 'setup.py' (shared setup file for all three of the above)

Once you have these files downloaded into a shared folder, running is simple:

1. start up python in your terminal by typing 'python'
2. Import the function from the file you want to run. (for central stimulation, this command would be 'from Wave_1 import Wave')
3. Run the function. It's first argument is the seed number ('seednum'), which is used for reproducibility (using the same seed should yield exactly equivalent results). The second argument is the batch name ('batch'), which is used to place recorded data into the correct folder. Make sure you have an empty folder with the same name in your directory. The third argument is the RATE, which is a scaling factor for the STDP coefficients. Our paper uses RATE = 0, 1, 2, 3, 4 to adjust this as an independent variable.

To recap, to run the code you will type (in terminal):

**python**

**from *filename* import Wave**

**Wave(*seed*, *batch_name*, *STDP_scaling_factor*)**

Remember, you will need a folder with the same name as your *batch_name* in order for the code to run and data to be saved in the desired location. 
