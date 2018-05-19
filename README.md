# Multicellularity
Visualisation of a model for multicellular systems based on secrete-and-sense cells.  

### Description
How do communicating cells organize themselves spatially and temporally? We modelled cells that secrete and sense the same signalling molecule as a discrete dynamical system. Through tuning properties such as secretion rate, effective distance between cells and noise level we were able to identify various phases, in some of which cells behave autonomously whereas in others they behave collectively. Our goal here is to provide a direct visualisation of the model for you to get acquainted with our model. 

More information on the model and research can be found [here](http://youklab.org/research.html). The papers in which the model was first published and further developed are:
*  T. Maire and H. Youk. [Molecular-level tuning of cellular autonomy controls collective behaviors of cell populations](http://www.youklab.org/papers/CellSystems2015_Maire.pdf). Cell Systems 1, 349–360 (2015).
* E. P. Olimpio*, Y. Dang*, and H. Youk. [Statistical dynamics of spatial-order formation by communicating cells](http://www.youklab.org/papers/iScience2018_Olimpio_Dang.pdf). iScience 2: 27-40 (2018).

### **For additional information, troubleshooting and a project plan please refer to the [Wiki page](https://github.com/YitengD/Multicellularity/wiki).**

## MATLAB app
### Installation
#### Requirements

MATLAB (R2017b, recent earlier versions should also work)

#### Instructions

Using the App Installer:
1. Download the file 'Multicellularity R*.mlappinstall', where * is the version number of the latest release.
2. Run the MATLAB Installer file. The program should unpack itself and install as an app in MATLAB. 
3. Run the app from the "Apps" menu in MATLAB.

Manually:
1. Download the file 'Multicellularity R*.zip', where * is the version number of the latest release. 
2. Extract its contents and set it to a folder in which you have write permission (to store files).
3. Run the file "app_R*_final.mlapp* to directly run the app.
4. To edit the app, open "app_R*_final.mlapp* from the MATLAB "Open folder" menu or after you start MATLAB App Designer.

### Usage
The GUI consists of (1) a menu at the top, (2) a left box with different tabs, with editable fields and switches under each tab, (3) a figure on the right which will output the lattice, (4) buttons on the bottom right for simulation control and (5) a message box for output messages. Let's go through the components one by one.

#### Visualisation

The graph on the right becomes a visualisation of the multicellular system once you run the simulation. The circles represent the cells and the color their expression level. Generally, a darker color corresponds to a higher expression level. For the binary system, OFF cells are white and ON cells are black. In the continuous system, the darker of the shade of grey, the higher the expression level of the cell. 
The color bar on the right shows how the colors match to the exact states of the cells. The caption above the figure displays the time (in units of time steps) and statistics of the current figure. Here, p is the mean expression level of the cells and I (optionally displayed) is a spatial index measuring the degree of spatial organisation of the cells (see the Olimpio et al. paper).

#### Running the simulation

The controls for running the simulation are based on media playback buttons and should be intuitive to use. Click on the "PLAY" button to run a new simulation. The "PAUSE" button will pause the simulation and resume once you press "PLAY" again. The simulation continues until the system reaches equilibrium (when none of the cells changes state upon updating), or when the maximum simulation time is reached (by default set to tmax = 5000). To stop a current simulation or reset a terminated simulation, press the "STOP" button. This erases the current simulation and all data associated with it, so if you intend to save the simulation, do so before you press "STOP". 
The four buttons at the bottom are for going through a simulation without running it for more time (Replay mode). The skip forward and backward buttons skip to the start and end of the current simulation. The seek forward and backward buttons move the simulation by one step in time. At any time, the simulation can be run again by pressing "PLAY". If the simulation did not reach equilibrium, it will continue running after the last frame of the replay has been reached.

#### Display signal
For two types of signals (see below), the simulator only displays the expression of one gene at a time. The "Display signal" button group allows you to toggle between the views of the two states of the cells. 

### Tabs
There are three tabs for simulating cells that communicate in different ways.

#### Tab 1: 1 type of signal
The first tab is for simulating a system where the cells communicate through one type of signalling molecule. Detailed descriptions of the variables can be found in the papers listed above. Below is a quick overview of all parameters of the system.

1. **Cell states**. This dial controls whether the cells have a binary response (ON/OFF) or a continuous response (secretion rate continuously varying). The first corresponds to a response curve (secretion rate vs. sensed concentration) that is a step function, whereas the latter corresponds to a Hill function with finite Hill coefficient.
2. **Gridsize**. Sets the size of the lattice. The constructed lattice will be a triangular lattice with both size consisting of 'gridsize' number of cells.
3. **a0**. Effective distance between the cells. This sets the interaction strength between the cells. A small value of a0 corresponds to strongly interacting cells, whereas for large a0 cells have little influence on each other.
4. **K**. Sensing threshold. In the binary system, this sets the threshold for the concentration cells must sense to turn ON. In the continuous system, this is the sensed concentration at which the cells are exactly halfway between their lowest and highest expression levels.
5. **Con**. Maximum secretion rate. This is the secretion rate attained by ON cells in the binary system, or the maximum secretion rate in the continuous system.
6. **Noise**. Stochasticity is modelled as a stochastic term in the sensing threshold K. The value of this field sets the strength of the fluctuations of K. 
7. **Hill coefficient**. (Only for continuous system) sets the steepness of the response function. A high Hill coefficient corresponds to a sharp response around the value of K, whereas a low value corresponds to a gradual response.
8. **p(t=0)**. Initial mean expression level (values between 0 and 1). In the binary system, this corresponds to the initial fraction of ON cells. 
9. **I(t=0)**. Initial spatial index (values between -1 and 1). A value of 0 corresponds to a random distribution of cell states across space, whereas a value close to 1 (-1) corresponds to cell states that are positively (negatively) correlated with each other in space.
10. **Initial state (cont. system)**. This option applies only to the continuous system. Specify the type of initial state to generate. By default, the simulator uses a Monte Carlo algorithm ('MC randomized' button) to generate a lattice with given p(t=0) that is not uniform. Alterntively, one can choose a binary initial state ('binary' button) where the cells are all either completely ON or completely OFF at the beginning, similar to the binary cells / infinite Hill coefficient case.

#### Tab 2: 2 types of signal, multiplicative interaction
This tab is for simulating cells that communicate through two types of signalling molecule, which we assume do not directly interfere with each other. However, both are able to influence the state of the cell, and can do so in different ways. This tab is for simulating multiplicative interactions, where the effects to the the two signaling molecules multiply. This corresponds to an AND logic gate in the case of step-function responses. Biophysically, multiplicative interactions can happen when there is a gene with multiple promoters, that requires activators to be bound to all promoters before the gene is transcribed. Concretely, an example of multiplicative interaction is a cell that upregulates its expression of gene 1 only if it senses a high concentration of both signal 1 and signal 2. If one of the two has a low level, the cell will not upregulate its expression level. 

In this case, there is a number of changes in the parameters that describe the system:
1. **lambda12**. This is the ratio between the diffusion length of the two signalling molecules. The diffusion length sets a typical length scale for the decay of the signalling molecule surrounding a given cell. All lengths in the system are expressed with respect to lambda1, the diffusion length of molecule 1. 
2. **K** = K_ij. Instead of a having a single threshold K (for infinite Hill coefficient), we now have a threshold for each interaction of the system. With two molecules, there are four interactions.
3. **Con** = Con_i. Likewise, the two signalling molecules can have different maximum production rates. In addition, there is also the 'OFF' production rate for each gene, Coff_i, which can be different from each other. For simplicity, however, we set Coff_i=1 for all genes.
4. **p1(t=0), p2(t=0)**. We define a mean expression level for each of the two genes. As initial condition we can then set the average gene expression for both genes. 
5. **I1(t=0), I2(t=0)**. Likewise, we define spatial indices for both types of molecules. The two spatial indices are independent of each other for any fixed configuration of the system.

#### Tab 3: 2 types of signal, additive interaction
The two types of molecules can also interact additively, where the effects of sensing the molecules 
add up to each other. Biophysically, additive interactions can happen when there are multiple copies of the same gene with different promoters that do not interact. Concretely, an example of additive interaction is a cell that upregulates its expression of gene 1 whenever it senses a high concentration of either signal 1 or signal 2. If both have a high level, the cell will also upregulate its expression level. 

The additive interaction has one change in parameter specification compared to the multiplicative interaction. 
1. **Con** = Con_ij. Instead of two parameters for Con, we now have four. There are still only two molecules being secreted, but we now have to specify the upregulation due to each of the signalling molecules separately. For instance, if both activators of gene 1 are present, the secretion rate of molecule 1 will be Con,11+Con,12. If only activator 1 is present, it would be Con,11+Coff,12 = Con,11+1. 

### Menu 
The menu on top has a number of options that are useful for plotting results, importing and exporting data and tuning settings of the simulations. 

#### File
* **Save trajectory**. Saves a simulation trajectory as a '.mat' file. The simulation data is stored as a cell array 'CellsHist' that contains the state of the cell at each time step.
* **Open trajectory**. Opens a saved trajectory for replaying. The opened trajectory will stay in memory until the "STOP" (reset) button is pushed or a change is made to the parameter tabs. 
* **Save figures (pdf)**. Saves all opened plots as '*.pdf' after user confirmation for each plot.
* **Save figures (eps)**. Saves all opened plots as '*.eps' after user confirmation for each plot.
* **Close all figures**. Closes all plots generated by the user (see 'Plot' menu).
* **Close app**. Closes all figures and close the app. It is recommended that you close the app in this way, because closing it by pressing the close button will not close all figures. 

#### Simulation
This contains the same controls as the buttons (described under 'Running the simulation') below the figure. It also show the shortcuts that can be used in place of pressing the buttons.

#### Plot
* **Plot p(t)**. Plots a trajectory of p(t) against time for the running simulation. If there are two signals, the plot will contain two lines on the same plot.
* **Plot I(t)**. Same as above, for I(t).
* **Plot Theta(t)**. Same as above, for Theta(t)/fN.
* **Plot (p(t),I(t))**. Plots a trajectory in the (p,I) plane. The starting points are marked by coloured circles, while the endpoints are marked by circles with a cross. If there are two signals, both trajectories will be plotted.
* **Plot (p(t),Theta(t))**. See above, for (p,Theta/fN).
* **Plot h(t)** (only for 1 signal, binary cell states). Plots the pseudo-energy h(t).

#### Options
* **Show I**. Selecting this option will display the value of I on top of the figure. This might cause the simulation to run slightly slower.
* **Initiate with I(t=0)** (only for binary cells). Use the value of I(t=0) - or I1(t=0) and I2(t=0) in case of two signals - for the initial state. If unchecked (default), the system will generate an initial lattice by randomly selecting ON cells according to the value of p(t=0), without regard to the initial spatial order. This typically gives a value for I of around zero, but to be precise (up to 0.01), choose this option. The simulation will then run an algorithm to try to generate a lattice with the input I(t=0). As this is not always possible, the messages will display whether the outcome has been achieved or not. Note that for very high I, this can be a slow process.

#### Help
* **Documentation**. Opens the GitHub page https://github.com/YitengD/Multicellularity.

### Troubleshooting
A list of known issues can be found under the issues page (https://github.com/YitengD/Multicellularity/issues). In general, if the simulation gets stuck, try pressing the "STOP" (reset) button to see if it resolves the problem. Restarting the app sometimes also does magic. If you are running the app on a laptop, the responses to user input can be slow and the simulation might not run smoothly. Note also that the app needs to be able to save data to your folders to run. If you installed using the installation package, the folder is likely the default folder in which MATLAB stores apps. If you installed manually, choose a folder where you have write permission to deploy the app. Feel free to report any other issues you find on the issues page. 

### Credits
Yiteng Dang

### License
MIT Licence

### Last update 19/05/2018
