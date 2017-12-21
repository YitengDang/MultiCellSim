# Multicellularity
Visualisation of a model for multicellular systems based on secrete-and-sense cells.  

### Description
How do communicating cells organize themselves spatially and temporally? We modelled cells that secrete and sense the same signalling molecule as a discrete dynamical system. Through tuning properties such as secretion rate, effective distance between cells and noise level we were able to identify various phases, in some of which cells behave autonomously whereas in others they behave collectively. Our goal here is to provide a direct visualisation of the model for you to get acquainted with our model. 

More information on the model and research can be found [here](http://youklab.org/research.html). The papers in which the model was first published and further developed are:
*  T. Maire and H. Youk. [Molecular-level tuning of cellular autonomy controls collective behaviors of cell populations](http://www.youklab.org/papers/CellSystems2015_Maire.pdf). Cell Systems 1, 349–360 (2015).
* E. P. Olimpio, Y. Dang, and H. Youk. [Statistical dynamics of spatial-order formation by communicating cells](https://arxiv.org/abs/1706.06481). arXiv: 1706.06481 (2017).

### **For additional information, troubleshooting and a project plan please refer to the [Wiki page](https://github.com/YitengD/Multicellularity/wiki).**

## MATLAB app
### Installation
#### Requirements
* MATLAB (R2017b, recent earlier versions should also work)

#### Instructions
1. Download the [latest release](https://github.com/YitengD/Multicellularity/releases). 
2. Run the MATLAB Installer file 'Multicellularity R1.0.mlappinstall'. The program should unpack itself and install as an app in MATLAB. 
3. Run the app from MATLAB. Go to the tab "Apps" to find the installed app.

#### Advice
Make or select a folder for running your simulation from. In MATLAB, select the folder as your current folder before you start the app and run simulations (and do not change folder while you are still running simulations). Make sure you have the right to write to that folder. MATLAB will create subfolders of the form "./data/dist_matrix_hex" to store data, so that the next time you run a simulation with the same gridsize (see below), the initialisation will be much faster.

### Usage
The GUI consists of a left column with a dial and numerical fields, a graph on the right and buttons on the bottom right. Let's go through the components one by one.

#### Visualisation
The graph on the right becomes a visualisation of the multicellular system once you press on the "RUN" button. The circles represent the cells and the color their expression level. Generally, a darker color corresponds to a higher expression level. For the binary system, OFF cells are white and ON cells are black. In the continuous system, the darker of the shade of grey, the higher the expression level of the cell. 
The caption above the figure displays the time (in units of time steps) and statistics of the current figure. For the binary system, Non is the number of ON cells. For the continuous system, p is the mean expression level of the cells (also see below).

#### Variable conditions
1. **Cell states**. This dial controls whether the cells have a binary response (ON/OFF) or a continuous response (secretion rate continuously varying). The first corresponds to a response curve (secretion rate vs. sensed concentration) that is a step function, whereas the latter corresponds to a Hill function with finite Hill coefficient.
2. **Gridsize**. Sets the size of the lattice. The constructed lattice will be a triangular lattice with both size consisting of 'gridsize' number of cells.
3. **a0**. Effective distance between the cells. This sets the interaction strength between the cells. A small value of a0 corresponds to strongly interacting cells, whereas for large a0 cells have little influence on each other.
4. **K**. Sensing threshold. In the binary system, this sets the threshold for the concentration cells must sense to turn ON. In the continuous system, this is the sensed concentration at which the cells are exactly halfway between their lowest and highest expression levels.
5. **Con**. Maximum secretion rate. This is the secretion rate attained by ON cells in the binary system, or the maximum secretion rate in the continuous system.
6. **Noise**. Stochasticity is modelled as a stochastic term in the sensing threshold K. The value of this field sets the strength of the fluctuations of K. 
7. **Hill coefficient**. (Only for continuous system) sets the steepness of the response function. A high Hill coefficient corresponds to a sharp response around the value of K, whereas a low value corresponds to a gradual response.
8. **Initial p**. Initial mean expression level (values between 0 and 1). In the binary system, this corresponds to the initial fraction of ON cells. 

Detailed descriptions of the variables can be found in the papers listed above.

#### Running the simulation
To run the simulation, simply click on "Run". To terminate a running simulation, click on "Stop". If you click on "Run" again, the app will run a new simulation rather than continue with the stopped simulation. The binary system will terminate when all the cells are in equilibrium (i.e. do not change their state any more). The continuous system will terminate only if you press stop, or when it has progressed 1000 time steps.

## JavaScript app
To be constructed

### Credits
Yiteng Dang

### License
MIT Licence.
