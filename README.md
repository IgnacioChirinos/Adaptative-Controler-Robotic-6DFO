# Adaptative-Controler-Robotic-6DF
**Description![](Aspose.Words.7c38fbed-de80-4a5e-8993-30358aeb2c3e.001.png)**

This project is an implementation of an adaptive controller for a 6 Degrees of Freedom (DOF)system. The controller is implemented in MATLABand its performance is evaluated through simulations. The main objective of the project is to show how well the controller performs in controlling the motion of a 6 DOFsystem.

1. First, the user inputs the desired trajectory for the robotic arm as well as the initial position of the arm.![](Aspose.Words.7c38fbed-de80-4a5e-8993-30358aeb2c3e.002.png)
1. The code initializes the state of the arm based on the initial position input by the user.![](Aspose.Words.7c38fbed-de80-4a5e-8993-30358aeb2c3e.003.png)
1. The adaptive controller is initialized with initial values of the gains and parameters.
1. The code enters a loop where it iteratively simulates the motion of the robotic arm using a numerical integration method (e.g., Euler's method).![](Aspose.Words.7c38fbed-de80-4a5e-8993-30358aeb2c3e.004.png)
1. At each iteration, the current state of the arm is used to compute the control input using the adaptive controller.![](Aspose.Words.7c38fbed-de80-4a5e-8993-30358aeb2c3e.005.png)
1. The control input is then applied to the arm, which updates its state.
1. The code continues this loop until the end of the trajectory is reached or until a certain stopping condition is met.![](Aspose.Words.7c38fbed-de80-4a5e-8993-30358aeb2c3e.003.png)
1. Finally,the results of the simulation are plotted to show the performance of the adaptive controller.

**Features![](Aspose.Words.7c38fbed-de80-4a5e-8993-30358aeb2c3e.006.png)**

- Implementation of an adaptive controller for a 6 DOFsystem.
- Evaluation of the controller's performance through simulations.
- Graphical representation of the simulation results.
- User-friendly interface to adjust the controller parameters.

**Requirements![](Aspose.Words.7c38fbed-de80-4a5e-8993-30358aeb2c3e.007.png)**

- MATLAB2019a or later versions.

**Installation![](Aspose.Words.7c38fbed-de80-4a5e-8993-30358aeb2c3e.008.png)**

1. Clone the repository to your local machine.
2. Open MATLABand navigate to the project directory.![](Aspose.Words.7c38fbed-de80-4a5e-8993-30358aeb2c3e.009.png)
2. Run theModelado\_Control\_Adaptativo.m script to launch the controller interface.

**Usage![](Aspose.Words.7c38fbed-de80-4a5e-8993-30358aeb2c3e.010.png)**

1. Launch the Modelado\_Control\_Adaptativo.m script in MATLABto open the controller interface.
1. Adjust the controller parameters as desired.
1. Run the simulation to evaluate the controller's performance.
1. Viewing the simulation results in the graphical user interface.
