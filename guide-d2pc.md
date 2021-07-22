## How to test D2PC for your own plant?

1. Choose your plant and get its discrete-time matrices (A,B,C).
   - D2PC is not yet able to handle the plant with non-zero D term.
   - If your (A,B,C) is in continuous-time, use `c2d` command of MATLAB to get discrete-time (A,B,C). 
   - This will be the real plant, which will be learned by D2PC and used for simulation.
2. Let the simulator know your plant.
   - Find the function `f_DefinePlant` at the end of the m-file (`run_me.m')
   - Give it a name and write a `case` statement, as done for other examples in `f_DefinePlant`.
   - Copy a plant description of other example in `f_DefinePlant', and modify it for your plant.
   - `Ts` is a sampling time (for continuous-time to discrete-time conversion), but it is also used for plotting the outcome by converting the discrete-time steps into continuous-time. Simply set 1 if you don't use it.
   - `An` is the level of noise for output measurement, and the way how to generate the noise is specified in `plant.gen_noise`.
   - Note that `plant.gen_noise` is Anonymous Functions of MATLAB. If you're not familiar with it, consult https://www.mathworks.com/help/matlab/matlab_prog/anonymous-functions.html.
3. 

