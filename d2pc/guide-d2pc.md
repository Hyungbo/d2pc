## How to test D2PC for your own plant?

0. Download `sim_d2pc.m`.
1. Choose your plant and get its discrete-time matrices (A,B,C).
   - D2PC is not yet able to handle the plant with non-zero D term.
   - If your (A,B,C) is in continuous-time, use `c2d` command of MATLAB to get discrete-time (A,B,C). 
   - This will be the real plant, which will be learned by D2PC and used for simulation.
2. Let the simulator know your plant.
   - Find the function `f_DefinePlant` at the end of the m-file (`sim_d2pc.m`)
   - Give it a name and write a `case` statement, as done for other examples in `f_DefinePlant`.
   - Copy a plant description of other example in `f_DefinePlant', and modify it for your plant.
   - `Ts` is a sampling time (for continuous-time to discrete-time conversion), but it is also used for plotting the outcome by converting the discrete-time steps into continuous-time. Simply set 1 if you don't use it.
   - `An` is the level of noise for output measurement, and the way how to generate the noise is specified in `plant.gen_noise`.
   - Note that `plant.gen_noise` is Anonymous Functions of MATLAB. If you're not familiar with it, consult https://www.mathworks.com/help/matlab/matlab_prog/anonymous-functions.html.
3. Specify your control specification in `sim_d2pc.m`.
   - `OP.r`: Target setpoint to which the plant's output follows.
   - `OP.Q`, `OP.R`, `OP.N`: Weight matrices in the cost function, and the length of optimization horizon
   - `OP.solver`: Set as `OSQP` if you have OSQP solver; `CVX` if you have CVXR (but its performance is not very good); or `handful` if you don't have any solver but are satisfied with an analytic solution that cannot handle your input/output constraints.
4. Set D2PC parameters. (There are only two parameters `NBar` and `Nd`.)
   - Increasing both of them will always given you better results (at the cost of more samples and computation).
   - As a rule of thumb, set twice or triple of your estimated order of the plant.
   - `Nd` specifies how many times your sample episodes are averaged. Typically, 1, or 10, or 50 and so on.
5. Run `sim_d2pc.m`.
6. The simulation outcome briefly displays, but you can draw them by yourself. The outcomes are stored in `plot_time`, `Xd2pc`, and `Ud2pc`.
   - `plot_time`: a row vector contains time tick of `Xd2pc` and `U2pc`
   - `Xd2pc`: each column is a state vector at each time, grows horizontally as time goes.
   - `Ud2pc`: each column is the input vector at each time, grows horizontally as time goes.





