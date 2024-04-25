# Characteristic Mapping Method for the Vlasov Poisson Equation

This is the repository that includes the Matlab source code of the characteristic mapping method and spectral method to solve the Vlasov Poisson equation.

If you use the basis of this code please cite our paper:

    @article{KrahYinBergmannNaveSchneider2024,
            title   ={A characteristic mapping method for Vlasov-Poisson with extreme resolution properties},
            author  ={Krah, Philipp and Yin, Xi-Yuan and Bergmann, Julius and Nave, Jean-Christophe and Schneider, Kai},
            journal ={accepted in CiCP},
            year    ={2024}
    }

Arxiv version: https://arxiv.org/abs/2311.09379

To run the simulations please execute:

    + main_CMM.m to start the CMM simulation
    + main_spectral.m to start the spectral simulation

In the main functions, you find the individual parameter files for the different test cases (see also the folder 'params/'):

    + linear Landau Damping: PARAMS_landau_damping
    + non-lin Landau Damping: PARAMS_non_lin_landau_damping
    + two-stream instability: PARAMS_two_stream  

In the parameter-files you can adjust the parameters according to your needs.