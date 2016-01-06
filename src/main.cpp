#include "main_list.h"
#include <time.h>
#include <ctime>


int main(){
    cout << "Welcome to SMART-UQ!" << endl;
    
    // // Accelerated kepler problem
    // main_AKP_i_canonical();
    // main_AKP_ni_sparse();
    // main_AKP_ni();
    // main_AKP_i_chebyshev();

    // // Van der Pol oscillator
    // main_vanderpol_ni();
    // main_vanderpol_ni_sparse();
    // main_vanderpol_i_chebyshev();
    // main_vanderpol_i_canonical();

    // // Lotka-Volterra predator prey    
    // main_lotka_volterra_i_canonical();
    // main_lotka_volterra_ni();
    // main_lotka_volterra_ni_sparse();
    // main_lotka_volterra_i_chebyshev();

    // // Spring-mass n-dof system
    main_spring_mass_ni();
    
    // // Other problems
    // main_collision_avoidance();
    // main_multiphase();

    // // Tests
    // main_test_dct();
    // main_test_canonical_multiplication();
    // main_test_multivariate();
    // main_test_accuracy();
    // main_test_sampling();
    // main_tests();


    return 0;
}