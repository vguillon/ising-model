#include <iostream>
#include <random>

#include "ising1d.h"
#include "ising2d.h"

int main(int argc, const char * argv[]) {
    
    /*
        THE MAIN PROGRAM BEGIN IN THE try-catch ENVIRONMENT BELOW
    */
    
    std::random_device rd;
    std::mt19937 gen(rd()); // Seeded with a true random number
    
    /*
        Path to datas files
    */
    /* 1D */
    const std::string datasThermalization_1D ("../../datas/Ising1D/datasThermalization/datasThermalization.txt");
    const std::string datasMagnetizationVsT_1D ("../../datas/Ising1D/datasMagnetizationVsT/datasMagnetizationVsT.txt");
    const std::string datasEnergyVsT_1D ("../../datas/Ising1D/datasEnergyVsT/datasEnergyVsT.txt");
    const std::string datasSusceptibility_1D ("../../datas/Ising1D/datasSusceptibility/datasSusceptibility.txt");
    const std::string datasSpecificHeat_1D ("../../datas/Ising1D/datasSpecificHeat/datasSpecificHeat.txt");
    
    /* 2D */
    const std::string datasInitialization_2D ("../../datas/Ising2D/datasInitialization/datasInitialization.txt");
    const std::string datasThermalization_M_E_2D ("../../datas/Ising2D/datasThermalization/datasThermalization.txt");
    const std::string datasThermalization_spinsConfig_2D ("../../datas/Ising2D/datasThermalization/datasSpinsConfig.txt");
    const std::string datasMagnetizationVsT_2D ("../../datas/Ising2D/datasMagnetizationVsT/datasMagnetizationVsT.txt");
    const std::string datasEnergyVsT_2D ("../../datas/Ising2D/datasEnergyVsT/datasEnergyVsT.txt");
    const std::string datasSpecificHeatVsT_2D ("../../datas/Ising2D/datasSpecificHeatVsT/datasSpecificHeatVsT.txt");
    const std::string datasSusceptibilityVsT_2D ("../../datas/Ising2D/datasSusceptibilityVsT/datasSusceptibilityVsT.txt");
    const std::string datasBinderCumulant ("../../datas/Ising2D/datasBinderCumulant/datasBinderCumulant.txt");
    const std::string datasReducedBinderCumulant ("../../datas/Ising2D/datasReducedBinderCumulant/reducedBinderCumulant.txt");
    
    /* For spins config */
    const std::string datasThermalizationSpinsConfig_2D ("../../datas/Ising2D/datasSpinsConfig/thermalizationSpinsConfig.txt");
    const std::string datasThermalizationSpinsConfig5Points ("../../datas/Ising2D/datasSpinsConfig/thermalizationSpinsConfig5Points.txt");
    const std::string datasSpinsConfig1 ("../../datas/Ising2D/datasSpinsConfig/SpinsConfig1.txt");
    const std::string datasSpinsConfig2 ("../../datas/Ising2D/datasSpinsConfig/SpinsConfig2.txt");
    const std::string datasSpinsConfig3 ("../../datas/Ising2D/datasSpinsConfig/SpinsConfig3.txt");
    const std::string datasSpinsConfig4 ("../../datas/Ising2D/datasSpinsConfig/SpinsConfig4.txt");
    const std::string datasSpinsConfig5 ("../../datas/Ising2D/datasSpinsConfig/SpinsConfig5.txt");
    const std::string datasSpinsConfig6 ("../../datas/Ising2D/datasSpinsConfig/SpinsConfig6.txt");
    
    std::vector<std::string> datasSpinsConfig {""};
    datasSpinsConfig.push_back(datasSpinsConfig1);
    datasSpinsConfig.push_back(datasSpinsConfig2);
    datasSpinsConfig.push_back(datasSpinsConfig3);
    datasSpinsConfig.push_back(datasSpinsConfig4);
    datasSpinsConfig.push_back(datasSpinsConfig5);
    datasSpinsConfig.push_back(datasSpinsConfig6);
    
    // Parameters of every non-member functions:
    // iterations, gen, L, T_min, T_max, h, J, increment, begin_average, datasFile
    try {
        /* 1D */
        //Ising1D ising1d (50, 0.2, 0.0, 1.0); // Instantiate,
        //ising1d.init(gen); // initialize,
        //ising1d.show_spins(); // Display spins configurationa after initialization,
        //ising1d.metropolis(50, gen); // run the Metropolis algorithm,
        //std::cout<<std::endl;
        //ising1d.show_spins(); // Display spins configurations after Metropolis
        //ising1d.metropolis(1000, gen, datasThermalization_1D); // To write the Metropolis datas in a file
        
        // Uncomment one the functions and run the program to write the datas into a file
        
        //magnetization_vsT_1D(150000, gen, 500, 0.01, 5.0, 1.0, 1.0, 0.05, 20000, datasMagnetizationVsT_1D);
        //energy_vsT_1D(150000, gen, 500, 0.01, 5.0, 0.0, 1.0, 0.05, 20000, datasEnergyVsT_1D);
        //susceptibility_vsT_1D(4000000, gen, 500, 0.01, 5.0, 1.0, 1.0, 0.05, 100000, datasSusceptibility_1D);
        //specific_heat_vsT_1D(300000, gen, 500, 0.01, 5.0, 0.0, 1.0, 0.05, 20000, datasSpecificHeat_1D);
        
        /* 2D */
        //Ising2D ising2d (128, 2.269, 0.0, 1.0); // Instantiate
        //ising2d.init(gen, datasInitialization_2D); // Initialize and write the datas into a file or,
        //ising2d.init(gen); // just initialize,
        //ising2d.show_spins(); // Dispaly spins configurations after initialization,
        //ising2d.metropolis(1000, gen); // run the Metropolis algorithm,
        //ising2d.metropolis(5, gen, datasThermalization_M_E_2D, datasThermalization_spinsConfig_2D); // Write the datas into a file
        //std::cout<<std::endl;
        //ising2d.show_spins(); // Display spins configurations after Metropolis
        //ising2d.spins_config(40000000, gen, datasThermalizationSpinsConfig_2D, datasThermalizationSpinsConfig5Points, datasSpinsConfig); // Write the spins configuration into a file
        
        // Uncomment one the functions and run the program to write the datas into a file
        
        //magnetization_vsT_2D(3000000, gen, 16, 0.01, 5.0, 0.0, 1.0, 0.05, 1000000, datasMagnetizationVsT_2D);
        //energy_vsT_2D(3000000, gen, 16, 0.01, 5.0, 0.0, 1.0, 0.05, 1000000, datasEnergyVsT_2D);
        //susceptibility_vsT_2D(5000000, gen, 16, 1.0, 4.0, 0.0, 1.0, 0.05, 2000000, datasSusceptibilityVsT_2D);
        //specific_heat_vsT_2D(5000000, gen, 16, 1.0, 4.0, 0.0, 1.0, 0.05, 2000000, datasSpecificHeatVsT_2D);
        //binder_cumulant(3000000, gen, 4, 2.0, 2.5, 0.0, 1.0, 0.005, 1000000, datasBinderCumulant);
        //reduced_binder_cumulant(50000000, gen, 16, 2.0, 3.0, 0.0, 1.0, 0.004, 20000000, datasReducedBinderCumulant);
    }
    catch (std::invalid_argument const & arg) {
        std::cout<<arg.what()<<std::endl;
    }
    catch (std::runtime_error const & run) {
        std::cout<<run.what()<<std::endl;
    }
    catch (std::exception const & e) {
        std::cout<<e.what()<<std::endl;
    }
    
    return 0;
}
