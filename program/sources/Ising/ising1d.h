#ifndef ising1d_h
#define ising1d_h

#include <random>
#include <string>

class Ising1D {
private:
    int * m_spins; // value of the spins
    int m_latticeSize;
    double m_T; // temperature
    double m_h; // magnetic field
    double m_J; // coupling constant
    double m_energy; // Instantaneous energy
    double m_magnetization; // Instantaneous magnetization
    double * m_datas_magnetization; // keep track of the magnetization
    double * m_datas_energy; // keep track of the energy
public:
    // Constructors
    Ising1D () = default;
    Ising1D (int latticeSize, double T, double h, double J);
    ~Ising1D ();
    
    // Methods
    void init (std::mt19937 & gen) noexcept; // init the value of the spins
    void show_spins () const noexcept; // Displays the spin in the command line
    void metropolis (int iterations, std::mt19937 & gen); // Apply metropolis scheme
    void metropolis (int iterations, std::mt19937 & gen, const std::string & filename); // overload in order to write in file
    double calculate_magnetization () const noexcept;
    double calculate_energy () const noexcept;
    
    // Getters
    double * get_datas_magnetization () const noexcept;
    double * get_datas_energy () const noexcept;
};

// Non member function
void magnetization_vsT_1D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename);

void energy_vsT_1D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename);

void specific_heat_vsT_1D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename);

void susceptibility_vsT_1D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename);

#endif /* ising1d_h */
