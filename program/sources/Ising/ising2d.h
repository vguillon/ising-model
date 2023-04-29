#ifndef ising2d_h
#define ising2d_h

#include <vector>
#include <random>
#include <string>

class Ising2D {
private:
    int ** m_spins; // value of the spins
    int m_latticeSize; // in 2D the surface of the lattice is latticeSize power 2
    double m_T; // temperature
    double m_h; // magnetic field
    double m_J; // coupling constant
    double m_energy; // Instananeous energy
    double m_magnetization; // Instantaneous magnetization
    double * m_datas_magnetization; // keep track of the magnetization
    double * m_datas_energy; // keep track of the energy
public:
    // Constructors
    Ising2D () = default;
    Ising2D (int latticeSize, double T, double h, double J);
    ~Ising2D ();
    
    // Methods
    void init (std::mt19937 & gen) noexcept; // init the value of the spins
    void init (std::mt19937 & gen, const std::string & filename); // overload in order to write in file
    void show_spins () const noexcept; // Displays the spins in the command line
    double get_neighbor (int i, int j) const noexcept; // Returns the value of the nearest neigbor
    void metropolis (int iterations, std::mt19937 & gen); // Apply Metropolis scheme
    void metropolis (int iterations, std::mt19937 & gen, const std::string & filename, const std::string & filename2); // overload in order to write in file
    void spins_config (int iterations, std::mt19937 & gen, const std::string & filename, const std::string & filename2, const std::vector<std::string> & filenames); // Write magnetization and spins configuration into a file
    double calculate_magnetization () const noexcept;
    double calculate_energy () const noexcept;
    
    // Getters
    double * get_datas_magnetization () const noexcept;
    double * get_datas_energy () const noexcept;
};

// Non-member functions
void magnetization_vsT_2D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename);

void energy_vsT_2D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename);

void specific_heat_vsT_2D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename);

void susceptibility_vsT_2D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename);

void binder_cumulant (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename);

void reduced_binder_cumulant (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename);

#endif /* ising2d_h */
