#include "ising1d.h"

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <memory>

                // 1D Ising
// Constructors
Ising1D::Ising1D (int latticeSize, double T, double h, double J) {
    if (latticeSize <= 0 || T <= 0.0) {
        throw std::invalid_argument("lattice size and/or temperature can't be < 0.");
    }
    
    m_latticeSize = latticeSize;
    m_T = T;
    m_h = h;
    m_J = J;
    m_energy = 0.0;
    m_magnetization = 0.0;
    m_datas_magnetization = nullptr;
    m_datas_energy = nullptr;
    m_spins = new int [latticeSize];
}

Ising1D::~Ising1D () {
    if (m_spins != nullptr) {
        delete [] m_spins;
        //std::cout<<"m_spins deleted."<<std::endl; // Just to check
    }
    if (m_datas_magnetization != nullptr) {
        delete [] m_datas_magnetization;
        //std::cout<<"m_datas_magnetization deleted."<<std::endl;
    }
    if (m_datas_energy != nullptr) {
        delete [] m_datas_energy;
        //std::cout<<"m_datas_energy deleted."<<std::endl;
    }
}

// Methods
void Ising1D::init (std::mt19937 & gen) noexcept {
    std::bernoulli_distribution b (0.5); // half a chance to be true
    for (int i {0}; i < m_latticeSize; i++) {
        if (b(gen)) {
            m_spins[i] = 1;
        }
        else {
            m_spins[i] = -1;
        }
    }
}

// Displays the spins in the command line
void Ising1D::show_spins () const noexcept {
    for (int i {0}; i < m_latticeSize; i++) {
        if (m_spins[i] == 1) {
            std::cout<<"+ ";
        }
        if (m_spins[i] == -1) {
            std::cout<<"- ";
        }
    }
    std::cout<<std::endl;
}

void Ising1D::metropolis (int iterations, std::mt19937 & gen) {
    if (iterations < 1) {
        throw std::invalid_argument("number of iterations can't be < 1.");
    }
    
    std::uniform_int_distribution<int> int_dist (0,m_latticeSize-1); // interval [0,L-1]
    std::uniform_real_distribution<double> real_dist (0,1); // interval [0,1[
    
    // keep track of the magnetization and energy
    m_datas_magnetization = new double [iterations];
    m_datas_energy = new double [iterations];
    
    // Apply metropolis scheme at fixed T
    for (int i {0}; i < iterations; i++) {
        if (i == 0) {
            // init magnetization and energy
            m_magnetization = this->calculate_magnetization();
            m_energy = this->calculate_energy();
        }
        int select_random { int_dist(gen) };
        // boundary condition
        double dE {0.0};
        if (select_random != m_latticeSize-1) {
            dE = 2.0*double(m_spins[select_random])*(m_h+m_J*double(m_spins[select_random+1]));
        }
        else {
            dE = 2.0*double(m_spins[select_random])*(m_h+m_J*double(m_spins[0]));
        }
        
        if (dE <= 0) {
            m_spins[select_random] *= -1;
            m_magnetization += 2.0*double(m_spins[select_random]);
        }
        else {
            if (real_dist(gen) < std::exp(-dE/m_T)) {
                m_spins[select_random] *= -1;
                m_magnetization += 2.0*double(m_spins[select_random]);
            }
        }
        m_energy = this->calculate_energy(); // Update energy
        m_datas_magnetization[i] = m_magnetization; // Stock magnetization
        m_datas_energy[i] = m_energy; // Stock energy
    }
}

// Same than above except that M and E are written in a file
void Ising1D::metropolis (int iterations, std::mt19937 & gen, const std::string & filename) {
    if (iterations < 1) {
        throw std::invalid_argument("number of iterations can't be < 1.");
    }
    
    std::ofstream write (filename);
    if (!write) {
        throw std::runtime_error("can't open the file.");
    }
    
    std::uniform_int_distribution<int> int_dist (0,m_latticeSize-1);
    std::uniform_real_distribution<double> real_dist (0,1); // interval [0,1[
    
    // keep track of the magnetization and energy (not needed here)
    m_datas_magnetization = new double [iterations];
    m_datas_energy = new double [iterations];
    
    // Apply metropolis scheme
    for (int i {0}; i < iterations; i++) {
        if (i == 0) {
            // init magnetization and energy
            m_magnetization = this->calculate_magnetization();
            m_energy = this->calculate_energy();
        }
        int select_random { int_dist(gen) };
        // boundary condition
        double dE {0.0};
        if (select_random != m_latticeSize-1) {
            dE = 2.0*double(m_spins[select_random])*(m_h+m_J*double(m_spins[select_random+1]));
        }
        else {
            dE = 2.0*double(m_spins[select_random])*(m_h+m_J*double(m_spins[0]));
        }
        
        if (dE <= 0) {
            m_spins[select_random] *= -1;
            m_magnetization += 2.0*double(m_spins[select_random]);
        }
        else {
            if (real_dist(gen) < std::exp(-dE/m_T)) {
                m_spins[select_random] *= -1;
                m_magnetization += 2.0*double(m_spins[select_random]);
            }
        }
        m_energy = this->calculate_energy();
        m_datas_magnetization[i] = m_magnetization; // not needed here
        m_datas_energy[i] = m_energy; // not needed here
        write<<i+1<<"\t\t"<<std::abs(m_magnetization)/m_latticeSize<<"\t\t"<<m_energy/m_latticeSize<<"\n";
    }
    
    write.close();
}

double Ising1D::calculate_magnetization () const noexcept {
    double mag {0.0};
    for (int i {0}; i < m_latticeSize; i++) {
        mag += double(m_spins[i]);
    }
    return mag;
}

double Ising1D::calculate_energy () const noexcept {
    double energy {0.0};
    for (int i {0}; i < m_latticeSize-1; i++) {
        energy += -m_J*double(m_spins[i]*m_spins[i+1]) - m_h*double(m_spins[i]);
    }
    // boundary condition
    energy += -m_J*double(m_spins[m_latticeSize-1]*m_spins[0]) - m_h*double(m_spins[m_latticeSize-1]);
    return energy;
}

// Getters
double * Ising1D::get_datas_magnetization () const noexcept {
    return m_datas_magnetization;
}

double * Ising1D::get_datas_energy () const noexcept {
    return m_datas_energy;
}

// Non-member functions
/* The functioning of the non-member functions is explained in detail in the report */
void magnetization_vsT_1D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename) {
    if (latticeSize <= 0 || T_max <= 0 || T_min <= 0) {
        throw std::invalid_argument("lattice size and/or temparature can't be <= 0.");
    }
    
    if (T_min >= T_max) {
        throw std::invalid_argument("T_min can't be >= T_max");
    }
    
    if (iterations < 1) {
        throw std::invalid_argument("iterations can't be < 1.");
    }
    
    if (increment <= 0 || increment >= 1) {
        throw std::invalid_argument("increment can't be <= 0 or >= 1.");
    }
    
    if (begin_average <= 0 || begin_average >= iterations) {
        throw std::invalid_argument("begin_average can't be <= 0 or >= iterations.");
    }
    
    std::ofstream write (filename);
    if (!write) {
        throw std::runtime_error("can't open the file.");
    }
    
    double inc {0};
    int max { int((T_max-T_min)/increment) };
    for (int i {0}; i < max; i++) {
        double T { T_min + inc };
        std::unique_ptr<Ising1D> ising1d { new Ising1D (latticeSize, T, h, J) };
        ising1d->init(gen); // init the lattice
        ising1d->metropolis(iterations, gen); // Apply metropolis scheme and then, average (see below)
        
        // Average
        double magnetization {0.0};
        // Uncomment the lines with an * at the end to display the statiscal error in the command line
        //double M_2 {0.0}; *
        for (int i {begin_average}; i < iterations; i++) {
            magnetization += (ising1d->get_datas_magnetization())[i]/(iterations-1-begin_average);
            //M_2 += std::pow((ising1d->get_datas_magnetization())[i], 2)/(iterations-1-begin_average); *
        }
        double L { double(latticeSize) };
        //double M2 { magnetization*magnetization/L }; *
        //M_2 /= L; *
        //double error { std::sqrt((M_2-M2)/(iterations-1-begin_average)) }; *
        //std::cout<<"M("<<T<<") = "<<std::abs(magnetization)/L<<" +- "<<error<<std::endl; *
        
        // magnetization per spin
        write<<T<<"\t\t"<<std::abs(magnetization)/L<<"\n";
        
        inc += increment;
    }
}

void energy_vsT_1D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename) {
    if (latticeSize <= 0 || T_max <= 0 || T_min <= 0) {
        throw std::invalid_argument("lattice size and/or temparature can't be <= 0.");
    }
    
    if (T_min >= T_max) {
        throw std::invalid_argument("T_min can't be >= T_max");
    }
    
    if (iterations < 1) {
        throw std::invalid_argument("iterations can't be < 1.");
    }
    
    if (increment <= 0 || increment >= 1) {
        throw std::invalid_argument("increment can't be <= 0 or >= 1.");
    }
    
    if (begin_average <= 0 || begin_average >= iterations) {
        throw std::invalid_argument("begin_average can't be <= 0 or >= iterations.");
    }
    
    std::ofstream write (filename);
    if (!write) {
        throw std::runtime_error("can't open the file.");
    }
    
    double inc {0};
    int max { int((T_max-T_min)/increment) };
    for (int i {0}; i < max; i++) {
        double T { T_min + inc };
        std::unique_ptr<Ising1D> ising1d { new Ising1D (latticeSize, T, h, J) };
        ising1d->init(gen);
        ising1d->metropolis(iterations, gen);
        
        // Average
        double energy {0.0};
        // Uncomment the lines with an * at the end to display the statiscal error in the command line
        //double E_2 {0.0}; *
        for (int i {begin_average}; i < iterations; i++) {
            energy += (ising1d->get_datas_energy())[i]/(iterations-1-begin_average);
            //E_2 += std::pow((ising1d->get_datas_energy())[i], 2)/(iterations-1-begin_average); *
        }
        double L { double(latticeSize) };
        //double E2 { energy*energy/L }; *
        //E_2 /= L; *
        //double error { std::sqrt((E_2-E2)/(iterations-1-begin_average)) }; *
        //std::cout<<"E("<<T<<") = "<<energy/L<<" +- "<<error<<std::endl; *
        
        // energy per spin
        write<<T<<"\t\t"<<energy/L<<"\n";
        
        inc += increment;
    }
}

void specific_heat_vsT_1D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename) {
    if (latticeSize <= 0 || T_max <= 0 || T_min <= 0) {
        throw std::invalid_argument("lattice size and/or temparature can't be <= 0.");
    }
    
    if (T_min >= T_max) {
        throw std::invalid_argument("T_min can't be >= T_max");
    }
    
    if (iterations < 1) {
        throw std::invalid_argument("iterations can't be < 1.");
    }
    
    if (increment <= 0 || increment >= 1) {
        throw std::invalid_argument("increment can't be <= 0 or >= 1.");
    }
    
    if (begin_average <= 0 || begin_average >= iterations) {
        throw std::invalid_argument("begin_average can't be <= 0 or >= iterations.");
    }
    
    std::ofstream write (filename);
    if (!write) {
        throw std::runtime_error("can't open the file.");
    }
    
    double inc {0.0};
    int max { int((T_max-T_min)/increment) };
    for (int i {0}; i < max; i++) {
        double T { T_min + inc };
        std::unique_ptr<Ising1D> ising1d { new Ising1D (latticeSize, T, h, J) };
        ising1d->init(gen);
        ising1d->metropolis(iterations, gen);
        
        // Average
        double energy {0.0};
        double E_2 {0.0};
        for (int i {begin_average}; i < iterations; i++) {
            energy += (ising1d->get_datas_energy())[i]/(iterations-1-begin_average);
            E_2 += std::pow((ising1d->get_datas_energy())[i], 2)/(iterations-1-begin_average);
        }
        double E2 { energy*energy };
        double L { double(latticeSize) };
        // specific heat per spin
        double specific_heat = (E_2 - E2)/(T*T*L); // fluctuation relation used here
        
        write<<T<<"\t\t"<<specific_heat<<"\n";
        
        inc += increment;
    }
}

void susceptibility_vsT_1D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename) {
    if (latticeSize <= 0 || T_max <= 0 || T_min <= 0) {
        throw std::invalid_argument("lattice size and/or temparature can't be <= 0.");
    }
    
    if (T_min >= T_max) {
        throw std::invalid_argument("T_min can't be >= T_max");
    }
    
    if (iterations < 1) {
        throw std::invalid_argument("iterations can't be < 1.");
    }
    
    if (increment <= 0 || increment >= 1) {
        throw std::invalid_argument("increment can't be <= 0 or >= 1.");
    }
    
    if (begin_average <= 0 || begin_average >= iterations) {
        throw std::invalid_argument("begin_average can't be <= 0 or >= iterations.");
    }
    
    std::ofstream write (filename);
    if (!write) {
        throw std::runtime_error("can't open the file.");
    }
    
    double inc {0.0};
    int max { int((T_max-T_min)/increment) };
    for (int i {0}; i < max; i++) {
        double T { T_min + inc };
        std::unique_ptr<Ising1D> ising1d { new Ising1D (latticeSize, T, h, J) };
        ising1d->init(gen);
        ising1d->metropolis(iterations, gen);
        
        // Average
        double magnetization {0.0};
        double M_2 {0.0};
        for (int i {begin_average}; i < iterations; i++) {
            magnetization += std::abs((ising1d->get_datas_magnetization())[i])/(iterations-1-begin_average);
            M_2 += std::pow((ising1d->get_datas_magnetization())[i], 2)/(iterations-1-begin_average);
        }
        double M2 { magnetization*magnetization };
        double L { double(latticeSize) };
        // susceptibility per spin
        double susceptibility = (M_2 - M2)/(T*L); // fluctuation relation used here
        
        write<<T<<"\t\t"<<susceptibility<<"\n";
        
        inc += increment;
    }
}
