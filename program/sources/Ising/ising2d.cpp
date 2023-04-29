#include "ising2d.h"

#include <iostream>
#include <fstream>
#include <memory>
#include <stdexcept>

            // 2D Ising
// Constructors
Ising2D::Ising2D (int latticeSize, double T, double h, double J) {
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
    m_spins = new int * [latticeSize];
    for (int i {0}; i < latticeSize; i++) {
        m_spins[i] = new int [latticeSize];
    }
}

Ising2D::~Ising2D () {
    if (m_spins != nullptr) {
        for (int i {0}; i < m_latticeSize; i++) {
            delete [] m_spins[i];
        }
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
void Ising2D::init(std::mt19937 & gen) noexcept {
    std::bernoulli_distribution b (0.5); // half a chance to be true
    for (int i {0}; i < m_latticeSize; i++) {
        for (int j {0}; j < m_latticeSize; j++) {
            if (b(gen)) {
                m_spins[i][j] = 1;
            }
            else {
                m_spins[i][j] = -1;
            }
        }
    }
}

// Same than above but write the spins configuration into a file
void Ising2D::init (std::mt19937 & gen, const std::string & filename) {
    std::ofstream write (filename);
    if (!write) {
        throw std::runtime_error("can't open the file.");
    }
    
    std::bernoulli_distribution b (0.5); // half a chance to be true
    for (int i {0}; i < m_latticeSize; i++) {
        for (int j {0}; j < m_latticeSize; j++) {
            if (b(gen)) {
                m_spins[i][j] = 1;
            }
            else {
                m_spins[i][j] = -1;
            }
            
            if (m_spins[i][j] == - 1) {
                write<<i+1<<"\t\t"<<j+1<<"\t\t"<<"1\n";
            }
            else {
                write<<i+1<<"\t\t"<<j+1<<"\t\t"<<"2\n";
            }
        }
    }
    
    write.close();
}

// Displays the spins in the command line
void Ising2D::show_spins () const noexcept {
    for (int i {0}; i < m_latticeSize; i++) {
        for (int j {0}; j < m_latticeSize; j++) {
            if (m_spins[i][j] == 1) {
                std::cout<<"+ ";
            }
            if (m_spins[i][j] == -1) {
                std::cout<<"- ";
            }
        }
        std::cout<<std::endl;
    }
}

/* Returns the value of the nearest neigbor. The functioning of the method is explained in
 detail in the the report */
double Ising2D::get_neighbor (int i, int j) const noexcept {
    double sum { 0.0 };
    
    // Consider every special cases
    if (i == 0) {
        sum = double(m_spins[i+1][j]) + double(m_spins[m_latticeSize-1][j]);
        if (j == 0) {
            return sum + double(m_spins[i][m_latticeSize-1]) + double(m_spins[i][j+1]);
        }
        else if (j == m_latticeSize - 1) {
            return sum + double(m_spins[i][j-1]) + double(m_spins[i][0]);
        }
        else {
            return sum + double(m_spins[i][j-1]) + double(m_spins[i][j+1]);
        }
    }
    
    if (i == m_latticeSize - 1) {
        sum = double(m_spins[0][j]) + double(m_spins[i-1][j]);
        if (j == 0) {
            return sum + double(m_spins[i][m_latticeSize-1]) + double(m_spins[i][j+1]);
        }
        else if (j == m_latticeSize - 1) {
            return sum + double(m_spins[i][j-1]) + double(m_spins[i][0]);
        }
        else {
            return sum + double(m_spins[i][j-1]) + double(m_spins[i][j+1]);
        }
    }
    
    if (j == 0) {
        return double(m_spins[i+1][j]) + double(m_spins[i-1][j]) + double(m_spins[i][m_latticeSize-1]) + double(m_spins[i][j+1]);
    }
    
    if (j == m_latticeSize - 1) {
        return double(m_spins[i+1][j]) + double(m_spins[i-1][j]) + double(m_spins[i][j-1]) + double(m_spins[i][0]);
    }
    
    // If not in a special case, just return the nearest neighbor
    return double(m_spins[i+1][j]) + double(m_spins[i-1][j]) + double(m_spins[i][j-1]) + double(m_spins[i][j+1]);
}

// Apply Metropolis scheme
void Ising2D::metropolis(int iterations, std::mt19937 & gen) {
    if (iterations < 1) {
        throw std::invalid_argument("number of iterations can't be < 1.");
    }
    
    std::uniform_int_distribution<int> int_dist (0,m_latticeSize-1); // interval [0,L-1]
    std::uniform_real_distribution<double> real_dist (0,1); // interval [0,1[
    
    // keep track of the magnetization and energy
    m_datas_magnetization = new double [iterations];
    m_datas_energy = new double [iterations];
    
    for (int i {0}; i < iterations; i++) {
        if (i == 0) {
            // init magnetization and energy
            m_magnetization = this->calculate_magnetization();
            m_energy = this->calculate_energy();
        }
        
        int random_i { int_dist(gen) };
        int random_j { int_dist(gen) };
        // Boundary conditions
        double dE { 2.0*double(m_spins[random_i][random_j])*(m_h+m_J*(this->get_neighbor(random_i, random_j))) };
        
        if (dE <= 0) {
            m_spins[random_i][random_j] *= -1;
            m_magnetization += 2.0*double(m_spins[random_i][random_j]); // update magnetization
            m_energy += dE; // update energy
        }
        else {
            if (real_dist(gen) < std::exp(-dE/m_T)) {
                m_spins[random_i][random_j] *= -1;
                m_magnetization += 2.0*double(m_spins[random_i][random_j]); // update magnetization
                m_energy += dE; // update energy
            }
        }
        m_datas_magnetization[i] = m_magnetization; // keep track of the magnetization
        m_datas_energy[i] = m_energy; // keep track of the energy
    }
}

// Same than above except that the magnetization, the energy as well as the spins configuration is written into a file
void Ising2D::metropolis (int iterations, std::mt19937 & gen, const std::string & filename, const std::string & filename2) {
    if (iterations < 1) {
        throw std::invalid_argument("number of iterations can't be < 1.");
    }
    
    std::ofstream write (filename);
    if (!write) {
        throw std::runtime_error("can't open the file.");
    }
    
    std::uniform_int_distribution<int> int_dist (0,m_latticeSize-1); // interval [0,L-1]
    std::uniform_real_distribution<double> real_dist (0,1); // interval [0,1[
    
    // keep track of the magnetization and energy
    m_datas_magnetization = new double [iterations];
    m_datas_energy = new double [iterations];
    
    for (int i {0}; i < iterations; i++) {
        if (i == 0) {
            // init magnetization and energy
            m_magnetization = this->calculate_magnetization();
            m_energy = this->calculate_energy();
        }
        int random_i { int_dist(gen) };
        int random_j { int_dist(gen) };
        double dE { 2.0*double(m_spins[random_i][random_j])*(m_h+m_J*(this->get_neighbor(random_i, random_j))) };
        
        if (dE <= 0) {
            m_spins[random_i][random_j] *= -1;
            m_magnetization += 2.0*double(m_spins[random_i][random_j]);
            m_energy += dE;
        }
        else {
            if (real_dist(gen) < std::exp(-dE/m_T)) {
                m_spins[random_i][random_j] *= -1;
                m_magnetization += 2.0*double(m_spins[random_i][random_j]);
                m_energy += dE;
            }
        }
        m_datas_magnetization[i] = m_magnetization;
        m_datas_energy[i] = m_energy;
        write<<i+1<<"\t\t"<<std::abs(m_magnetization)/(m_latticeSize*m_latticeSize)<<"\t\t"<<m_energy/(m_latticeSize*m_latticeSize)<<"\n";
    }
    
    write.close();
    
    // Write spins configuration into a file at the end of the process
    std::ofstream write2 (filename2);
    if (!write2) {
        throw std::runtime_error("can't open the file.");
    }
    
    for (int i {0}; i < m_latticeSize; i++) {
        for (int j {0}; j < m_latticeSize; j++) {
            if (m_spins[i][j] == - 1) {
                write2<<i+1<<"\t\t"<<j+1<<"\t\t"<<"1\n";
            }
            else {
                write2<<i+1<<"\t\t"<<j+1<<"\t\t"<<"2\n";
            }
        }
    }
    
    write2.close();
}

// Same than above except than 6 spins configurations are written into a file
void Ising2D::spins_config (int iterations, std::mt19937 & gen, const std::string & filename, const std::string & filename2, const std::vector<std::string> & filenames) {
    if (iterations < 1) {
        throw std::invalid_argument("number of iterations can't be < 1.");
    }
    
    std::ofstream write (filename); // datas thermalization
    if (!write) {
        throw std::runtime_error("can't open the file.");
    }
    std::ofstream write2 (filename2); // datas thermalization for just 6 points
    if (!write2) {
        throw std::runtime_error("can't open the file.");
    }
    
    // Spins config for 6 points
    std::ofstream writes1 (filenames[0]); // writeS -> the s stands for Spins
    if (!writes1) {
        throw std::runtime_error("can't open the file.");
    }
    std::ofstream writes2 (filenames[1]);
    if (!writes2) {
        throw std::runtime_error("can't open the file.");
    }
    std::ofstream writes3 (filenames[2]);
    if (!writes3) {
        throw std::runtime_error("can't open the file.");
    }
    std::ofstream writes4 (filenames[3]);
    if (!writes4) {
        throw std::runtime_error("can't open the file.");
    }
    std::ofstream writes5 (filenames[4]);
    if (!writes5) {
        throw std::runtime_error("can't open the file.");
    }
    std::ofstream writes6 (filenames[5]);
    if (!writes6) {
        throw std::runtime_error("can't open the file.");
    }
    
    std::uniform_int_distribution<int> int_dist (0,m_latticeSize-1); // interval [0,L-1]
    std::uniform_real_distribution<double> real_dist (0,1); // interval [0,1[
    
    for (int i {0}; i < iterations; i++) {
        if (i == 0) {
            // init magnetization and energy
            m_magnetization = this->calculate_magnetization();
            m_energy = this->calculate_energy();
        }
        int random_i { int_dist(gen) };
        int random_j { int_dist(gen) };
        double dE { 2.0*double(m_spins[random_i][random_j])*(m_h+m_J*(this->get_neighbor(random_i, random_j))) };
        
        if (dE <= 0) {
            m_spins[random_i][random_j] *= -1;
            m_magnetization += 2.0*double(m_spins[random_i][random_j]);
            m_energy += dE;
        }
        else {
            if (real_dist(gen) < std::exp(-dE/m_T)) {
                m_spins[random_i][random_j] *= -1;
                m_magnetization += 2.0*double(m_spins[random_i][random_j]);
                m_energy += dE;
            }
        }
        write<<i+1<<"\t\t"<<std::abs(m_magnetization)/(m_latticeSize*m_latticeSize)<<"\t\t"<<m_energy/(m_latticeSize*m_latticeSize)<<"\n";
        
        // The lines below write the spins configurations into files
        if (i == 0) {
            write2<<i+1<<"\t\t"<<std::abs(m_magnetization)/(m_latticeSize*m_latticeSize)<<"\t\t"<<m_energy/(m_latticeSize*m_latticeSize)<<"\n";
            for (int j {0}; j < m_latticeSize; j++) {
                for (int k {0}; k < m_latticeSize; k++) {
                    if (m_spins[j][k] == - 1) {
                        writes1<<j+1<<"\t\t"<<k+1<<"\t\t"<<"1\n";
                    }
                    else {
                        writes1<<j+1<<"\t\t"<<k+1<<"\t\t"<<"2\n";
                    }
                }
            }
            writes1.close();
        }
        if (i == int(iterations/1000)) {
            write2<<i+1<<"\t\t"<<std::abs(m_magnetization)/(m_latticeSize*m_latticeSize)<<"\t\t"<<m_energy/(m_latticeSize*m_latticeSize)<<"\n";
            for (int j {0}; j < m_latticeSize; j++) {
                for (int k {0}; k < m_latticeSize; k++) {
                    if (m_spins[j][k] == - 1) {
                        writes2<<j+1<<"\t\t"<<k+1<<"\t\t"<<"1\n";
                    }
                    else {
                        writes2<<j+1<<"\t\t"<<k+1<<"\t\t"<<"2\n";
                    }
                }
            }
            writes2.close();
        }
        if (i == int(3*iterations/100)) {
            write2<<i+1<<"\t\t"<<std::abs(m_magnetization)/(m_latticeSize*m_latticeSize)<<"\t\t"<<m_energy/(m_latticeSize*m_latticeSize)<<"\n";
            for (int j {0}; j < m_latticeSize; j++) {
                for (int k {0}; k < m_latticeSize; k++) {
                    if (m_spins[j][k] == - 1) {
                        writes3<<j+1<<"\t\t"<<k+1<<"\t\t"<<"1\n";
                    }
                    else {
                        writes3<<j+1<<"\t\t"<<k+1<<"\t\t"<<"2\n";
                    }
                }
            }
            writes3.close();
        }
        if (i == int(3*iterations/10)) {
            write2<<i+1<<"\t\t"<<std::abs(m_magnetization)/(m_latticeSize*m_latticeSize)<<"\t\t"<<m_energy/(m_latticeSize*m_latticeSize)<<"\n";
            for (int j {0}; j < m_latticeSize; j++) {
                for (int k {0}; k < m_latticeSize; k++) {
                    if (m_spins[j][k] == - 1) {
                        writes4<<j+1<<"\t\t"<<k+1<<"\t\t"<<"1\n";
                    }
                    else {
                        writes4<<j+1<<"\t\t"<<k+1<<"\t\t"<<"2\n";
                    }
                }
            }
            writes4.close();
        }
        if (i == int(iterations/2)) {
            write2<<i+1<<"\t\t"<<std::abs(m_magnetization)/(m_latticeSize*m_latticeSize)<<"\t\t"<<m_energy/(m_latticeSize*m_latticeSize)<<"\n";
            for (int j {0}; j < m_latticeSize; j++) {
                for (int k {0}; k < m_latticeSize; k++) {
                    if (m_spins[j][k] == - 1) {
                        writes5<<j+1<<"\t\t"<<k+1<<"\t\t"<<"1\n";
                    }
                    else {
                        writes5<<j+1<<"\t\t"<<k+1<<"\t\t"<<"2\n";
                    }
                }
            }
            writes5.close();
        }
        if (i == iterations - 1) {
            write2<<i+1<<"\t\t"<<std::abs(m_magnetization)/(m_latticeSize*m_latticeSize)<<"\t\t"<<m_energy/(m_latticeSize*m_latticeSize)<<"\n";
            for (int j {0}; j < m_latticeSize; j++) {
                for (int k {0}; k < m_latticeSize; k++) {
                    if (m_spins[j][k] == - 1) {
                        writes6<<j+1<<"\t\t"<<k+1<<"\t\t"<<"1\n";
                    }
                    else {
                        writes6<<j+1<<"\t\t"<<k+1<<"\t\t"<<"2\n";
                    }
                }
            }
            writes6.close();
        }
    }
    
    write.close();
}

double Ising2D::calculate_magnetization () const noexcept {
    double mag {0.0};
    for (int i {0}; i < m_latticeSize; i++) {
        for (int j {0}; j < m_latticeSize; j++) {
            mag += double(m_spins[i][j]);
        }
    }
    return mag;
}

double Ising2D::calculate_energy () const noexcept {
    double energy {0.0};
    for (int i {0}; i < m_latticeSize; i++) {
        for (int j {0}; j < m_latticeSize; j++) {
            energy += -m_J*double(m_spins[i][j])*(this->get_neighbor(i, j)) - m_h*double(m_spins[i][j]);
        }
    }
    return energy;
}

// Getters
double * Ising2D::get_datas_magnetization () const noexcept {
    return m_datas_magnetization;
}

double * Ising2D::get_datas_energy () const noexcept {
    return m_datas_energy;
}

// Non-member functions
/* Each of the above functions work the same. The detail of implementation are explained in
 detail in the report */
void magnetization_vsT_2D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename) {
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
        std::unique_ptr<Ising2D> ising2d { new Ising2D (latticeSize, T, h, J) };
        ising2d->init(gen);
        ising2d->metropolis(iterations, gen);
        
        // Average
        double magnetization {0.0};
        for (int i {begin_average}; i < iterations; i++) {
            magnetization += (ising2d->get_datas_magnetization())[i]/(iterations-1-begin_average);
        }
        double L { double(latticeSize*latticeSize) };
        // magnetization per spin
        write<<T<<"\t\t"<<std::abs(magnetization)/L<<"\n";
        
        inc += increment;
    }
}

void energy_vsT_2D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename) {
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
        std::unique_ptr<Ising2D> ising2d { new Ising2D (latticeSize, T, h, J) };
        ising2d->init(gen);
        ising2d->metropolis(iterations, gen);
        
        // Average
        double energy {0.0};
        for (int i {begin_average}; i < iterations; i++) {
            energy += (ising2d->get_datas_energy())[i]/(iterations-1-begin_average);
        }
        double L { double(latticeSize*latticeSize) };
        // energy per spin
        write<<T<<"\t\t"<<energy/L<<"\n";
        
        inc += increment;
    }
}

void specific_heat_vsT_2D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename) {
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
        std::unique_ptr<Ising2D> ising2d { new Ising2D (latticeSize, T, h, J) };
        ising2d->init(gen);
        ising2d->metropolis(iterations, gen);
        
        // Average
        double energy {0.0};
        double E_2 {0.0};
        for (int i {begin_average}; i < iterations; i++) {
            energy += (ising2d->get_datas_energy())[i]/(iterations-1-begin_average);
            E_2 += std::pow((ising2d->get_datas_energy())[i], 2)/(iterations-1-begin_average);
        }
        double E2 { energy*energy };
        double L { double(latticeSize*latticeSize) };
        // specific heat per spin
        double specific_heat = (E_2 - E2)/(T*T*L); // fluctuation relation used here
        
        write<<T<<"\t\t"<<specific_heat<<"\n";
        
        inc += increment;
    }
}

void susceptibility_vsT_2D (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename) {
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
        std::unique_ptr<Ising2D> ising2d { new Ising2D (latticeSize, T, h, J) };
        ising2d->init(gen);
        ising2d->metropolis(iterations, gen);
        
        // Average
        double magnetization {0.0};
        double M_2 {0.0};
        for (int i {begin_average}; i < iterations; i++) {
            magnetization += std::abs((ising2d->get_datas_magnetization())[i])/(iterations-1-begin_average);
            M_2 += std::pow((ising2d->get_datas_magnetization())[i], 2)/(iterations-1-begin_average);
        }
        double M2 { magnetization*magnetization };
        double L { double(latticeSize*latticeSize) };
        // susceptibility per spin
        double susceptibility = (M_2 - M2)/(T*L); // fluctuation relation used here
        
        write<<T<<"\t\t"<<susceptibility<<"\n";
        
        inc += increment;
    }
}

void binder_cumulant (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename) {
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
        std::unique_ptr<Ising2D> ising2d { new Ising2D (latticeSize, T, h, J) };
        ising2d->init(gen);
        ising2d->metropolis(iterations, gen);
        
        // Average
        double M_2 {0.0};
        double M_4 {0.0};
        for (int i {begin_average}; i < iterations; i++) {
            M_2 += std::pow((ising2d->get_datas_magnetization())[i], 2)/(iterations-1-begin_average);
            M_4 += std::pow((ising2d->get_datas_magnetization())[i], 4)/(iterations-1-begin_average);
        }
        double L { double(latticeSize*latticeSize) };
        double M2 { M_2*M_2 };
        M2 /= L;
        M_4 /= L;
        // binder cumulant per spin
        double binder { 1.0 - M_4/(3.0*M2) };
        
        write<<T<<"\t\t"<<binder<<"\n";
        
        inc += increment;
    }
}

// Same as the Binder Cumulant but for the energy (usually used to observe first order phase transition)
void reduced_binder_cumulant (int iterations, std::mt19937 & gen, int latticeSize, double T_min, double T_max, double h, double J, double increment, int begin_average, const std::string & filename) {
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
        std::unique_ptr<Ising2D> ising2d { new Ising2D (latticeSize, T, h, J) };
        ising2d->init(gen);
        ising2d->metropolis(iterations, gen);
        
        // Average
        double E_2 {0.0};
        double E_4 {0.0};
        for (int i {begin_average}; i < iterations; i++) {
            E_2 += std::pow((ising2d->get_datas_energy())[i], 2)/(iterations-1-begin_average);
            E_4 += std::pow((ising2d->get_datas_energy())[i], 4)/(iterations-1-begin_average);
        }
        double L { double(latticeSize*latticeSize) };
        double E2 { E_2*E_2 };
        E2 /= L;
        E_4 /= L;
        // reduced binder cumulant per spin
        double reduced_binder { 1.0 - E_4/(3.0*E2) };
        
        write<<T<<"\t\t"<<reduced_binder<<"\n";
        
        inc += increment;
    }
}


