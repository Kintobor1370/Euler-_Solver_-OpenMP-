#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <utility> 
#include <fstream>
#include <cassert>
#include <tuple>
#include <algorithm>
#include <thread>
#include <time.h>
#include <windows.h>
#include <conio.h>
#include <omp.h>

#include "types.h"

using namespace std;

#define KEY_UP 72
#define KEY_DOWN 80
#define KEY_ENT 13
#define KEY_ESC 27

const CalcType PI = 3.1415927;

enum KernelType { Normal, };

template <class T> void sgn(T &x)
{
    if (x > 0) x = 1;
    if (x < 0) x = -1;
}


// Make time from clock_t
class Time
{
    string str_seconds;
    string str_minutes;
    string str_hours;

public:
    Time(clock_t time)
    {
        int seconds = (int) time / CLOCKS_PER_SEC;
        int hours = seconds / 3600;
        int minutes = (seconds % 3600) / 60;
        seconds = seconds % 60;

        str_seconds = seconds / 10 == 0 ? "0" + to_string(seconds) : to_string(seconds);
        str_minutes = minutes / 10 == 0 ? "0" + to_string(minutes) : to_string(minutes);
        str_hours = to_string(hours);
    }

    friend ostream& operator << (ostream& out, const Time &t)
	{
        out << t.str_hours << ':' << t.str_minutes << ':' << t.str_seconds;
		return out;
	}
};


/**************************************************************************************************************************
**                                            LOGISTIC EQUATION IMPLEMENTATION                                           **
***************************************************************************************************************************
**       The code below realizes the logistic equation for the first and second moments on plane (x, y), as well as      **
**                 properties of the spatial model, kernels generation and all the means of interaction                  **
**                 between parameters and results of the equation and further components of the program.                 **
**************************************************************************************************************************/

//_________________________________________________CONVERTION FUNCTIONS____________________________________________________
Vector MatrixOfMatrixToVector(MatrixOfMatrix mat)
{
    Vector result;
    for (int i=0; i<mat.Rows(); i++)
        for (int j=0; j<mat.Cols(); j++)
        {
            Matrix mat_elem = mat(i, j);
            
            for (int k=0; k<mat_elem.Rows(); k++)
                for (int l=0; l<mat_elem.Cols(); l++)
                    result.add(mat_elem[k][l]);
        }
    return result;
}

MatrixOfMatrix VectorToMatrixOfMatrix(Vector vec, int mat_of_mat_size, int mat_size)
{
    MatrixOfMatrix result(mat_of_mat_size, mat_of_mat_size);
    int index = 0;

    for (int i=0; i<result.Rows(); i++)
        for (int j=0; j<result.Cols(); j++)
        {
            Matrix mat_elem(mat_size, mat_size);

            for (int k=0; k<mat_elem.Rows(); k++)
                for (int l=0; l<mat_elem.Cols(); l++)
                {
                    mat_elem[k][l] = vec[index];
                    index++;
                }
            result(i, j) = mat_elem;
        }
    return result;
}


//___________________________________________________MODEL PROPERTIES______________________________________________________
struct Properties
{
    int SpeciesNumber;                                                  // Number of species                
        
    Vector b;                                                           // Birth rate oer unit time in species i  (b)
    Vector d;                                                           // Death rate per unit time in species i  (d)
    Matrix d_prime;                                                     // Linear dependence of death rate of species i on neighbours of species  (d')
        
    VectorOfMatrix M;                                                   // Rate of birth-movement of species i to distance (x,y)  (m(x))
    MatrixOfMatrix W;                                                   // Death-weighting for neighbours j of species i at distance (x,y)  (w(x))

    Matrix WMaxBin;                                                     // Maximum number of bins over which species i d-interacts with species j
    Vector MMaxBin;                                                     // Maximum number of bins over which birth-movement of species i

    Vector SigmaM;
    Matrix SigmaW;

    CalcType CutoffDistance;                                            // The largest distance between two individuals for which we compute the pair density
    Vector MidPointsKernels;                                            // The mid-points of bins that correspond to kernels indices
    CalcType MidPointsKernelsWidth;
    Vector MidPointsMoments;                                            // The mid-points of bins that correspond to second moments indices

    //..............................GENERATE CALCULATION PARAMETERS
    void GenerateParameters(CalcType sigma_w, CalcType sigma_m)
    {
        CalcType cutoff_distance = (sigma_w + sigma_m) * 3;
        CalcType max_radius = std::max(sigma_w, sigma_m) * 3;

        int num_for_moments = 5;
        auto moments_partition = Vector::abs(-Vector::log(-Vector::setLinSpaced(num_for_moments, 0, num_for_moments - 1) / num_for_moments + 1) * 0.01 * num_for_moments);
        while (moments_partition.last() < cutoff_distance - 0.01)
        {
            num_for_moments++;
            moments_partition = Vector::abs(-Vector::log(-Vector::setLinSpaced(num_for_moments, 0, num_for_moments - 1) / num_for_moments + 1) * 0.01 * num_for_moments);
        }
        moments_partition[moments_partition.size() - 1] = cutoff_distance;

        int num_for_kernels = (num_for_moments - 1) * 2 - 1;
        if (std::min(sigma_w, sigma_m) < 0.015)
        {
            num_for_kernels = (num_for_kernels - 2) * 4 - 1;
        }
        cout << "moments: " << num_for_moments << "\nkernels: " << num_for_kernels << "\n\n";
        assert((num_for_kernels % 2 == 1) && (num_for_kernels >= 5));
        auto kernels_partition = Vector::setLinSpaced(num_for_kernels, 0, max_radius);

        CutoffDistance  = moments_partition.last();
        MidPointsMoments = moments_partition;
        MidPointsKernels = kernels_partition;
        MidPointsKernelsWidth = kernels_partition[2] - kernels_partition[0];
    }

    //..............................GENERATE KERNEL
    pair<Matrix, CalcType> GenerateKernel(
        CalcType s,                                                     // Sigma 
        CalcType R,                                                     // Max radius 
        KernelType type
    ) {
        Vector m = MidPointsKernels;
        CalcType width = MidPointsKernelsWidth;

        Matrix result(m.size(), m.size());
        Vector m_square = m * m;
        int max_non_zero_num = 0; 

        // Calculate
        for (int x=0; x<m.size(); x++)
        {
            for (int y=0; y<m.size(); y++)
            {
                const auto r = std::hypot(m[x], m[y]);
                if(r <= R)
                {
                    switch (type)
                    {
                        case Normal:
                            auto norm = s * std::sqrt(2 * PI) * (1 - std::exp((-R * R) / (2 * s * s)));
                            auto weight = std::exp((-r * r) / (2 * s * s));
                            result[x][y] = weight / norm;
                            break;
                    }
                    if (result[x][y] != 0)
                        max_non_zero_num = std::max({max_non_zero_num, x+1, y+1});
                }
                else result[x][y] = 0;
            }
        }
        if (m.size() % 2 == 1 && max_non_zero_num % 2 == 0 && max_non_zero_num != 2) max_non_zero_num++;
        if (max_non_zero_num == 1) max_non_zero_num = 2;

        // Renormalize
        CalcType integral = 0;
        for (int x=1-max_non_zero_num; x<max_non_zero_num-1; x+=2)
            for (int y=1-max_non_zero_num; y<max_non_zero_num-1; y+=2)
                integral += result[abs(x + 1)][abs(y + 1)];
        
        if (integral == 0)
            cout << result << "\n\n";
        assert(integral != 0);
        integral *= width * width;
        result = result / integral;

        return make_pair(result, max_non_zero_num);
    }

    void CheckKernels()
    {
        cout << "    Checking kernels...\n";
        for (int i=0; i<SpeciesNumber; i++)
        {
            for (int j=0; j<SpeciesNumber; j++)
            {
                cout << "      W_" << i + 1 << "_" << j + 1 << ": integral = ";
                
                CalcType integral = 0;
                //for (int x=1-MidPointsKernels.size(); x<MidPointsKernels.size()-1; x+=2)
                    //for (int y=1-MidPointsKernels.size(); y<MidPointsKernels.size()-1; y+=2)
                for (int x=1-WMaxBin[i][j]; x<WMaxBin[i][j]-1; x+=2)
                    for (int y=1-WMaxBin[i][j]; y<WMaxBin[i][j]-1; y+=2)
                        integral += W(i, j)[abs(x + 1)][abs(y + 1)];
                
                integral *= std::pow(MidPointsKernelsWidth, 2);
                cout << integral;
                assert((integral < 1.00000001) && (integral > 0.999999999));
                cout << ".  Cool!!\n";
            }
        }

        for (int i=0; i<SpeciesNumber; i++)
        {
            cout << "      M_" << i + 1 << ":   integral = ";
                
            CalcType integral = 0;
            //for (int x=1-MidPointsKernels.size(); x<MidPointsKernels.size()-1; x+=2)
                //for (int y=1-MidPointsKernels.size(); y<MidPointsKernels.size()-1; y+=2)
            for (int x=1-MMaxBin[i]; x<MMaxBin[i]-1; x+=2)
                for (int y=1-MMaxBin[i]; y<MMaxBin[i]-1; y+=2)
                    integral += M(i)[abs(x + 1)][abs(y + 1)];
                
            integral *= std::pow(MidPointsKernelsWidth, 2);
            cout << integral;
            assert((integral < 1.00000001) && (integral > 0.999999999));
            cout << ".  Cool!!\n";
        }
        cout << "    Check complete!\n\n";
    }

    //..............................GENERATE MODEL PROPERTIES
    static Properties Generate(
        int num,                                                            // Number of species
        Matrix sW,                                                          // Sigma W
        Vector sM,                                                          // Sigma M
        Vector b,                                                           
        Vector d,
        Matrix dd,
        KernelType w_type = Normal,                                         // W kernel type
        KernelType m_type = Normal                                          // M kernel type
    ) {
        Properties props;
        props.SpeciesNumber = num;
        props.SigmaM = sM;
        props.SigmaW = sW;
        props.b = b;
        props.d = d;
        props.d_prime = dd;
        props.GenerateParameters(sW.max_value(), sM.max_value());

        MatrixOfMatrix w(num, num);
        Matrix w_max_bin(num, num);
        
        Matrix max_radius_W;
        if (w_type == Normal)
            max_radius_W = sW * 3;
        
        for (int i=0; i<num; i++)
        {
            for (int j=0; j<num; j++)
            {
                auto kernel_W = props.GenerateKernel(sW[i][j], max_radius_W[i][j], w_type);
                w(i, j) = kernel_W.first;
                w_max_bin[i][j] = kernel_W.second;
            }
        }    

        VectorOfMatrix m(num);
        Vector m_max_bin(num);
        
        Vector max_radius_M;
        if (m_type == Normal)
            max_radius_M = sM * 3;

        for (int i=0; i<num; i++)
        {
            auto kernel_M = props.GenerateKernel(sM[i], max_radius_M[i], m_type);
            m(i) = kernel_M.first;
            m_max_bin[i] = kernel_M.second;
        }

        props.M = m;
        props.W = w;
        props.WMaxBin = w_max_bin;
        props.MMaxBin = m_max_bin;
        return props;
    }

    //..............................GENERATE MODEL PROPERTIES FOR ONE SPECIES
    static Properties GenerateOneSpecies(
        CalcType sW,
        CalcType sM,
        CalcType b,
        CalcType d,
        CalcType dd,
        KernelType w_type,
        KernelType m_type
    ) {
        return Generate(
            1, 
            Matrix(1, 1, sW), 
            Vector(1, sM), 
            Vector(1, b),
            Vector(1, d),
            Matrix(1, 1, dd),
            w_type,
            m_type
        );
    }
};


//_______________________________EQUATIONS FOR THE DYNAMICS OF SINGLET AND PAIR DENSITIES__________________________________
class Equation
{
    Properties props;                                                   // Model properties
    Vector N;                                                           // Vector of current singlet densities
    MatrixOfMatrix C;                                                   // Current pair densities

private:
    // Make C (Matrix of Matrix) from Vector of values
    MatrixOfMatrix MakeC(Vector C_vals)
    {
        assert(C_vals.size() == props.SpeciesNumber * props.SpeciesNumber * props.MidPointsMoments.size() * props.MidPointsMoments.size());

        auto result = VectorToMatrixOfMatrix(C_vals, props.SpeciesNumber, props.MidPointsMoments.size());
        return result;
    }

    // Make Vector of N (Vector) and C (Matrix of Matrix) values
    Vector MakeVector(Vector N, MatrixOfMatrix C)
    {
        Vector result = N;
        Vector C_To_Vector = MatrixOfMatrixToVector(C);
        assert(result.size() == props.SpeciesNumber);
        assert(C_To_Vector.size() == props.SpeciesNumber * props.SpeciesNumber * props.MidPointsMoments.size() * props.MidPointsMoments.size());
        
        result.add(C_To_Vector);
        return result;
    }

    // Calculte integral w(x)C(x)dx
    Matrix CalculateIntegralOfWAndC(MatrixOfMatrix C, Vector N)
    {
        Matrix result(props.SpeciesNumber, props.SpeciesNumber);

        for (int i=0; i<props.SpeciesNumber; i++) 
            for (int j=0; j<props.SpeciesNumber; j++)
                result[i][j] = MakeIntegralWithPairDensity(0, 0, props.WMaxBin[i][j], props.W(i, j), C(i, j), N[i] * N[j]);

        return result;
    }

    // Calculte integral m(x)C(x)dx
    MatrixOfMatrix CalculateIntegralOfMAndC(MatrixOfMatrix C, Vector N)
    {
        MatrixOfMatrix result(props.SpeciesNumber, props.SpeciesNumber);
        
        for (int i=0; i<props.SpeciesNumber; i++)
        {
            for (int j=0; j<props.SpeciesNumber; j++)
            {
                Matrix result_elem(props.MidPointsMoments.size(), props.MidPointsMoments.size());

                for (int x_num=0; x_num<props.MidPointsMoments.size(); x_num++)
                    for (int y_num=0; y_num<props.MidPointsMoments.size(); y_num++)
                        result_elem[x_num][y_num] = MakeIntegralWithPairDensity(x_num, y_num, props.MMaxBin[i], props.M(i), C(i, j), N[i] * N[j]);

                result(i, j) = result_elem;
            }
        }
        return result;
    }

    CalcType MakeIntegralWithPairDensity(int x_shift, int y_shift, int bin, Matrix kernel, Matrix C, CalcType asymptotic_value)
    {
        Vector m = props.MidPointsKernels;
        CalcType integral = 0;
        //for (int i=1-m.size(); i<m.size()-1; i+=2)
        for (int i=1-bin; i<bin-1; i+=2)
        {
            //for (int j=1-m.size(); j<m.size()-1; j+=2)
            for (int j=1-bin; j<bin-1; j+=2)
            {
                assert(i + 1 < bin);
                assert(j + 1 < bin);

                const CalcType sum_x = abs(props.MidPointsMoments[x_shift] + (i >= 0 ? props.MidPointsKernels[i+1] : -props.MidPointsKernels[-i-1]));
                const CalcType sum_y = abs(props.MidPointsMoments[y_shift] + (j >= 0 ? props.MidPointsKernels[j+1] : -props.MidPointsKernels[-j-1]));

                auto pair_density = BilinearInterpolation(C, props.MidPointsMoments, sum_x, sum_y, props.CutoffDistance, asymptotic_value);
                integral += kernel[abs(i + 1)][abs(j + 1)] * pair_density;
            }
        }
        integral *= pow(props.MidPointsKernelsWidth, 2);
        return integral;
    }

    CalcType BilinearInterpolation(Matrix values, Vector mid_points, CalcType x, CalcType y, CalcType max_distance, CalcType asymptotic_value)
    {
        assert(values.Rows() == mid_points.size());
        assert(values.Cols() == mid_points.size());

        CalcType distance = std::hypot(x, y);
        if (distance > max_distance)
            return asymptotic_value;
        
        // Find segments [a, b] and [c, d] so that x will be in [a, b] and y will be in [c, d]
        int i = 1, j = 1;
        while (x > mid_points[i]) i++;
        while (y > mid_points[j]) j++;
        
        CalcType a = mid_points[i-1];
        CalcType b = mid_points[i];
        CalcType c = mid_points[j-1];
        CalcType d = mid_points[j];

        // Compute weights
        CalcType wx = (x - a) / (b - a);
        CalcType wy = (y - c) / (d - c);

        // Compute weighted value
        CalcType value = values[i-1][j-1] * (1 - wx) * (1 - wy) + 
                         values[i][j-1] * wx * (1 - wy) +
                         values[i-1][j] * (1 - wx) * wy +
                         values[i][j] * wx * wy;
        return value;
    }

public:
    Equation(Properties p) : props(p) {}

    // Make Vector of equilibrium values for N (Vector) and C (Matrix of Matrix). Here C=N*N
    Vector GenerateEquilibriumValuesForNAndC(CalcType value)
    {
        int N_size = props.SpeciesNumber;
        int C_size = props.SpeciesNumber * props.SpeciesNumber * props.MidPointsMoments.size() * props.MidPointsMoments.size();
        Vector result(N_size + C_size);
        
        for (int i=0; i<N_size; i++)
            result[i] = value;
        for (int i=N_size; i<result.size(); i++)
            result[i] = value * value;

        return result;
    }

    // Start Calculations
    Vector operator () (CalcType times, Vector vals)
    {
        N = vals.head(props.SpeciesNumber);
        vals.erase_first(props.SpeciesNumber);
        C = MakeC(vals);

        Vector dN(props.SpeciesNumber);                                 // Vector of derivatives of current first moments
        MatrixOfMatrix dC(props.SpeciesNumber, props.SpeciesNumber);    // Derivatives of current second moments

        Matrix WC = CalculateIntegralOfWAndC(C, N);                     // Integral of W and C
        MatrixOfMatrix MC = CalculateIntegralOfMAndC(C, N);             // Integral of M and C

        //.......................................CALCULATE FIRST MOMENT
        dN = N * (props.b - props.d);                                   // Contribution of birth and death, density independent
        Vector d_prime_sum(props.SpeciesNumber, 0);                     // Death contribution, density dependent
        for (int i=0; i<props.SpeciesNumber; i++)
            for (int j=0; j<props.SpeciesNumber; j++)
                d_prime_sum[i] = d_prime_sum[i] + props.d_prime[i][j] * WC[i][j];
        dN = dN - d_prime_sum;

        //.......................................CALCULATE SECOND MOMENT
        for (int i=0; i<props.SpeciesNumber; i++)
            for (int j=0; j<props.SpeciesNumber; j++)
            {
                Matrix dC_elem(props.MidPointsMoments.size(), props.MidPointsMoments.size());

                for (int x_num=0; x_num<props.MidPointsMoments.size(); x_num++)
                    for (int y_num=0; y_num<props.MidPointsMoments.size(); y_num++)
                    {
                        // Birth contribution, density independent
                        CalcType result = props.b[i] * MC(i, j)[x_num][y_num] + props.b[j] * MC(j, i)[x_num][y_num]; 

                        // Birth contribution, Kronecker symbols
                        if (i == j)
                        {
                            result += 2 * props.b[i] * N[i] * BilinearInterpolation(
                                props.M(i),
                                props.MidPointsKernels,
                                props.MidPointsMoments[x_num],
                                props.MidPointsMoments[y_num],
                                props.SigmaM[i] * 3,
                                0
                            );
                        }

                        // Death contribution, density independent
                        result -= props.d[i] * C(i, j)[x_num][y_num];
                        result -= props.d[j] * C(j, i)[x_num][y_num];

                        // Simple closure
                        for (int k=0; k<props.SpeciesNumber; k++)
                        {
                            result -= props.d_prime[i][k] * WC[i][k] * C(i, j)[x_num][y_num] / N[i];
                            result -= props.d_prime[j][k] * WC[j][k] * C(j, i)[x_num][y_num] / N[j];
                        }
                        
                        result -= props.d_prime[i][j] * C(i, j)[x_num][y_num] * BilinearInterpolation(
                            props.W(i, j),
                            props.MidPointsKernels,
                            props.MidPointsMoments[x_num],
                            props.MidPointsMoments[y_num],
                            props.SigmaW[i][j] * 3,
                            0
                        );
                        result -= props.d_prime[j][i] * C(j, i)[x_num][y_num] * BilinearInterpolation(
                            props.W(j, i),
                            props.MidPointsKernels,
                            props.MidPointsMoments[x_num],
                            props.MidPointsMoments[y_num],
                            props.SigmaW[j][i] * 3,
                            0
                        );
                        
                        dC_elem[x_num][y_num] = result;
                    }
                dC(i, j) = dC_elem;
            }
        return MakeVector(dN, dC);
    }
};


/**************************************************************************************************************************
**                                                         SOLVERS                                                       **
***************************************************************************************************************************
**         The code below realizes solvers for the logistic equation using vairous numerical methods, as well as         **
**                 a structure of all calculations' final results that will be used for creating plots.                  **
**************************************************************************************************************************/

//_____________________________________________________SOLVER RESULTS______________________________________________________
struct SolverResults
{
    vector<CalcType> Times;                                             // Time points
    vector<Vector> Values;                                              // Values of the solution at time
    vector<Vector> Derivatives;                                         // Derivatives

    void clear()                                                        // Clear solver results
    {
        Times.clear();
        Values.clear();
        Derivatives.clear();
    }
};


//___________________________________________________EULER METHOD SOLVER___________________________________________________
class EulerSolver
{
    SolverResults solver_results;

private:
    void Add(Vector value, CalcType time)
    {
        solver_results.Values.emplace_back(std::move(value));
        solver_results.Times.emplace_back(time);
    }

    void AddDerivative(Vector derivative)
    { solver_results.Derivatives.emplace_back(std::move(derivative)); }

public:
    SolverResults Solve(Equation func, CalcType time_max, CalcType step, Vector y_0, bool debug_print_enabled = false)
    {
        CalcType t = 0;                                                 // Time points                
        Vector y = y_0;                                                 // Value points
        Vector dy = func(t, y);                                         // Derivative points
        
        solver_results.clear();
        Add(y, t);
        AddDerivative(dy);
        if (debug_print_enabled)
                cout << "   t = " << t << ";   N = " << y.first() << '\n';

        while (t < time_max)
        {
            t += step;
            y = y + dy * step;
            for (auto& value : y)
                if (value < 0)
                    value = 0;
            Add(y, t);
            dy = func(t, y);
            AddDerivative(dy);
            
            if (debug_print_enabled)
                cout << "   t = " << t << ";   N = " << y.first() << '\n';
        }
    
        return solver_results;
    }
};


/**************************************************************************************************************************
**                                                    INPUT AND OUTPUT                                                   **
***************************************************************************************************************************
**          The code below realizes means of retreiving spatial model parameters from input data ".txt" file,          **
**                       saving the calculations' results in a file and visualising them in a plot.                      **
**************************************************************************************************************************/

//_____________________________________________________FILE READER_________________________________________________________
class FileReader
{
    ifstream file;
    string filename;

private:
    CalcType GetValue()
    {
        string result;
        char new_symbol;
        
        do
        { new_symbol = file.get(); }
        while ((!file.eof()) && (new_symbol == '\n'));
        
        while ((!file.eof()) && (new_symbol != '\n') && (new_symbol != ';') && (new_symbol != ','))
        {
            if(new_symbol != ' ')
                result += new_symbol;
            new_symbol = file.get();
        }
        return stod(result);
    }

    int getInt() { return (int) GetValue(); }

    CalcType getCalcType() { return GetValue(); }

    Vector getVector(int size)
    {
        Vector vector_from_file;
        for (int i=0; i<size; i++)
            vector_from_file.add(GetValue());
        
        return vector_from_file;
    }

    Matrix getMatrix(int rows, int cols)
    {
        Matrix matrix_from_file(rows, cols);
        for (int i=0; i<rows; i++)
            for (int j=0; j<cols; j++)
                matrix_from_file[i][j] = GetValue();

        return matrix_from_file;
    }

    void debug_print(Properties props)
    {
        cout << "MODEL PROPERTIES:\n"
             << "   Species number: " << props.SpeciesNumber << '\n'
             << "   Cutoff distance : " << props.CutoffDistance << '\n'
             << "   Sigma W:  " << props.SigmaW
             << "   Sigma M:  " << props.SigmaM
             << "         b:  " << props.b
             << "         d:  " << props.d
             << "         d': " << props.d_prime
             << "\n\n";
    }

public:
    FileReader(string str) : filename(str) { file.open(filename); }

    ~FileReader() { file.close(); }

    Properties GetProperties(bool debug_print_enabled = false)
    {
        int species = getInt();
    
        Matrix sW = getMatrix(species, species);
        Vector sM = getVector(species);
        
        Vector b = getVector(species);
        Vector d = getVector(species);
        Matrix dd = getMatrix(species, species);

        auto properties = Properties::Generate(species, sW, sM, b, d, dd);
        if (debug_print_enabled)
            debug_print(properties);
        properties.CheckKernels();
        return properties;
    }

    Properties GetWithoutSigmaW(CalcType sigmaW, bool debug_print_enabled = false)
    {
        int species = getInt();
    
        Matrix sW = Matrix(species, species, sigmaW);
        Vector sM = getVector(species);
        
        Vector b = getVector(species);
        Vector d = getVector(species);
        Matrix dd = getMatrix(species, species);

        auto properties = Properties::Generate(species, sW, sM, b, d, dd);
        if (debug_print_enabled)
            debug_print(properties);
        return properties;
    }

    Properties GetWithoutSigmaM(CalcType sigmaM, bool debug_print_enabled = false)
    {
        int species = getInt();
    
        Matrix sW = getMatrix(species, species);
        Vector sM = Vector(species, sigmaM);
        
        Vector b = getVector(species);
        Vector d = getVector(species);
        Matrix dd = getMatrix(species, species);

        auto properties = Properties::Generate(species, sW, sM, b, d, dd);
        if (debug_print_enabled)
            debug_print(properties);
        return properties;
    }
};


//_____________________________________________________FILE WRITER_________________________________________________________
class FileWriter
{
    ofstream file;
    string filename;

public:
    FileWriter(string str) : filename(str) { file.open(filename); }

    ~FileWriter() { file.close(); }

    void WriteResult(string xLabel, string yLabel, vector<CalcType> xValues, vector<CalcType> yValues)
    {
        assert(xValues.size() == yValues.size());

        for (int i=0; i<xValues.size(); i++)
        { file << xLabel + " = " << xValues.at(i) << ";   " + yLabel + " = " << yValues.at(i) << '\n'; }
    }

    void WritePlotData(string number, string title, string xLabel, string yLabel, vector<CalcType> xValues, vector<CalcType> yValues, int number_of_plots = 1)
    {
        assert(number_of_plots > 0);
        assert(number_of_plots * xValues.size() == yValues.size());

        file.clear();
        file << number << '\n' << title << '\n' << xLabel << '\n' << yLabel << '\n' << xValues.size() <<'\n';
        for (auto x : xValues)
            file << x << '\n';
        for (auto y : yValues)
            file << y << '\n';
    }
};


//_____________________________________________________EXPERIMENTS_________________________________________________________
class Experiments
{
    Properties Props;
    vector<Properties> PropsVector;
    string reader_path;
    string writer_path;
    string visdata_path;
    bool enable_debug_print;

private:
    static bool SortBy(const pair<CalcType, CalcType> &a, const pair<CalcType, CalcType> &b)
    { return (a.first < b.first); }

    // Comparator for sorting vector of Solver Results
    static bool SortBySigma(const pair<CalcType, SolverResults> &a, const pair<CalcType, SolverResults> &b)
    { return (a.first < b.first); }

    // Calculate first moment
    SolverResults CalculateFirstMoment(CalcType t_max)
    {
        FileReader Reader(reader_path);

        Props = Reader.GetProperties(enable_debug_print = true);
        Equation equation(Props);
        EulerSolver solver;

        cout << "\n==========================CALCULATIONS START==========================\n";
        clock_t Start = clock();
        
        SolverResults result = solver.Solve(equation, t_max, 0.1, equation.GenerateEquilibriumValuesForNAndC(200), enable_debug_print = true);
        
        clock_t End = clock();
        cout << "\n===========================CALCULATIONS END===========================\n";
        
        Time time(End - Start);
        cout << "TIME:  " << time << '\n';

        return result;
    }

    // Calculate first moment with varying sigma (sM or sW)
    vector<SolverResults> CalculateWithVaryingSigma(CalcType t_max, char SigmaType, Vector SigmaVec)
    {
        vector<pair<CalcType, SolverResults>> SigmasAndResults;
        
        for (int i=0; i<SigmaVec.size(); i++)
        {
            FileReader Reader(reader_path);
            cout << "ENTRY " << i << ":\n";
            switch(SigmaType)
            {
                case 'W':
                    PropsVector.push_back(Reader.GetWithoutSigmaW(SigmaVec[i], enable_debug_print = true));
                    break;
                case 'M':
                    PropsVector.push_back(Reader.GetWithoutSigmaM(SigmaVec[i], enable_debug_print = true));
                    break;
                default:
                    cerr << "===============================ERROR!!!===============================\n" << endl;
                    break;
            }
        }

        int done_iterations = 0;
        cout << "\n==========================CALCULATIONS START==========================\n";
        clock_t Start = clock();

        #pragma omp parallel for default(none) shared(t_max, SigmaVec, SigmasAndResults, done_iterations, cout)
        for (int i=0; i<SigmaVec.size(); i++)
        {
            cout << (int) (done_iterations * 100 / SigmaVec.size()) << "%    THREAD " << omp_get_thread_num() << ":     ENTRY " << i << ": START\n";
            
            auto new_props = PropsVector.at(i);
            Equation equation(new_props);
            EulerSolver solver;

            auto new_result = solver.Solve(equation, t_max, 0.1, equation.GenerateEquilibriumValuesForNAndC(200), enable_debug_print = false);
            SigmasAndResults.push_back(make_pair(SigmaVec[i], new_result));
            
            done_iterations++;
            cout << (int) (done_iterations * 100 / SigmaVec.size()) << "%    THREAD " << omp_get_thread_num() << ":     ENTRY " << i << ": FINISH    "
                 << "x = " << SigmaVec[i] << ";  N = " << new_result.Values.at(new_result.Values.size() - 1).first() << '\n';
        }
        clock_t End = clock();
        cout << "\n===========================CALCULATIONS END===========================\n";
        
        Time time(End - Start);
        cout << "TIME:  " << time << '\n';

        sort(SigmasAndResults.begin(), SigmasAndResults.end(), SortBySigma);
        vector<SolverResults> results;
        for (auto result : SigmasAndResults)
            results.push_back(result.second);
        return results;
    }

    void FunctionTest(SolverResults Result, string plot_title)
    {
        assert(Props.SpeciesNumber == 1);

        auto ConvergedNAndC = Result.Values.at(Result.Values.size() - 1);
        CalcType N = ConvergedNAndC[0];
        ConvergedNAndC.erase_first(1);
        Vector C_vals = ConvergedNAndC;
        
        auto C = VectorToMatrixOfMatrix(C_vals, 1, Props.MidPointsMoments.size());
        vector<pair<CalcType, CalcType>> C_r;
        bool unique;

        cout << Props.MidPointsMoments.size() << "\n\n";
        for (int x=0; x<Props.MidPointsMoments.size(); x++)
        {
            for (int y=0; y<Props.MidPointsMoments.size(); y++)
            {
                unique = true;
                CalcType r = std::hypot(Props.MidPointsMoments[x], Props.MidPointsMoments[y]);
                cout << "HERE1\n";
                for (auto item : C_r)
                {
                    if (r == item.first)
                        unique = false;
                }
                cout << "HERE2\n";
                if (unique)
                    C_r.push_back(make_pair(r, C(0, 0)[x][y]));
                cout << "HERE3\n";
            }
        }
        sort(C_r.begin(), C_r.end(), SortBy);

        vector<CalcType> x;
        vector<CalcType> y;
        for (auto item : C_r)
        {
            CalcType func = -std::log(std::abs(item.second - N * N));
            x.push_back(item.first);
            y.push_back(func);
        }

        FileWriter VisData(visdata_path);
        VisData.WritePlotData("14.4", plot_title, "r", "-ln|C(r) - N^2|", x, y, 1);
        Visualise_Windows();
    }

    // Launch Visualiser.py     !! THIS FUNCTION WORKS ONLY ON WINDOWS !!
    void Visualise_Windows() { WinExec("python Visualiser.py", 1); }

public:
    Experiments() : writer_path("Data\\Results.txt"), visdata_path("Data\\VisData.txt") {}

    ~Experiments() { PropsVector.clear(); }

    void Experiment_14_4(int input_data)
    {
        string plot_title;
        switch(input_data)
        {
            case 0:
                reader_path = "Data\\14\\4\\a.txt";
                plot_title = "(a)";
                break;
            case 1:
                reader_path = "Data\\14\\4\\b.txt";
                plot_title = "(b)";
                break;
            case 2:
                reader_path = "Data\\14\\4\\c.txt";
                plot_title = "(c)";
                break;
        }
        PropsVector.clear();
        SolverResults Results = CalculateFirstMoment(100);
    
        Vector N;
        Vector C_sqrt;
        for (auto Values : Results.Values)
        {
            CalcType N_val = Values.first();
            N.add(N_val);
            
            CalcType C_square_root = std::sqrt(Values[Props.MidPointsMoments.size()]);
            C_sqrt.add(C_square_root);
        }
        Vector result = N;
        result.add(C_sqrt);

        FileWriter VisData(visdata_path);
        VisData.WritePlotData("14.4", plot_title, "Time", "N", Results.Times, result.get_vector(), 2);
        Visualise_Windows();
    }

    void Experiment_14_6(char SigmaType, string RadiusType)
    {
        reader_path = "Data\\14\\6\\a.txt";
        string plot_title = "(a)";
        if(SigmaType == 'M')
        {
            reader_path = "Data\\14\\6\\b.txt";
            plot_title = "(b)";
        }
        PropsVector.clear();
        Vector SigmaVec = Vector::setLinSpaced(160, 0.00000001, 0.2);
        auto SolverResults = CalculateWithVaryingSigma(250, SigmaType, SigmaVec);
        
        Vector ConvergedNVec;
        Vector ConvergedCSqrtVec;
        for (auto Result : SolverResults)
        {
            auto ConvergedNAndC = Result.Values.at(Result.Values.size() - 1);
            CalcType ConvergedN = ConvergedNAndC.first();
            ConvergedNVec.add(ConvergedN);

            CalcType ConvergedCSqrt = std::sqrt(ConvergedNAndC[Props.MidPointsMoments.size()]);
            ConvergedCSqrtVec.add(ConvergedCSqrt);
        }
        Vector result = ConvergedNVec;
        result.add(ConvergedCSqrtVec);
        
        FileWriter VisData(visdata_path);
        VisData.WritePlotData("14.6", plot_title, RadiusType, "Equilibrium density", SigmaVec.get_vector(), result.get_vector(), 2);
        Visualise_Windows();
    }

    void Custom_Experiment()
    {
        PropsVector.clear();
        reader_path = "Data\\Custom\\experiment.txt";
        SolverResults Results = CalculateFirstMoment(100);
    
        Vector N;
        Vector C_sqrt;
        for (auto Values : Results.Values)
        {
            CalcType N_val = Values.first();
            N.add(N_val);
            
            CalcType C_square_root = std::sqrt(Values[Props.MidPointsMoments.size()]);
            C_sqrt.add(C_square_root);
        }

        Vector result = N;
        result.add(C_sqrt);

        FileWriter VisData(visdata_path);
        VisData.WritePlotData("14.4", " ", "Time", "N", Results.Times, result.get_vector(), 2);
        Visualise_Windows();

        //FunctionTest(Results, " ");
    }

    void KernelTest()
    {
        FileReader Reader("Data\\14\\4\\a.txt");
        auto Props = Reader.GetProperties(enable_debug_print = true);
        
        auto kernel = Props.W(0, 0);
        //auto kernel = Props.M(0);
    
        vector<pair<CalcType, CalcType>> result;
        bool unique;

        for (int x=0; x<Props.MidPointsKernels.size(); x++)
        {
            for (int y=0; y<Props.MidPointsKernels.size(); y++)
            {
                unique = true;
                CalcType r = std::hypot(Props.MidPointsKernels[x], Props.MidPointsKernels[y]);
                for (auto item : result)
                {
                    cout << "x = " << Props.MidPointsKernels[x] << "    y = " << Props.MidPointsKernels[y] << "     r = " << r << "    kernel = " << kernel[x][y] << '\n';
                    if (r == item.first)
                        unique = false;
                }
                if (unique)
                    result.push_back(make_pair(r, kernel[x][y]));
            }
        }
        sort(result.begin(), result.end(), SortBy);

        vector<CalcType> x;
        vector<CalcType> y;
        for (auto item : result)
        {
            x.push_back(item.first);
            y.push_back(item.second);
        }
        cout << "MAX BIN = " << Props.WMaxBin << '\n';
        //cout << "MAX BIN = " << Props.MMaxBin << '\n';
        FileWriter VisData(visdata_path);
        VisData.WritePlotData("14.4", "Kernel test", "x", "distribution", x, y);
        Visualise_Windows();
    }
};


/**************************************************************************************************************************
**                                                        MAIN BODY                                                      **
**************************************************************************************************************************/
int main()
{
    Experiments Run;
    //Run.KernelTest();
    //Run.DistanceTest();
    string commands[] =
    {
        "  > 14.4.a\n    14.4.b\n    14.4.c\n    14.6.a\n    14.6.b\n    Custom experiment\n",
        "    14.4.a\n  > 14.4.b\n    14.4.c\n    14.6.a\n    14.6.b\n    Custom experiment\n",
        "    14.4.a\n    14.4.b\n  > 14.4.c\n    14.6.a\n    14.6.b\n    Custom experiment\n",
        "    14.4.a\n    14.4.b\n    14.4.c\n  > 14.6.a\n    14.6.b\n    Custom experiment\n",
        "    14.4.a\n    14.4.b\n    14.4.c\n    14.6.a\n  > 14.6.b\n    Custom experiment\n",
        "    14.4.a\n    14.4.b\n    14.4.c\n    14.6.a\n    14.6.b\n  > Custom experiment\n"
    };
    int current_index = 0;

    system("cls");
    cout << "Welcome!\nPlease choose experiment:\n" + commands[current_index] + "\n(Navigate with up arrow & down arrow\n Press ENTER to confirm\n Press ESC to quit)\n";
    char input = getch();

    while(input != KEY_ESC)
    {
        switch (input)
        {
            case KEY_UP:
                current_index = (current_index + 5) % 6;
                system("cls");
                cout << "Welcome!\nPlease choose experiment:\n" + commands[current_index] + "\n(Navigate with up arrow & down arrow.\n Press ENTER to confirm\n Press ESC to quit)\n";
                input = getch();
                break;

            case KEY_DOWN:
                current_index = (current_index + 1) % 6;
                system("cls");
                cout << "Welcome!\nPlease choose experiment:\n" + commands[current_index] + "\n(Navigate with up arrow & down arrow.\n Press ENTER to confirm\n Press ESC to quit)\n";
                input = getch();
                break;

            case KEY_ENT:
                system("cls");
                switch(current_index)
                {
                    case 0: case 1: case 2:
                        Run.Experiment_14_4(current_index);
                        break;
                    case 3:
                        Run.Experiment_14_6('W', "Competition radius");
                        break;
                    case 4:
                        Run.Experiment_14_6('M', "Dispersal radius");
                        break;
                    case 5:
                        cout << "Please navigate to \'Data\\Custom\' folder and fill the \'experiment.txt\' file according to the example.\nPress ENTER to confirm and start the experiment\nPress ESC to quit\n";
                        input = getch();
                        if (input == KEY_ENT)
                        {
                            system("cls");
                            Run.Custom_Experiment();
                        }
                        else if (input != KEY_ESC)
                            input = getch();
                        break;
                }

            default:
                input = getch();
                break;
        }
    }
    system("cls");

    return 0;
}