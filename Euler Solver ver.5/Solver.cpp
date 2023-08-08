#define NOMINMAX

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


/**************************************************************************************************************************
**                                            LOGISTIC EQUATION IMPLEMENTATION                                           **
***************************************************************************************************************************
**       The code below realizes the logistic equation for the first and second moments on plane (x, y), as well as      **
**                 properties of the spatial model, kernels generation and all the means of interaction                  **
**                 between parameters and results of the equation and further components of the program.                 **
**************************************************************************************************************************/

//___________________________________________________MODEL PROPERTIES______________________________________________________
struct Properties
{
    int SpeciesNumber;                                                  // Number of species                
        
    Vector b;                                                           // Birth rate oer unit time in species i  (b)
    Vector d;                                                           // Death rate per unit time in species i  (d)
    Matrix d_prime;                                                     // Linear dependence of death rate of species i on neighbours of species  (d')
        
    VectorOfMatrix M;                                                   // Rate of birth-movement of species i to distance (x,y)  (m(x))
    MatrixOfMatrix W;                                                   // Death-weighting for neighbours j of species i at distance (x,y)  (w(x))

    Matrix MaxNumForW;                                                  // Maximum number of non-zero bins for 'w' kernel
    Vector MaxNumForM;                                                  // Maximum number of non-zero bins for 'm' kernel

    Vector SigmaM;
    Matrix SigmaW;

    Vector MomentsGrid;                                                 // The non-uniform grid that corresponds to pair densities
    CalcType CutoffDistance;                                            // The largest distance between two individuals for which we compute the pair density
    Vector KernelsGrid;                                                 // The uniform grid that corresponds to kernels
    CalcType KernelsGridStep;                                           // The step of uniform grid for kernels

    //..............................GENERATE GRIDS FOR KERNELS AND PAIR DENSITIES
    void GenerateGrids(CalcType sigma_w, CalcType sigma_m)
    {
        CalcType cutoff_distance = (sigma_w + sigma_m) * 5;
        CalcType max_radius = std::max(sigma_w, sigma_m) * 5;

        int num_for_moments = 5;
        auto moments_grid = Vector::abs(-Vector::log(-Vector::setLinSpaced(num_for_moments, 0, num_for_moments - 1) / num_for_moments + 1) * 0.01 * num_for_moments);
        while (moments_grid.last() < cutoff_distance - 0.01)
        {
            num_for_moments++;
            moments_grid = Vector::abs(-Vector::log(-Vector::setLinSpaced(num_for_moments, 0, num_for_moments - 1) / num_for_moments + 1) * 0.01 * num_for_moments);
        }
        CalcType parameter = cutoff_distance / moments_grid.last();
        moments_grid = moments_grid * parameter;
        assert(std::abs(moments_grid.last() - cutoff_distance) < 1e-08);

        int num_for_kernels = (num_for_moments - 1) * 2 - 1;
        if (std::min(sigma_w, sigma_m) < 0.015)
        {
            num_for_kernels = (num_for_kernels - 2) * 4 - 1;
            num_for_kernels = (num_for_kernels - 1) * 2 - 1;
        }
        cout << "number of mid-points for pair density grid: " << num_for_moments
             << "\nnumber of mid-points for kernels grid:      " << num_for_kernels << "\n\n";

        assert ((num_for_kernels % 2 == 1) && (num_for_kernels >= 5));
        auto kernels_grid = Vector::setLinSpaced(num_for_kernels, 0, max_radius);

        CutoffDistance  = moments_grid.last();
        MomentsGrid = moments_grid;
        KernelsGrid = kernels_grid;
        KernelsGridStep = kernels_grid[2] - kernels_grid[0];
    }

    //..............................GENERATE KERNEL
    pair<Matrix, CalcType> GenerateKernel(
        CalcType s,                                                     // Sigma 
        CalcType R,                                                     // Max radius 
        KernelType type
    ) {
        Vector grid = KernelsGrid;
        CalcType step = KernelsGridStep;

        Matrix kernel_values(grid.size());
        int max_non_zero_num = 0; 

        // Calculate
        for (int x=0; x<grid.size(); x++)
        {
            for (int y=0; y<grid.size(); y++)
            {
                const auto r = std::hypot(grid[x], grid[y]);
                if(r <= R)
                {
                    switch (type)
                    {
                        case Normal:
                            auto norm = s * std::sqrt(2 * PI) * (1 - std::exp((-R * R) / (2 * s * s)));
                            auto weight = std::exp((-r * r) / (2 * s * s));
                            kernel_values[x][y] = weight / norm;
                            break;
                    }
                    if (kernel_values[x][y] != 0)
                        max_non_zero_num = std::max({max_non_zero_num, x+1, y+1});
                }
                else kernel_values[x][y] = 0;
            }
        }
        if (grid.size() % 2 == 1 && max_non_zero_num % 2 == 0 && max_non_zero_num != 2) max_non_zero_num++;
        if (max_non_zero_num == 1) max_non_zero_num = 2;

        // Renormalize
        CalcType integral = 0;
        for (int x=1-max_non_zero_num; x<max_non_zero_num-1; x+=2)
            for (int y=1-max_non_zero_num; y<max_non_zero_num-1; y+=2)
                integral += kernel_values[abs(x + 1)][abs(y + 1)];
        
        assert (integral != 0);
        integral *= step * step;
        kernel_values = kernel_values / integral;

        return make_pair(kernel_values, max_non_zero_num);
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
                for (int x=1-MaxNumForW[i][j]; x<MaxNumForW[i][j]-1; x+=2)
                    for (int y=1-MaxNumForW[i][j]; y<MaxNumForW[i][j]-1; y+=2)
                        integral += W(i, j)[abs(x + 1)][abs(y + 1)];
                
                integral *= std::pow(KernelsGridStep, 2);
                cout << integral;
                assert(std::abs(integral - 1) < 1e-08);
                cout << ".  Cool!!\n";
            }
        }

        for (int i=0; i<SpeciesNumber; i++)
        {
            cout << "      M_" << i + 1 << ":   integral = ";
                
            CalcType integral = 0;
            for (int x=1-MaxNumForM[i]; x<MaxNumForM[i]-1; x+=2)
                for (int y=1-MaxNumForM[i]; y<MaxNumForM[i]-1; y+=2)
                    integral += M(i)[abs(x + 1)][abs(y + 1)];
                
            integral *= std::pow(KernelsGridStep, 2);
            cout << integral;
            assert(std::abs(integral - 1) < 1e-08);
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
        props.GenerateGrids(sW.max_value(), sM.max_value());

        MatrixOfMatrix w(num);
        Matrix max_num_w(num);
        
        Matrix max_radius_W;
        if (w_type == Normal)
            max_radius_W = sW * 3;
        
        for (int i=0; i<num; i++)
        {
            for (int j=0; j<num; j++)
            {
                auto kernel_W = props.GenerateKernel(sW[i][j], max_radius_W[i][j], w_type);
                w(i, j) = kernel_W.first;
                max_num_w[i][j] = kernel_W.second;
            }
        }    

        VectorOfMatrix m(num);
        Vector max_num_m(num);
        
        Vector max_radius_M;
        if (m_type == Normal)
            max_radius_M = sM * 3;

        for (int i=0; i<num; i++)
        {
            auto kernel_M = props.GenerateKernel(sM[i], max_radius_M[i], m_type);
            m(i) = kernel_M.first;
            max_num_m[i] = kernel_M.second;
        }

        props.M = m;
        props.W = w;
        props.MaxNumForW = max_num_w;
        props.MaxNumForM = max_num_m;
        return props;
    }
};


//_______________________________EQUATIONS FOR THE DYNAMICS OF SINGLET AND PAIR DENSITIES__________________________________
class Equation
{
    Properties props;                                                   // Model properties
    Vector N;                                                           // Vector of current singlet densities
    MatrixOfMatrix C;                                                   // Current pair densities

private:
    // Make Vector of N (Vector) and C (Matrix of Matrix) values
    Vector MakeVector(Vector N, MatrixOfMatrix C)
    {
        Vector result = N;
        result.add(MatrixOfMatrixToVector(C));

        assert (result.size() == props.SpeciesNumber + props.SpeciesNumber * props.SpeciesNumber * props.MomentsGrid.size() * props.MomentsGrid.size());
        return result;
    }

    // Calculte integral w(x)C(x)dx
    Matrix CalculateIntegralOfWAndC(MatrixOfMatrix C, Vector N)
    {
        Matrix result(props.SpeciesNumber);

        for (int i=0; i<props.SpeciesNumber; i++) 
            for (int j=0; j<props.SpeciesNumber; j++)
                result[i][j] = MakeIntegralWithPairDensity(0, 0, props.MaxNumForW[i][j], props.W(i, j), C(i, j), N[i] * N[j]);

        return result;
    }

    // Calculte integral m(x)C(x)dx
    MatrixOfMatrix CalculateIntegralOfMAndC(MatrixOfMatrix C, Vector N)
    {
        MatrixOfMatrix result(props.SpeciesNumber);
        
        for (int i=0; i<props.SpeciesNumber; i++)
        {
            for (int j=0; j<props.SpeciesNumber; j++)
            {
                Matrix result_elem(props.MomentsGrid.size());

                for (int x_num=0; x_num<props.MomentsGrid.size(); x_num++)
                    for (int y_num=0; y_num<props.MomentsGrid.size(); y_num++)
                        result_elem[x_num][y_num] = MakeIntegralWithPairDensity(x_num, y_num, props.MaxNumForM[i], props.M(i), C(i, j), N[i] * N[j]);

                result(i, j) = result_elem;
            }
        }
        return result;
    }

    CalcType MakeIntegralWithPairDensity(int x_shift, int y_shift, int max_non_zero_num, Matrix kernel, Matrix C, CalcType asymptotic_value)
    {
        CalcType integral = 0;
        for (int i=1-max_non_zero_num; i<max_non_zero_num-1; i+=2)
        {
            for (int j=1-max_non_zero_num; j<max_non_zero_num-1; j+=2)
            {
                assert (i + 1 < max_non_zero_num);
                assert (j + 1 < max_non_zero_num);

                const CalcType sum_x = abs(props.MomentsGrid[x_shift] + (i >= 0 ? props.KernelsGrid[i+1] : -props.KernelsGrid[-i-1]));
                const CalcType sum_y = abs(props.MomentsGrid[y_shift] + (j >= 0 ? props.KernelsGrid[j+1] : -props.KernelsGrid[-j-1]));

                auto pair_density = BilinearInterpolation(props.MomentsGrid, C, sum_x, sum_y, props.CutoffDistance, asymptotic_value);
                integral += kernel[abs(i + 1)][abs(j + 1)] * pair_density;
            }
        }
        integral *= pow(props.KernelsGridStep, 2);
        return integral;
    }

    CalcType BilinearInterpolation(Vector mid_points, Matrix values, CalcType x, CalcType y, CalcType max_distance, CalcType asymptotic_value)
    {
        assert (mid_points.size() == values.size());

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
        int C_size = props.SpeciesNumber * props.SpeciesNumber * props.MomentsGrid.size() * props.MomentsGrid.size();
        
        Vector result;  
        for (int i=0; i<N_size; i++)
            result.add(value);
        for (int i=0; i<C_size; i++)
            result.add(value * value);

        return result;
    }

    // Start Calculations
    Vector operator () (CalcType times, Vector vals)
    {
        N = vals.head(props.SpeciesNumber);
        vals.erase_first(props.SpeciesNumber);
        C = VectorToMatrixOfMatrix(vals, props.SpeciesNumber, props.MomentsGrid.size());

        Vector dN(props.SpeciesNumber);                                 // Vector of derivatives of current first moments
        MatrixOfMatrix dC(props.SpeciesNumber);                         // Derivatives of current second moments

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
                Matrix dC_elem(props.MomentsGrid.size());

                for (int x_num=0; x_num<props.MomentsGrid.size(); x_num++)
                {
                    for (int y_num=0; y_num<props.MomentsGrid.size(); y_num++)
                    {
                        // Birth contribution, density independent
                        CalcType result = props.b[i] * MC(i, j)[x_num][y_num] + props.b[j] * MC(j, i)[x_num][y_num]; 

                        // Birth contribution, Kronecker symbols
                        if (i == j)
                        {
                            result += 2 * props.b[i] * N[i] * BilinearInterpolation(
                                props.KernelsGrid,
                                props.M(i),
                                props.MomentsGrid[x_num],
                                props.MomentsGrid[y_num],
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
                            props.KernelsGrid,
                            props.W(i, j),
                            props.MomentsGrid[x_num],
                            props.MomentsGrid[y_num],
                            props.SigmaW[i][j] * 3,
                            0
                        );
                        result -= props.d_prime[j][i] * C(j, i)[x_num][y_num] * BilinearInterpolation(
                            props.KernelsGrid,
                            props.W(j, i),
                            props.MomentsGrid[x_num],
                            props.MomentsGrid[y_num],
                            props.SigmaW[j][i] * 3,
                            0
                        );
                        
                        dC_elem[x_num][y_num] = result;
                    }
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
    Vector Times;                                                       // Time points
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
        solver_results.Times.add(time);
    }

    void AddDerivative(Vector derivative)
    { solver_results.Derivatives.emplace_back(std::move(derivative)); }

public:
    // Solve given equation numericaly for fixed amount of iterations
    SolverResults Solve(
        Equation current_equation,
        CalcType time_max,
        CalcType step,
        Vector initial_value,
        bool debug_print_enabled = false
    ) {
        CalcType time = 0;
        Vector value = initial_value;
        Vector derivative = current_equation(time, value);

        Add(value, time);
        AddDerivative(derivative);
        if (debug_print_enabled)
                cout << "   t = " << time << ";   N = " << value.first() << '\n';
        
        time += step;
        while (time < time_max)
        {
            value = value + derivative * step;
            for (auto& val : value)
                if (val < 0)
                    val = 0;
            Add(value, time);
            derivative = current_equation(time, value);
            AddDerivative(derivative);
            
            if (debug_print_enabled)
                cout << "   t = " << time << ";   N = " << value.first() << '\n';

            time += step;
        }
    
        return solver_results;
    }
};


//________________________________________________RUNGE-KUTTA 4 METHOD SOLVER_______________________________________________
class RungeKutta4thSolver
{
    SolverResults solver_results;

private:
    void Add(Vector value, CalcType time)
    {
        solver_results.Values.emplace_back(std::move(value));
        solver_results.Times.add(time);
    }

    void AddDerivative(Vector derivative)
    { solver_results.Derivatives.emplace_back(std::move(derivative)); }

public:
    // Solve given equation numericaly for fixed amount of iterations
    SolverResults Solve(
        Equation current_equation,
        CalcType time_max,
        CalcType step,
        Vector initial_value,
        bool debug_print_enabled = false
    ) {
        CalcType time = 0;
        Vector value = initial_value;
        
        solver_results.clear();
        Add(value, time);
        while (time < time_max)
        {
            AddDerivative(current_equation(time, value));
            
            Vector c_0 = current_equation(time, value);
            time += step / 2;
            Vector c_1 = current_equation(time, value + (c_0 * step / 2));
            Vector c_2 = current_equation(time, value + (c_1 * step / 2));
            time += step / 2;
            Vector c_3 = current_equation(time, value + (c_2 * step));
            Vector c = c_0;
            c = c + (c_1 * 2);
            c = c + (c_2 * 2);
            c = c + c_3;
            value = value + (c * step / 6);

            for (auto& val : value)
                if (val < 0)
                    val = 0;
            Add(value, time);

            if (debug_print_enabled)
                cout << "   t = " << time << ";   N = " << value.first() << '\n';
        }

        AddDerivative(current_equation(time, value));
        return solver_results;
    }
};


/**************************************************************************************************************************
**                                                    INPUT AND OUTPUT                                                   **
***************************************************************************************************************************
**          The code below realizes means of retreiving spatial model parameters from input data ".txt" file,            **
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

    Vector getVector(int vec_size)
    {
        Vector vector_from_file;
        for (int i=0; i<vec_size; i++)
            vector_from_file.add(GetValue());
        
        return vector_from_file;
    }

    Matrix getMatrix(int mat_size)
    {
        Matrix matrix_from_file(mat_size);
        for (int i=0; i<mat_size; i++)
            for (int j=0; j<mat_size; j++)
                matrix_from_file[i][j] = GetValue();

        return matrix_from_file;
    }

    void DebugPrint(Properties props)
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
    
        auto sW = getMatrix(species);
        auto sM = getVector(species);
        
        auto b = getVector(species);
        auto d = getVector(species);
        auto dd = getMatrix(species);

        auto properties = Properties::Generate(species, sW, sM, b, d, dd);
        if (debug_print_enabled)
            DebugPrint(properties);
        properties.CheckKernels();
        return properties;
    }

    Properties GetWithoutSigmaW(CalcType sigmaW, bool debug_print_enabled = false)
    {
        int species = getInt();
    
        Matrix sW(species, sigmaW);
        auto sM = getVector(species);
        
        auto b = getVector(species);
        auto d = getVector(species);
        auto dd = getMatrix(species);

        auto properties = Properties::Generate(species, sW, sM, b, d, dd);
        if (debug_print_enabled)
            DebugPrint(properties);
        return properties;
    }

    Properties GetWithoutSigmaM(CalcType sigmaM, bool debug_print_enabled = false)
    {
        int species = getInt();
    
        auto sW = getMatrix(species);
        Vector sM(species, sigmaM);
        
        auto b = getVector(species);
        auto d = getVector(species);
        auto dd = getMatrix(species);

        auto properties = Properties::Generate(species, sW, sM, b, d, dd);
        if (debug_print_enabled)
            DebugPrint(properties);
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

    void WritePlotData(
        string experiment_number,
        string plot_title,
        string x_label,
        string y_label,
        Vector x_values,
        Vector y_values,
        int number_of_plots = 1
    ) {
        assert (number_of_plots > 0);
        assert (number_of_plots * x_values.size() == y_values.size());

        file.clear();
        file << experiment_number << '\n' << plot_title << '\n' << x_label << '\n' << y_label << '\n' << x_values.size() <<'\n';
        for (auto x : x_values.get_vector())
            file << x << '\n';
        for (auto y : y_values.get_vector())
            file << y << '\n';
    }
};


//_____________________________________________________EXPERIMENTS_________________________________________________________
class Experiments
{
    Properties Props;
    vector<Properties> PropsVector;
    string reader_file_path;
    string visdata_path;
    bool enable_debug_print;

private:
    static bool SortByFirst(const pair<CalcType, CalcType> &a, const pair<CalcType, CalcType> &b)
    { return (a.first < b.first); }

    // Comparator for sorting vector of Solver Results
    static bool SortBySigma(const pair<CalcType, SolverResults> &a, const pair<CalcType, SolverResults> &b)
    { return (a.first < b.first); }

    // Calculate first moment
    SolverResults CalculateFirstMoment(CalcType t_max)
    {
        FileReader Reader(reader_file_path);

        Props = Reader.GetProperties(enable_debug_print = true);
        Equation equation(Props);
        EulerSolver solver;

        cout << "\n==========================CALCULATIONS START==========================\n";
        clock_t Start = clock();
        
        SolverResults result = solver.Solve(equation, t_max, 0.1, equation.GenerateEquilibriumValuesForNAndC(200), enable_debug_print = true);
        
        clock_t End = clock();
        cout << "\n===========================CALCULATIONS END===========================\n";
        
        Time calculation_time(End - Start);
        cout << "TIME:  " << calculation_time << '\n';

        return result;
    }

    // Calculate first moment with varying sigma (sM or sW)
    vector<SolverResults> CalculateWithVaryingSigma(CalcType t_max, char SigmaType, Vector SigmaVec)
    {
        vector<pair<CalcType, SolverResults>> SigmasAndResults;
        
        for (int i=0; i<SigmaVec.size(); i++)
        {
            FileReader Reader(reader_file_path);
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
                    cout << "===============================ERROR!!================================\n" << endl;
                    break;
            }
        }

        int done_iterations = 0;
        omp_lock_t print_lock;
        omp_init_lock(&print_lock);

        cout << "\n==========================CALCULATIONS START==========================\n";
        clock_t Start = clock();

        #pragma omp parallel for default(none) shared(t_max, SigmaVec, SigmasAndResults, print_lock, done_iterations, cout)
        for (int i=0; i<SigmaVec.size(); i++)
        {
            auto thread_num = omp_get_thread_num();

            omp_set_lock(&print_lock);
            cout << (int) (done_iterations * 100 / SigmaVec.size()) << "%    THREAD " << thread_num << ":     ENTRY " << i << ": START\n";
            omp_unset_lock(&print_lock);
            
            auto new_props = PropsVector.at(i);
            Equation equation(new_props);
            EulerSolver solver;

            auto new_result = solver.Solve(equation, t_max, 0.1, equation.GenerateEquilibriumValuesForNAndC(200), enable_debug_print = false);
            SigmasAndResults.push_back(make_pair(SigmaVec[i], new_result));
            
            done_iterations++;
            CalcType N = new_result.Values.at(new_result.Values.size() - 1).first();
            
            omp_set_lock(&print_lock);
            cout << (int) (done_iterations * 100 / SigmaVec.size()) << "%    THREAD " << thread_num << ":     ENTRY " << i <<
                 ": FINISH    sigma = " << SigmaVec[i] << ";  N = " << N << '\n';
            omp_unset_lock(&print_lock);
        }
        clock_t End = clock();
        cout << "\n===========================CALCULATIONS END===========================\n";
        
        omp_destroy_lock(&print_lock);
        Time calculation_time(End - Start);
        cout << "TIME:  " << calculation_time << '\n';

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
        
        auto C = VectorToMatrixOfMatrix(C_vals, 1, Props.MomentsGrid.size());
        vector<pair<CalcType, CalcType>> C_r;
        bool unique;

        cout << Props.MomentsGrid.size() << "\n\n";
        for (int x=0; x<Props.MomentsGrid.size(); x++)
        {
            for (int y=0; y<Props.MomentsGrid.size(); y++)
            {
                unique = true;
                CalcType r = std::hypot(Props.MomentsGrid[x], Props.MomentsGrid[y]);
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
        sort(C_r.begin(), C_r.end(), SortByFirst);

        Vector x;
        Vector y;
        for (auto item : C_r)
        {
            CalcType func = -std::log(std::abs(item.second - N * N));
            x.add(item.first);
            y.add(func);
        }

        FileWriter VisData(visdata_path);
        VisData.WritePlotData("14.4", plot_title, "r", "-ln|C(r) - N^2|", x, y, 1);
        Visualise_Windows();
    }

    // Launch Visualiser.py     !! THIS FUNCTION WORKS ONLY ON WINDOWS !!
    void Visualise_Windows() { WinExec("python Visualiser.py", 1); }

public:
    Experiments() : visdata_path("Data\\VisData.txt") {}

    ~Experiments() { PropsVector.clear(); }

    void Experiment_14_4(int input_data)
    {
        string plot_title;
        switch(input_data)
        {
            case 0:
                reader_file_path = "Data\\14\\4\\a.txt";
                plot_title = "(a)";
                break;
            case 1:
                reader_file_path = "Data\\14\\4\\b.txt";
                plot_title = "(b)";
                break;
            case 2:
                reader_file_path = "Data\\14\\4\\c.txt";
                plot_title = "(c)";
                break;
        }
        PropsVector.clear();
        auto solver_results = CalculateFirstMoment(100);
    
        Vector N;
        Vector CSquareRoot;
        for (auto Values : solver_results.Values)
        {
            CalcType N_value = Values.first();
            N.add(N_value);
            
            CalcType C_at_cutoff_distance = Values[Props.MomentsGrid.size()];
            CSquareRoot.add(std::sqrt(C_at_cutoff_distance));
        }
        Vector values = N;
        values.add(CSquareRoot);

        FileWriter VisData(visdata_path);
        VisData.WritePlotData("14.4", plot_title, "Time", "N", solver_results.Times, values, 2);
        Visualise_Windows();
    }

    void Experiment_14_6(char SigmaType, string RadiusType)
    {
        reader_file_path = "Data\\14\\6\\a.txt";
        string plot_title = "(a)";
        if(SigmaType == 'M')
        {
            reader_file_path = "Data\\14\\6\\b.txt";
            plot_title = "(b)";
        }
        PropsVector.clear();
        auto SigmaVec = Vector::setLinSpaced(160, 1e-08, 0.2);
        auto solver_results = CalculateWithVaryingSigma(200, SigmaType, SigmaVec);
        
        Vector ConvergedN;
        Vector ConvergedCSquareRoot;
        assert (PropsVector.size() == solver_results.size());
        assert (PropsVector.size() == SigmaVec.size());
        for (int i=0; i<PropsVector.size(); i++)
        {
            auto values = solver_results.at(i).Values;
            auto converged_N_and_C = values.at(values.size() - 1);
            CalcType converged_N_value = converged_N_and_C.first();
            ConvergedN.add(converged_N_value);

            CalcType converged_C_at_cutoff_distance = converged_N_and_C[PropsVector[i].MomentsGrid.size()];
            ConvergedCSquareRoot.add(std::sqrt(converged_C_at_cutoff_distance));
        }
        Vector values = ConvergedN;
        values.add(ConvergedCSquareRoot);
        
        FileWriter VisData(visdata_path);
        VisData.WritePlotData("14.6", plot_title, RadiusType, "Equilibrium density", SigmaVec, values, 2);
        Visualise_Windows();
    }

    void Custom_Experiment()
    {
        PropsVector.clear();
        reader_file_path = "Data\\Custom\\experiment.txt";
        auto solver_results = CalculateFirstMoment(100);
    
        Vector N;
        Vector CSquareRoot;
        for (auto Values : solver_results.Values)
        {
            CalcType N_value = Values.first();
            N.add(N_value);

            CalcType C_at_cutoff_distance = Values[Props.MomentsGrid.size()];
            CSquareRoot.add(std::sqrt(C_at_cutoff_distance));
        }
        Vector values = N;
        values.add(CSquareRoot);

        FileWriter VisData(visdata_path);
        VisData.WritePlotData("14.4", " ", "Time", "N", solver_results.Times, values, 2);
        Visualise_Windows();
    }

    void KernelTest()
    {
        FileReader Reader("Data\\14\\4\\a.txt");
        auto Props = Reader.GetProperties(enable_debug_print = true);
        
        auto kernel = Props.W(0, 0);
        //auto kernel = Props.M(0);
    
        vector<pair<CalcType, CalcType>> result;
        bool unique;

        for (int x=0; x<Props.KernelsGrid.size(); x++)
        {
            for (int y=0; y<Props.KernelsGrid.size(); y++)
            {
                unique = true;
                CalcType r = std::hypot(Props.KernelsGrid[x], Props.KernelsGrid[y]);
                for (auto item : result)
                {
                    cout << "x = " << Props.KernelsGrid[x] << "    y = " << Props.KernelsGrid[y] << "     r = " << r << "    kernel = " << kernel[x][y] << '\n';
                    if (r == item.first)
                        unique = false;
                }
                if (unique)
                    result.push_back(make_pair(r, kernel[x][y]));
            }
        }
        sort(result.begin(), result.end(), SortByFirst);

        Vector x;
        Vector y;
        for (auto item : result)
        {
            x.add(item.first);
            y.add(item.second);
        }
        cout << "MAX NON-ZERO NUMBERS: " << Props.MaxNumForW << '\n';
        //cout << "MAX NON-ZERO NUMBERS: " << Props.MaxNumForM << '\n';
        
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
