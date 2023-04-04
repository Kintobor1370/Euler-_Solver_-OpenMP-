#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <utility> 
#include <fstream>
#include <cassert>
#include <algorithm>
#include <thread>
#include <mutex>
#include <time.h>
#include <windows.h>
#include <conio.h>
#include <omp.h>
#include <tuple>

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


//________________________________________________CALCULATION PARAMETERS___________________________________________________
struct Parameters
{
    Vector MidPointsKernels;                                            // The mid-points of bins that correspond to kernels indices
    CalcType MidPointsKernelsWidth;
    Vector MidPointsMoments;                                            // The mid-points of bins that correspond to second moments indices

public:
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
        int max_bin = 0; 

        // Compute
        for (int x=0; x<m.size(); x++)
        {
            for (int y=0; y<m.size(); y++)
            {
                const auto r = hypot(m[x], m[y]);
                if(r < R)
                {
                    switch(type)
                    {
                        case Normal:
                            auto norm = 2 * PI * s * s * (1 - exp((-R * R) / (2 * s * s)));
                            auto weight = exp((-r * r) / (2 * s * s));
                            result[x][y] = weight / norm;
                            break;
                    }
                    max_bin = std::max({max_bin, x, y});
                }
                else
                    result[x][y] = 0;
            }
        }

        // Renormalize
        CalcType integral = 0;
        for (int x=1-m.size(); x<m.size(); x++)
            for (int y=1-m.size(); y<m.size(); y++)
                integral += result[abs(x)][abs(y)];
            
        integral *= pow(width, 2);
        result = result / integral;

        return make_pair(result, max_bin);
    }

    //..............................GENERATE CALCULATION PARAMETERS
    static Parameters Generate(
        int moments_num,                                                    // Number of mid-points for moments
        CalcType moments_width,                                             // Bin width for moments
        int kernels_num,                                                    // Number of mid-points for kernels
        CalcType kernels_width                                              // Bin width for kernels
    ) {
        Vector kernels_partition = Vector::setLinSpaced(kernels_num, 0, kernels_num - 1) * kernels_width;   
        Vector moments_partition = Vector::abs(-Vector::log(-Vector::setLinSpaced(moments_num, 0, moments_num - 1) / moments_num + 1) * moments_width * moments_num);

        Parameters params;
        params.MidPointsKernels = kernels_partition;
        params.MidPointsKernelsWidth = kernels_width;
        params.MidPointsMoments = moments_partition;
        return params;
    }
};


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

public:
    //..............................GENERATE MODEL PROPERTIES
    static Properties Generate(
        int num,                                                            // Number of species
        Matrix sW,                                                          // Sigma W
        Vector sM,                                                          // Sigma M
        Vector b,                                                           
        Vector d,
        Matrix dd,
        Parameters params,                                                  // Calculation parameters
        KernelType w_type = Normal,                                         // W kernel type
        KernelType m_type = Normal                                          // M kernel type
    ) {            
        MatrixOfMatrix w(num, num);
        Matrix w_max_bin(num, num);
        Matrix max_radius_W;
        if (w_type == Normal)
            max_radius_W = sW * 3;
        
        for (int i=0; i<num; i++)
        {
            for (int j=0; j<num; j++)
            {
                auto kernel_W = params.GenerateKernel(sW[i][j], max_radius_W[i][j], w_type);
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
            auto kernel_M = params.GenerateKernel(sM[i], max_radius_M[i], m_type);
            m(i) = kernel_M.first;
            m_max_bin[i] = kernel_M.second;
        }

        Properties props;
        props.SpeciesNumber = num;
        props.b = b;
        props.d = d;
        props.d_prime = dd;
        props.M = m;
        props.W = w;
        props.WMaxBin = w_max_bin;
        props.MMaxBin = m_max_bin;
        props.SigmaM = sM;
        props.SigmaW = sW;
        return props;
    }

    //..............................GENERATE MODEL PROPERTIES FOR ONE SPECIES
    static Properties GenerateOneSpecies(
        CalcType sW,
        CalcType sM,
        CalcType b,
        CalcType d,
        CalcType dd,
        Parameters params,
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
            params,
            w_type,
            m_type
        );
    }
};


//__________________________________EQUATIONS FOR THE FIRST AND SECOND MOMENT DYNAMICS_____________________________________
class Equation
{
    Properties props;                                                   // Model properties
    Parameters params;                                                  // Calculation parameters
    Vector N;                                                           // Vector of current first moments
    MatrixOfMatrix C;                                                   // Current second moments

private:
    // Make C (Matrix of Matrix) from Vector of values
    MatrixOfMatrix MakeC(Vector C_vals)
    {
        assert(C_vals.size() == props.SpeciesNumber * props.SpeciesNumber * params.MidPointsMoments.size() * params.MidPointsMoments.size());

        MatrixOfMatrix result = VectorToMatrixOfMatrix(C_vals, props.SpeciesNumber, params.MidPointsMoments.size()); 
        return result;
    }

    // Make Vector of N (Vector) and C (Matrix of Matrix) values
    Vector MakeVector(Vector N, MatrixOfMatrix C)
    {
        Vector result = N;
        Vector C_To_Vector = MatrixOfMatrixToVector(C);
        assert(result.size() == props.SpeciesNumber);
        assert(C_To_Vector.size() == props.SpeciesNumber * props.SpeciesNumber * params.MidPointsMoments.size() * params.MidPointsMoments.size());
        
        result.add(C_To_Vector);
        return result;
    }

    // Calculte integral W(xi)C(xi)dxi
    Matrix CalculateIntegralOfWAndC(MatrixOfMatrix C, Vector N)
    {
        Matrix result(props.SpeciesNumber, props.SpeciesNumber);

        for (int i=0; i<props.SpeciesNumber; i++) 
            for (int j=0; j<props.SpeciesNumber; j++)
                result[i][j] = IntegrateWithPairDensity(0, 0, props.WMaxBin[i][j], props.W(i, j), C(i, j), N[i] * N[j]);

        return result;
    }

    MatrixOfMatrix CalculateIntegralOfMAndC(MatrixOfMatrix C, Vector N)
    {
        MatrixOfMatrix result(props.SpeciesNumber, props.SpeciesNumber);
        
        for (int i=0; i<props.SpeciesNumber; i++)
            for (int j=0; j<props.SpeciesNumber; j++)
            {
                Matrix result_elem(params.MidPointsMoments.size(), params.MidPointsMoments.size());

                for (int xBin=0; xBin<params.MidPointsMoments.size(); xBin++)
                    for (int yBin=0; yBin<params.MidPointsMoments.size(); yBin++)
                        result_elem[xBin][yBin] = IntegrateWithPairDensity(xBin, yBin, props.MMaxBin[i], props.M(i), C(i, j), N[i] * N[j]);

                result(i, j) = result_elem;
            }
        return result;
    }

    CalcType IntegrateWithPairDensity(int x_shift, int y_shift, int binCounts, Matrix kernels, Matrix C, CalcType asymptotic_value)
    {
        CalcType integral = 0;
        
        for (int i=-binCounts; i<=binCounts; i++)
        {
            for (int j=-binCounts; j<=binCounts; j++)
            {
                const CalcType sum_x = abs(params.MidPointsMoments[x_shift] + (i > 0 ? params.MidPointsKernels[i] : -params.MidPointsKernels[-i]));
                const CalcType sum_y = abs(params.MidPointsMoments[y_shift] + (j > 0 ? params.MidPointsKernels[j] : -params.MidPointsKernels[-j]));

                CalcType density = InterpolateDensity(C, params.MidPointsMoments, sum_x, sum_y, asymptotic_value);
                integral += kernels[abs(i)][abs(j)] * density;
            }
        }
        integral *= pow(params.MidPointsKernelsWidth, 2);
        return integral;
    }

    CalcType InterpolateDensity(Matrix density, Vector midPoints, CalcType x, CalcType y, CalcType asymptotic_value)
    {
        assert(density.Rows() == midPoints.size());
        assert(density.Cols() == midPoints.size());

        if (x >= midPoints[midPoints.size() - 1] || y >= midPoints[midPoints.size() - 1]) 
            return asymptotic_value;
        assert(std::is_sorted(midPoints.begin(), midPoints.end()));

        // Find cell corned by mid-points
        const auto xMaxIt = std::upper_bound(midPoints.begin(), midPoints.end(), x);
        const auto xMinIt = xMaxIt - 1;
        const auto yMaxIt = std::upper_bound(midPoints.begin(), midPoints.end(), y);
        const auto yMinIt = yMaxIt - 1;

        assert(xMaxIt != midPoints.end());
        assert(yMaxIt != midPoints.end());

        // Get index found cells
        const auto xMaxBin = std::distance(midPoints.begin(), xMaxIt);
        const auto xMinBin = std::distance(midPoints.begin(), xMinIt);
        const auto yMaxBin = std::distance(midPoints.begin(), yMaxIt);
        const auto yMinBin = std::distance(midPoints.begin(), yMinIt);

        // Compute weights
        CalcType wx = (x - *xMinIt) / (*xMaxIt - *xMinIt);
        CalcType wy = (y - *yMinIt) / (*yMaxIt - *yMinIt);

        // Compute weighted value
        return  density[xMinBin][yMinBin] * (1. - wx) * (1. - wy) +
                density[xMaxBin][yMinBin] * wx * (1. - wy) +
                density[xMinBin][yMaxBin] * (1. - wx) * wy +
                density[xMaxBin][yMaxBin] * wx * wy;
    }

public:
    Equation(Properties a, Parameters b) : props(a), params(b) {}

    // Make Vector of equilibrium values for N (Vector) and C (Matrix of Matrix). Here C=N*N
    Vector GenerateEquilibriumValuesForNAndC(CalcType N_vals)
    {
        int N_size = props.SpeciesNumber;
        int C_size = props.SpeciesNumber * props.SpeciesNumber * params.MidPointsMoments.size() * params.MidPointsMoments.size();
        Vector result(N_size + C_size);
        
        for (int i=0; i<N_size; i++)
            result[i] = N_vals;
        for (int i=N_size; i<result.size(); i++)
            result[i] = N_vals * N_vals;

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
        Vector dPrimeSum;                                               // Death contribution, density dependent
        for(int j=0; j<props.SpeciesNumber; j++)
            dPrimeSum = dPrimeSum + props.d_prime.Col(j) * WC.Col(j);
        
        dN = dN - dPrimeSum;

        //.......................................CALCULATE SECOND MOMENT
        for (int i=0; i<props.SpeciesNumber; i++)
            for (int j=0; j<props.SpeciesNumber; j++)
            {
                Matrix dC_elem(params.MidPointsMoments.size(), params.MidPointsMoments.size());

                for (int xBin=0; xBin<params.MidPointsMoments.size(); xBin++)
                    for (int yBin=0; yBin<params.MidPointsMoments.size(); yBin++)
                    {
                        CalcType result;

                        // Birth contribution, density independent
                        result = props.b[i] * MC(i, j)[xBin][yBin] + props.b[j] * MC(j, i)[xBin][yBin]; 

                        // Birth contribution, Kronecker symbols
                        if (i == j)
                        {
                            result += 2 * props.b[i] * N[i] * InterpolateDensity(
                                props.M(i),
                                params.MidPointsKernels,
                                params.MidPointsMoments[xBin],
                                params.MidPointsMoments[yBin],
                                0
                            );
                        }

                        // Death contribution, density independent
                        result -= props.d[i] * C(i, j)[xBin][yBin];
                        result -= props.d[j] * C(j, i)[xBin][yBin];

                        // Simple closure
                        for (int k=0; k<props.SpeciesNumber; k++)
                        {
                            result -= props.d_prime[i][k] * C(i, j)[xBin][yBin] / N[i] * WC[i][k];
                            result -= props.d_prime[j][k] * C(j, i)[xBin][yBin] / N[j] * WC[j][k];
                        }
                        
                        result -= props.d_prime[i][j] * C(i, j)[xBin][yBin] * InterpolateDensity(
                            props.W(i, j),
                            params.MidPointsKernels,
                            params.MidPointsMoments[xBin],
                            params.MidPointsMoments[yBin],
                            0
                        );
                        result -= props.d_prime[j][i] * C(j, i)[xBin][yBin] * InterpolateDensity(
                            props.W(j, i),
                            params.MidPointsKernels,
                            params.MidPointsMoments[xBin],
                            params.MidPointsMoments[yBin],
                            0
                        );
                        
                        dC_elem[xBin][yBin] = result;
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
        CalcType t = 0;
        Vector y = y_0;
        
        solver_results.clear();
        Add(y, t);
        while (t <= time_max)
        {
            Vector y_new = func(t, y);
            
            AddDerivative(y_new);
            y = y + y_new * step;
            for (auto& v : y)
                if (v < 0)
                    v = 0;
            
            t += step;
            Add(y, t);

            if (debug_print_enabled)
                cout << "   t = " << t << ";   N = " << y[0] << '\n';
        }

        AddDerivative(func(t, y));
        return solver_results;
    }
};


//________________________________________________RUNGE-KUTTA 4 METHOD SOLVER_______________________________________________
class R_K_4_Solver
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
        CalcType t = 0;
        Vector y = y_0;
        
        solver_results.clear();
        Add(y, t);
        while (t <= time_max)
        {
            AddDerivative(func(t, y));
            
            Vector c_0 = func(t, y);
            t += step / 2;
            Vector c_1 = func(t, y + (c_0 * step / 2));
            Vector c_2 = func(t, y + (c_1 * step / 2));
            t += step / 2;
            Vector c_3 = func(t, y + (c_2 * step));
            Vector c = c_0;
            c = c + (c_1 * 2);
            c = c + (c_2 * 2);
            c = c + c_3;
            y = y + (c * step / 6);

            for (auto& v : y)
                if (v < 0)
                    v = 0;
            Add(y, t);

            if (debug_print_enabled)
                cout << "   t = " << t << ";   N = " << y[0] << '\n';
        }

        AddDerivative(func(t, y));
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


public:
    FileReader(string str) : filename(str) { file.open(filename); }

    ~FileReader() { file.close(); }

    Parameters GetParameters(bool debug_print_enabled = false)
    {
        int m_num = getInt();
        CalcType m_width = getCalcType();
        int k_num = getInt();
        CalcType k_width = getCalcType();

        if (debug_print_enabled)
            cout << "\nCALCULATION PARAMETERS:"
                 << "\n   Number of points for moments: " << m_num
                 << "\n   Moments bin width: " << m_width
                 << "\n   Number of points for kernels: " << k_num
                 << "\n   Kernels bin width: " <<k_width
                 << "\n\n";

        return Parameters::Generate(m_num, m_width, k_num, k_width);
    }

    Properties GetProperties(Parameters params, bool debug_print_enabled = false)
    {
        int species = getInt();
    
        Matrix sW = getMatrix(species, species);
        Vector sM = getVector(species);
        
        Vector b = getVector(species);
        Vector d = getVector(species);
        Matrix dd = getMatrix(species, species);

        if (debug_print_enabled)
            cout << "MODEL PROPERTIES:\n"
                 << "   Species number: " << species << '\n'
                << "   Sigma W:  " << sW
                << "   Sigma M:  " << sM
                << "         b:  " << b
                << "         d:  " << d
                << "         d': " << dd
                << "\n\n";
        return Properties::Generate(species, sW, sM, b, d, dd, params);
    }

    Properties GetWithoutSigmaW(CalcType sigmaW, Parameters params, bool debug_print_enabled = false)
    {
        int species = getInt();
    
        Matrix sW = Matrix(species, species, sigmaW);
        Vector sM = getVector(species);
        
        Vector b = getVector(species);
        Vector d = getVector(species);
        Matrix dd = getMatrix(species, species);

        if (debug_print_enabled)
            cout << "MODEL PROPERTIES:\n"
                 << "   Species number: " << species << '\n'
                << "   Sigma W:  " << sW
                << "   Sigma M:  " << sM
                << "         b:  " << b
                << "         d:  " << d
                << "         d': " << dd
                << "\n\n";

        return Properties::Generate(species, sW, sM, b, d, dd, params);
    }

    Properties GetWithoutSigmaM(CalcType sigmaM, Parameters params, bool debug_print_enabled = false)
    {
        int species = getInt();
    
        Matrix sW = getMatrix(species, species);
        Vector sM = Vector(species, sigmaM);
        
        Vector b = getVector(species);
        Vector d = getVector(species);
        Matrix dd = getMatrix(species, species);

        if (debug_print_enabled)
            cout << "MODEL PROPERTIES:\n"
                 << "   Species number: " << species << '\n'
                << "   Sigma W:  " << sW
                << "   Sigma M:  " << sM
                << "         b:  " << b
                << "         d:  " << d
                << "         d': " << dd
                << "\n\n";
        
        return Properties::Generate(species, sW, sM, b, d, dd, params);
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

    void WriteVData(string ExperimentNum, string xLabel, string yLabel, vector<CalcType> xValues, vector<CalcType> yValues)
    {
        assert(xValues.size() == yValues.size());

        file.clear();
        file << ExperimentNum << '\n'  << xLabel << '\n' << yLabel << '\n' << xValues.size() <<'\n';
        for (auto x : xValues)
            file << x << '\n';
        for (auto y : yValues)
            file << y << '\n';
    }
};


//_____________________________________________________EXPERIMENTS_________________________________________________________
class Experiments
{
    Parameters Params;
    vector<Properties> PropsVector;
    string reader_path;
    string writer_path;
    string vdata_path;
    bool enable_debug_print;

private:
    // Read calculation parameters from "Parameters.txt"
    void ReadParams()
    {
        FileReader params_reader("Data\\Parameters.txt");
        Params = params_reader.GetParameters(enable_debug_print = true);
    }

    // Comparator for sorting vector of Solver Results
    static bool SortBySigma(const pair<CalcType, SolverResults> &a, const pair<CalcType, SolverResults> &b)
    { return (a.first < b.first); }

    // Display recorded time of calculations
    static tuple<string, string, string> MakeTime(clock_t start, clock_t end)
    {
        int seconds = (int) (end - start) / CLOCKS_PER_SEC;
        int hours = seconds / 3600;
        int minutes = (seconds % 3600) / 60;
        seconds = seconds % 60;

        string str_minutes = minutes / 10 == 0 ? "0" + to_string(minutes) : to_string(minutes);
        string str_seconds = seconds / 10 == 0 ? "0" + to_string(seconds) : to_string(seconds);

        return make_tuple(to_string(hours), str_minutes, str_seconds);
    }

    // Calculate first moment
    SolverResults CalculateFirstMoment(CalcType t_max)
    {
        FileReader Reader(reader_path);
        ReadParams();

        auto Props = Reader.GetProperties(Params, enable_debug_print = true);
        Equation equation(Props, Params);   
        EulerSolver solver;

        cout << "\n==========================CALCULATIONS START==========================\n";
        clock_t Start = clock();
        
        SolverResults result = solver.Solve(equation, t_max, 0.1, equation.GenerateEquilibriumValuesForNAndC(200), enable_debug_print = true);
        
        clock_t End = clock();
        cout << "\n===========================CALCULATIONS END===========================\n";
        
        auto Time = MakeTime(Start, End);
        cout << "TIME:  " << get<0>(Time) << ':' << get<1>(Time) << ':' << get<2>(Time) << '\n';

        return result;
    }

    // Calculate first moment with varying sigma (sM or sW)
    vector<SolverResults> CalculateWithVaryingSigma(CalcType t_max, char SigmaType, Vector SigmaVec)
    {
        vector<pair<CalcType, SolverResults>> SigmasAndResults;
        ReadParams();
        for (int i=0; i<SigmaVec.size(); i++)
        {
            FileReader Reader(reader_path);
            cout << "ENTRY " << i << ":\n";
            switch(SigmaType)
            {
                case 'W':
                    PropsVector.push_back(Reader.GetWithoutSigmaW(SigmaVec[i], Params, enable_debug_print = true));
                    break;
                case 'M':
                    PropsVector.push_back(Reader.GetWithoutSigmaM(SigmaVec[i], Params, enable_debug_print = true));
                    break;
                default:
                    cerr << "===============================ERROR!!!===============================\n" << endl;
                    break;
            }
        }

        cout << "\n==========================CALCULATIONS START==========================\n";
        clock_t Start = clock();

        #pragma omp parallel for default(none) shared(t_max, SigmaVec, SigmasAndResults, cout)
        for (int i=0; i<SigmaVec.size(); i++)
        {
            cout << "THREAD " << omp_get_thread_num() << ":     ENTRY " << i << ": START\n";
            
            auto new_props = PropsVector.at(i);
            Equation equation(new_props, Params);
            EulerSolver solver;

            auto new_result = solver.Solve(equation, t_max, 0.1, equation.GenerateEquilibriumValuesForNAndC(200), enable_debug_print = true);
            SigmasAndResults.push_back(make_pair(SigmaVec[i], new_result));
            
            cout << "THREAD " << omp_get_thread_num() << ":     ENTRY " << i << ": FINISH\n";
        }
        clock_t End = clock();
        cout << "\n===========================CALCULATIONS END===========================\n";
        
        auto Time = MakeTime(Start, End);
        cout << "TIME:  " << get<0>(Time) << ':' << get<1>(Time) << ':' << get<2>(Time) << '\n';

        sort(SigmasAndResults.begin(), SigmasAndResults.end(), SortBySigma);
        vector<SolverResults> results;
        for (auto result : SigmasAndResults)
            results.push_back(result.second);
        return results;
    }

    // Launch Visualiser.py 
    // !! WORKS ONLY ON WINDOWS !!
    void Visualise_Windows() { WinExec("python Visualiser.py", 1); }

public:
    Experiments() : writer_path("Data\\Results.txt"), vdata_path("Data\\VData.txt") {}

    ~Experiments() { PropsVector.clear(); }

    void Experiment_14_4(int input_data)
    {
        switch(input_data)
        {
            case 0:
                reader_path = "Data\\14.4\\a.txt";
                break;
            case 1:
                reader_path = "Data\\14.4\\b.txt";
                break;
            case 2:
                reader_path = "Data\\14.4\\c.txt";
                break;
        }
        PropsVector.clear();
        SolverResults Results = CalculateFirstMoment(100);
    
        vector<CalcType> N;
        for (auto item : Results.Values)
            N.push_back(item[0]);
        
        FileWriter Writer(writer_path);
        Writer.WriteResult("t", "N", Results.Times, N);

        FileWriter VData(vdata_path);
        VData.WriteVData("14.4", "Time", "N", Results.Times, N);
        Visualise_Windows();
    }

    void Experiment_14_6(char SigmaType, string RadiusType)
    {
        reader_path = "Data\\14.6\\a.txt";
        if(SigmaType == 'M')
            reader_path = "Data\\14.6\\b.txt";

        PropsVector.clear();
        Vector SigmaVec = Vector::setLinSpaced(160, 0.00000001, 0.2);
        
        vector<SolverResults> solverResults = CalculateWithVaryingSigma(300, SigmaType, SigmaVec);
        vector<CalcType> ConvergedNVec;
        
        for (auto Result : solverResults)
        {
            auto FinalValue = Result.Values.at(Result.Values.size()-1);
            CalcType ConvergedN = FinalValue[0]; 
            ConvergedNVec.push_back(ConvergedN);
        }

        FileWriter Writer(writer_path);
        Writer.WriteResult("s_" + SigmaType, "N", SigmaVec.get_vector(), ConvergedNVec);
        
        FileWriter VData(vdata_path);
        VData.WriteVData("14.6", RadiusType, "Equilibrium density", SigmaVec.get_vector(), ConvergedNVec);
        Visualise_Windows();
    }
};


/**************************************************************************************************************************
**                                                        MAIN BODY                                                      **
**************************************************************************************************************************/
int main()
{
    Experiments Run;
    string commands[] =
    {
        " >  14.4.a\n    14.4.b\n    14.4.c\n    14.6.a\n    14.6.b\n",
        "    14.4.a\n >  14.4.b\n    14.4.c\n    14.6.a\n    14.6.b\n",
        "    14.4.a\n    14.4.b\n >  14.4.c\n    14.6.a\n    14.6.b\n",
        "    14.4.a\n    14.4.b\n    14.4.c\n >  14.6.a\n    14.6.b\n",
        "    14.4.a\n    14.4.b\n    14.4.c\n    14.6.a\n >  14.6.b\n"
    };

    system("cls");
    int current_index = 0;
    cout << "Choose experiment:\n" + commands[current_index] + "(Navigate with up arrow & down arrow\n Press ENTER to confirm\n Press ESC to quit)";
    char input = getch();

    while(input != KEY_ESC)
    {
        switch (input)
        {
            case KEY_UP:
                current_index = (current_index + 4) % 5;
                system("cls");
                cout << "Choose experiment:\n" + commands[current_index] + "(Navigate with up arrow & down arrow.\n Press ENTER to confirm\n Press ESC to quit)";
                input = getch();
                break;

            case KEY_DOWN:
                current_index = (current_index + 1) % 5;
                system("cls");
                cout << "Choose experiment:\n" + commands[current_index] + "(Navigate with up arrow & down arrow.\n Press ENTER to confirm\n Press ESC to quit)";
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
                }

            default:
                input = getch();
                break;
        }
    }
    system("cls");

    return 0;
}