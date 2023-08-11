/**********************************************************************************************
**                                      TYPES HEADER FILE                                    **
**  This file contains classes for mathematical data structures used in Solver.cpp program.  **
**********************************************************************************************/

#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;

template<typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

//___________________________________MAIN CALCULATION TYPE_____________________________________
using CalcType = double;


//_________________________________VECTOR OF CALCULATION TYPE__________________________________
class Vector
{
    vector<CalcType> vec;

public:
    Vector(vector<CalcType> base_vec) : vec(base_vec) {}

    Vector(const int n=0) { vec.resize(n); }

    Vector(const int n, CalcType val)
    {
        for (int i=0; i<n; i++)
            vec.push_back(val);
    }

    //....................Returns the original vector (required for drawing plots)
    vector<CalcType> get_vector() { return vec; }

    //....................Vector indexing operator
    CalcType& operator [] (const int index) { return vec.at(index); }

    //....................Returns true if vector is empty
    bool empty() { return vec.empty(); }

    //....................Returns Vector's size
    int size() { return vec.size(); }

    //....................Modifies Vector's size
    void resize(int n) { vec.resize(n); }

    //....................Clears Vector
    void clear() { vec.clear(); }
    
    //....................Sorts Vector
    void sort() { std::sort(vec.begin(), vec.end()); }

    //....................Adds an element to the back of Vector
    void add(CalcType new_elem) { vec.push_back(new_elem); }

    //....................Adds a Vector to the back of Vector
    void add(Vector new_vec)
    {
        for (int i=0; i<new_vec.size(); i++)
            vec.push_back(new_vec[i]);
    }

    //....................Adds a unique element to the back of Vector. Returns true if the element was added
    bool add_unique(CalcType new_elem)
    {
        if (std::find(vec.begin(), vec.end(), new_elem) == vec.end())
        {
            vec.push_back(new_elem);
            return true;
        }
        return false;
    }
    
    //....................Returns first element of Vector
    CalcType first() { return vec.at(0); }

    //....................Returns last element of Vector
    CalcType last() { return vec.at(vec.size() - 1); }

    //....................Returns first n elements of Vector
    Vector head(int n)                  
    {
        Vector res;
        for (int i=0; i<n; i++)
            res.add(vec.at(i));
        return res; 
    }

    //....................Returns last n elements of Vector
    Vector tail(int n)
    {
        Vector res;
        for (int i=n; i>0; i--)
            res.add(vec.at(vec.size() - i));
        return res; 
    }

    //....................Erases first n elements of Vector
    void erase_first(int n)
    {
        for (int i=0; i<n; i++)
            vec.erase(vec.begin());
    }

    //....................Returns vector element of the lowest value
    CalcType min_value() { return *min_element(vec.begin(), vec.end()); }

    //....................Returns vector element of the highest value
    CalcType max_value() { return *max_element(vec.begin(), vec.end()); }

    //....................Returns a sum of vector elements
    CalcType sum() { return std::accumulate(vec.begin(), vec.end(), 0.00000000); }
    
    //....................Returns an average value of vector elements
    CalcType avg_value() { return vec.empty() ? 0 : std::accumulate(vec.begin(), vec.end(), 0.00000000) / vec.size(); }

    //....................Sets a linearly spaced Vector
    static Vector setLinSpaced(CalcType num, CalcType min_val, CalcType max_val)
    {
        Vector result;
        CalcType step = (max_val - min_val) / (num - 1);

        for (int i=0; i<num-1; i++)
            result.add(min_val + i * step);
        result.add(max_val);

        return result;
    }

    //....................Natural logarithm of vector elements
    static Vector log(Vector arg)
    {
        Vector result;
        for (int i=0; i<arg.size(); i++)
            result.add(std::log(arg[i]));
        
        return result;
    }

    //....................Absolute values of vector elements
    static Vector abs(Vector arg)
    {
        Vector result;
        for (int i=0; i<arg.size(); i++)
            result.add(std::abs(arg[i]));

        return result;
    }

    //....................Vector + calculation type value
    Vector operator + (CalcType add_val)
    {
        Vector sum;

        for (int i=0; i<vec.size(); i++)
            sum.add(vec.at(i) + add_val);

        return sum;
    }

    //....................Vector - calculation type value
    Vector operator - (CalcType sub_val)
    {
        Vector sub_result;

        for (int i=0; i<vec.size(); i++)
            sub_result.add(vec.at(i) - sub_val);

        return sub_result;
    }
    
    //....................-Vector
    Vector operator - ()
    {
        Vector minus_vec;

        for (int i=0; i<vec.size(); i++)
            minus_vec.add(-1 * vec.at(i));
        
        return minus_vec;
    }

    //....................Vector * calculation type value
    Vector operator * (CalcType com_val)
    {
        Vector composition;

        for (int i=0; i<vec.size(); i++)
            composition.add(vec.at(i) * com_val);

        return composition;
    }

    //....................Vector / calculation type value
    Vector operator / (CalcType div_val)
    {
        Vector quotient;

        for (int i=0; i<vec.size(); i++)
            quotient.add(vec.at(i) / div_val);

        return quotient;
    }

    //....................Vector + Vector
    Vector operator + (Vector add_vec)
    {
        int max_size = vec.size();
        max_size = max(max_size, add_vec.size());
        Vector sum;

        vec.resize(max_size);
        add_vec.resize(max_size);
        for (int i=0; i<max_size; i++)
            sum.add(vec.at(i) + add_vec[i]);

        return sum;
    }

    //....................Vector - Vector
    Vector operator - (Vector sub_vec)
    {
        int max_size = vec.size();
        max_size = max(max_size, sub_vec.size());
        Vector sub_result;

        vec.resize(max_size);
        sub_vec.resize(max_size);
        for (int i=0; i<max_size; i++)
            sub_result.add(vec.at(i) - sub_vec[i]);

        return sub_result;
    }

    //....................Vector * Vector
    Vector operator * (Vector com_vec)
    {
        int max_size = vec.size();
        max_size = max(max_size, com_vec.size());
        Vector composition;

        vec.resize(max_size);
        com_vec.resize(max_size);
        for (int i=0; i<max_size; i++)
            composition.add(vec.at(i) * com_vec[i]);

        return composition;
    }

    //....................Vector / Vector
    Vector operator / (Vector div_vec)
    {
        int max_size = vec.size();
        max_size = max(max_size, div_vec.size());
        Vector quotient;

        vec.resize(max_size);
        div_vec.resize(max_size);
        for (int i=0; i<max_size; i++)
            quotient.add(vec.at(i) / div_vec[i]);

        return quotient;
    }

    vector<CalcType>::iterator begin() { return vec.begin(); }

    vector<CalcType>::iterator end() { return vec.end(); }

    friend ostream& operator << (ostream& out, const Vector &v)
	{
		out << v.vec.at(0);
		for(int i=1; i<v.vec.size(); i++)
			out << ", " << v.vec.at(i);
		out << '\n';
		return out;
	}

    pair<int, CalcType> MaxStep()
    {
        CalcType max_step = vec.at(1) - vec.at(0);
        int index = 1;
        for (int i=2; i<vec.size(); i++)
        {
            CalcType step = vec.at(i) - vec.at(i-1);
            if (step > max_step)
            {
                max_step = step;
                index = i;
            }
        }
        return make_pair(index, max_step);
    }
};


//______________________________SQUARE MATRIX OF CALCULATION TYPE_____________________________
class Matrix
{
    vector<Vector> mat;

public:
    Matrix(const int size=0) { resize(size); }

    Matrix(const int size, CalcType val)
    {
        resize(size);
        for (int i=0; i<size; i++)
            for (int j=0; j<size; j++)
                mat.at(i)[j] = val;
    }

    Vector& operator [] (const int row_num) { return mat.at(row_num); }

    int size() { return mat.size(); }

    void resize(const int size)
    {
        mat.resize(size);
        for (int i=0; i<size; i++)
            mat[i].resize(size);
    }

    Vector Row(const int row_num) { return mat.at(row_num); }

    Vector Col(const int col_num)
    {
        Vector col;
        for (int i=0; i<mat.size(); i++)
            col.add(mat.at(i)[col_num]);
        return col;
    }

    void resize(const int r, const int c)
    {
        mat.resize(r);
        for (int i=0; i<r; i++)
            mat[i].resize(c);
    }

    //....................Returns an element of Matrix which has the lowest value
    CalcType min_value()
    {
        vector<CalcType> min_values_vec;
        for (int i=0; i<mat.size(); i++)
            min_values_vec.push_back(mat[i].min_value());
        
        return *max_element(min_values_vec.begin(), min_values_vec.end());
    }

    //....................Returns an element of Matrix which has the highest value
    CalcType max_value()
    {
        vector<CalcType> max_values_vec;
        for (int i=0; i<mat.size(); i++)
            max_values_vec.push_back(mat[i].max_value());
        
        return *max_element(max_values_vec.begin(), max_values_vec.end());
    }

    //....................Matrix + calculation type value
    Matrix operator + (CalcType add_val)
    {
        Matrix sum(mat.size(), mat.at(0).size());

        for (int i=0; i<sum.size(); i++)
            for (int j=0; j<sum.size(); j++)
                sum[i][j] = mat.at(i)[j] + add_val;

        return sum;
    }

    //....................Matrix - calculation type value
    Matrix operator - (CalcType sub_val)
    {
        Matrix sub_result(mat.size(), mat.at(0).size());
        
        for (int i=0; i<sub_result.size(); i++)
            for (int j=0; j<sub_result.size(); j++)
                sub_result[i][j] = mat.at(i)[j] - sub_val;

        return sub_result;
    }
    
    //....................-Matrix
    Matrix operator - ()
    {
        Matrix minus_mat(mat.size(), mat.at(0).size());
        
        for (int i=0; i<minus_mat.size(); i++)
            for (int j=0; j<minus_mat.size(); j++)
                minus_mat[i][j] = -1 * mat.at(i)[j];
        
        return minus_mat;
    }

    //....................Matrix * calculation type value
    Matrix operator * (CalcType com_val)
    {
        Matrix composition(mat.size(), mat.at(0).size());
        
        for (int i=0; i<composition.size(); i++)
            for (int j=0; j<composition.size(); j++)
                composition[i][j] = mat.at(i)[j] * com_val;

        return composition;
    }

    //....................Matrix / calculation type value
    Matrix operator / (CalcType div_val)
    {
        Matrix quotient(mat.size(), mat.at(0).size());
        
        for (int i=0; i<quotient.size(); i++)
            for (int j=0; j<quotient.size(); j++)
                quotient[i][j] = mat.at(i)[j] / div_val;

        return quotient;
    }

    friend ostream& operator << (ostream& out, const Matrix &m)
	{
        for (int i=0; i<m.mat.size(); i++)
            out << m.mat.at(i);
		return out;
	}
};


//____________________________VECTOR OF MATRIX OF CALCULATION TYPE_____________________________
class VectorOfMatrix
{
    vector<Matrix> vecofmat;

public:
    VectorOfMatrix(const int n=0) { vecofmat.resize(n); }

    //....................Returns the size of Vector of Matrix
    int size() { return vecofmat.size(); }

    //....................Vector of Matrix indexing operator
    Matrix& operator () (const int index) { return vecofmat.at(index); }
};


//____________________________MATRIX OF MATRIX OF CALCULATION TYPE_____________________________
class MatrixOfMatrix
{
    vector<vector<Matrix>> matofmat;

public:
    MatrixOfMatrix(const int size=0)
    {
        matofmat.resize(size);
        for (int i=0; i<size; i++)
            matofmat.at(i).resize(size);
    }

    //....................Returns numbers of rows & columns in Matrix of Matrix
    int size() { return matofmat.size(); }

    //....................Matrix of Matrix indexing operator
    Matrix& operator () (const int row_index, const int col_index) { return matofmat.at(row_index).at(col_index); }
};


//______________________FUNCTION OF CALCULATION TYPE ON 1D, 2D OR 3D SPACE_____________________
class Function
{
    Vector x;
    Vector y;

private:
    static bool SortByX(const pair<CalcType, CalcType> &a, const pair<CalcType, CalcType> &b)
    { return (a.first < b.first); }

public:
    Function(const int size=0)
    {
        x.resize(size);
        y.resize(size);
    }

    Function(Vector x_grid, Vector y_values) : x(x_grid), y(y_values) {}

    Function(Vector x_grid, Matrix y_values)
    {
        vector<pair<CalcType, CalcType>> x_and_y;
        for (int i=0; i<x_grid.size(); i++)
        {
            for (int j=i; j<x_grid.size(); j++)
            {
                assert (std::abs(y_values[i][j] - y_values[j][i]) < 1e-08);
                CalcType x_hypot = std::hypot(x_grid[i], x_grid[j]);
                x_and_y.push_back(make_pair(x_hypot, y_values[i][j]));
            }
        }
        sort(x_and_y.begin(), x_and_y.end(), SortByX);
        for (auto pair : x_and_y)
        {
            x.add(pair.first);
            y.add(pair.second);
        }
    }

    Function(Vector x_grid, Matrix y_values, CalcType max_x_value)
    {
        vector<pair<CalcType, CalcType>> x_and_y;
        for (int i=0; i<x_grid.size(); i++)
        {
            for (int j=i; j<x_grid.size(); j++)
            {
                assert (std::abs(y_values[i][j] - y_values[j][i]) < 1e-08);
                CalcType x_hypot = std::hypot(x_grid[i], x_grid[j]);
                x_and_y.push_back(make_pair(x_hypot, y_values[i][j]));
            }
        }
        sort(x_and_y.begin(), x_and_y.end(), SortByX);
        
        int i = 0; 
        while (x_and_y.at(i).first <= max_x_value)
        {
            x.add(x_and_y.at(i).first);
            y.add(x_and_y.at(i).second);
            i++;
        }
    }

    int size() { return x.size(); }

    Vector X() { return x; }

    Vector Y() { return y; }

    static Function trim(Function func, CalcType final_x_value)
    {
        auto x_values = func.X().get_vector();
        auto it = std::find(x_values.begin(), x_values.end(), final_x_value);
        assert (it != x_values.end());
        int final_index = std::distance(x_values.begin(), it);
        
        Vector x_values_trimmed;
        Vector y_values_trimmed;
        for (int i=0; i<=final_index; i++)
        {
            x_values_trimmed.add(func.X()[i]);
            y_values_trimmed.add(func.Y()[i]);
        }

        Function func_trimmed(x_values_trimmed, y_values_trimmed);
        return func_trimmed;
    }

    CalcType findY(CalcType x_value)
    {
        assert ((x_value >= x[0]) && (x_value <= x.last()));
        
        auto it = std::find(x.get_vector().begin(), x.get_vector().end(), x_value);
        assert (it != x.get_vector().end());
        
        int index = std::distance(x.get_vector().begin(), it);
        return y[index];
    }

    pair<bool, CalcType> ConstantSign()
    {
        bool is_constant = true;
        int sign_init = sgn(y[1] - y[0]);
        CalcType val = 0;
        for (int i=2; i<y.size(); i++)
        {
            int sign = sgn(y[i] - y[i-1]);
            if ((sign != sign_init) && (sign != 0))
            {
                is_constant = false;
                val = x[i];
            }
        }
        return make_pair(is_constant, val);
    }

    CalcType Interpolate(CalcType x_value)
    {
        int i = 1;
        while (x_value > x[i])
            i++;

        //cout << "X[i-1] = " << x[i-1] << "   X[i] = " << x[i] << "\nY[i-1] = " << y[i-1] << "   Y[i] = " << y[i] << "\n\n";
        
        CalcType weight = (x_value - x[i-1]) / (x[i] - x[i-1]);
        CalcType y_value = y[i-1] + weight * (y[i] - y[i-1]);
        return y_value;
    }

    CalcType Interpolate(CalcType x_value, CalcType asymptotic_value)
    {
        if (x_value > x.last())
            return asymptotic_value;

        int i = 1;
        while (x_value > x[i])
            i++;
        
        CalcType weight = (x_value - x[i-1]) / (x[i] - x[i-1]);
        CalcType y_value = y[i-1] + weight * (y[i] - y[i-1]);
        return y_value;
    }

    CalcType IntegralRectangles()
    {
        assert (x.size() % 2 == 1);

        CalcType integral = 0;
        CalcType step = (2 * (x.last() - x[0])) / x.size();
        for (int i=1-x.size(); i<x.size()-1; i+=2)
            integral += y[abs(i + 1)];
        
        assert (integral != 0);
        integral *= step;
        return integral;
    }

    //....................All Y values + calculation type value
    Function operator + (CalcType add_val)
    {
        Vector new_y;

        for (int i=0; i<y.size(); i++)
            new_y.add(y[i] + add_val);

        Function new_func(x, new_y);
        return new_func;
    }

    //....................Vector - calculation type value
    Function operator - (CalcType sub_val)
    {
        Vector new_y;

        for (int i=0; i<y.size(); i++)
            new_y.add(y[i] - sub_val);

        Function new_func(x, new_y);
        return new_func;
    }
    
    //....................-Vector
    Function operator - ()
    {
        Vector new_y;

        for (int i=0; i<y.size(); i++)
            new_y.add(-1 * y[i]);

        Function new_func(x, new_y);
        return new_func;
    }

    //....................Vector * calculation type value
    Function operator * (CalcType com_val)
    {
        Vector new_y;

        for (int i=0; i<y.size(); i++)
            new_y.add(y[i] * com_val);

        Function new_func(x, new_y);
        return new_func;
    }

    //....................Vector / calculation type value
    Function operator / (CalcType div_val)
    {
        Vector new_y;

        for (int i=0; i<y.size(); i++)
            new_y.add(y[i] / div_val);

        Function new_func(x, new_y);
        return new_func;
    }

    friend ostream& operator << (ostream& out, Function &f)
	{
		for(int i=0; i<f.size(); i++)
			out << "X: " << f.x[i] << "   Y: " << f.y[i] << '\n';
		return out;
	}
};


//_________________VECTOR OF FUNCTION OF CALCULATION TYPE ON 1D, 2D OR 3D SPACE________________
class VectorOfFunction
{
    vector<Function> vecoffunc;

public:
    VectorOfFunction(int size=0) { vecoffunc.resize(size); }

    //....................Returns the size of Vector of Function
    int size() { return vecoffunc.size(); }

    //....................Vector of Matrix indexing operator
    Function& operator () (const int index) { return vecoffunc.at(index); }
};


//_____________SQUARE MATRIX OF FUNCTION OF CALCULATION TYPE ON 1D, 2D OR 3D SPACE_____________
class MatrixOfFunction
{
    vector<vector<Function>> matoffunc;

public:
    MatrixOfFunction(int size=0)
    {
        matoffunc.resize(size);
        for (int i=0; i<size; i++)
            matoffunc.at(i).resize(size);
    }

    //....................Returns numbers of rows & columns in Matrix of Function
    int size() { return matoffunc.size(); }

    //....................Matrix of Matrix indexing operator
    Function& operator () (int row_index, int col_index) { return matoffunc.at(row_index).at(col_index); }
};


//__________________________________________________CONVERTION FUNCTIONS__________________________________________________
Vector MatrixToVector(Matrix mat, bool symmetric=false)
{
    Vector result;
    for (int i=0; i<mat.size(); i++)
        for (int j=symmetric?i:0; j<mat.size(); j++)
            result.add(mat[i][j]);
    return result;
}

Vector MatrixOfMatrixToVector(MatrixOfMatrix mat_of_mat, bool symmetric=false)
{
    Vector result;
    for (int i=0; i<mat_of_mat.size(); i++)
        for (int j=0; j<mat_of_mat.size(); j++)
            for (int k=0; k<mat_of_mat(i, j).size(); k++)
                for (int l=symmetric?k:0; l<mat_of_mat(i, j).size(); l++)
                    result.add(mat_of_mat(i, j)[k][l]);
    return result;
}

Matrix VectorToMatrix(Vector vec, int mat_size, bool symmetric=false)
{
    int size = mat_size * mat_size;
    if (symmetric)
        for (int i=1; i<mat_size; i++)
            size -= i;
    assert (vec.size() == size);
    
    Matrix result(mat_size);
    int index = 0;

    for (int i=0; i<result.size(); i++)
    {
        for (int j=symmetric?i:0; j<result.size(); j++)
        {
            result[i][j] = vec[index];
            if (symmetric)
                result[j][i] = vec[index];
            index++;
        }
    }
    return result;
}

MatrixOfMatrix VectorToMatrixOfMatrix(Vector vec, int mat_of_mat_size, int mat_size, bool symmetric=false)
{
    int size = mat_size * mat_size;
    if (symmetric)
        for (int i=1; i<mat_size; i++)
            size -= i;
    size *= mat_of_mat_size * mat_of_mat_size;
    assert (vec.size() == size);
    
    MatrixOfMatrix result(mat_of_mat_size);
    int index = 0;

    for (int i=0; i<result.size(); i++)
    {
        for (int j=0; j<result.size(); j++)
        {
            Matrix mat(mat_size);
            for (int k=0; k<mat.size(); k++)
            {
                for (int l=symmetric?k:0; l<mat.size(); l++)
                {
                    mat[k][l] = vec[index];
                    if (symmetric)
                        mat[l][k] = vec[index];
                    index++;
                }
            }
            result(i, j) = mat;
        }
    }
    return result;
}

MatrixOfFunction VectorToMatrixOfFunction(Vector x_vec, Vector y_vec, int mat_of_func_size)
{
    assert (y_vec.size() == mat_of_func_size * mat_of_func_size * x_vec.size());
    MatrixOfFunction result(mat_of_func_size);
    int index = 0;

    for (int i=0; i<result.size(); i++)
    {
        for (int j=0; j<result.size(); j++)
        {
            Vector y_values;

            for (int k=0; k<x_vec.size(); k++)
            {
                y_values.add(y_vec[index]);
                index++;
            }
            Function new_func(x_vec, y_values);
            result(i, j) = new_func;
        }
    }
    return result;
}

bool IsSymmetric(Matrix mat)
{
    bool is_symmetric = true;
    for (int i=0; i<mat.size(); i++)
        for (int j=0; j<mat.size(); j++)
            if (mat[i][j] != mat[j][i])
                is_symmetric = false;
    return is_symmetric;
}


//________________________________________________TIME FROM CLOCK_T_________________________________________________
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

        str_seconds = seconds / 10 == 0 ? "0" + to_string(std::abs(seconds)) : to_string(std::abs(seconds));
        str_minutes = minutes / 10 == 0 ? "0" + to_string(std::abs(minutes)) : to_string(std::abs(minutes));
        str_hours = to_string(std::abs(hours));
    }

    friend ostream& operator << (ostream& out, const Time &t)
	{
        out << t.str_hours << ':' << t.str_minutes << ':' << t.str_seconds;
		return out;
	}
};