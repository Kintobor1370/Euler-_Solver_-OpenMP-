/**********************************************************************************************
**                                      TYPES HEADER FILE                                    **
**  This file contains classes for mathematical data structures used in Solver.cpp program.  **
**********************************************************************************************/

#include <vector>
#include <numeric>
using namespace std;

//___________________________________MAIN CALCULATION TYPE_____________________________________
using CalcType = double;


//_________________________________VECTOR OF CALCULATION TYPE__________________________________
class Vector
{
    vector<CalcType> vec;

public:
    Vector(int n = 0) { vec.resize(n); }

    Vector(int n, CalcType val)
    {
        vec.resize(n);
        for (int i=0; i<n; i++)
            vec.at(i) = val;
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

    //....................Adds an element to the back of Vector
    void add(CalcType new_elem) { vec.push_back(new_elem); }

    //....................Adds a Vector to the back of Vector
    void add(Vector new_vec)
    {
        for (int i=0; i<new_vec.size(); i++)
            vec.push_back(new_vec[i]);
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
};


//_________________________________MATRIX OF CALCULATION TYPE__________________________________
class Matrix
{
    vector<Vector> mat;

public:
    Matrix(const int r = 0, const int c = 0) { resize(r, c); }

    Matrix(const int r, const int c, CalcType val)
    {
        resize(r, c);
        for (int i=0; i<r; i++)
            for (int j=0; j<c; j++)
                mat.at(i)[j] = val;
    }

    Vector& operator [] (const int row_num) { return mat.at(row_num); }

    int Rows() { return mat.size(); }

    int Cols() { return mat.at(0).size(); }

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

    //....................Returns an element of Matrix which has the highest value
    CalcType max_value()
    {
        CalcType max_val = mat.at(0)[0];

        for (int i=0; i<mat.size(); i++)
            for (int j=0; j<mat.at(0).size(); j++)
                if (mat.at(i)[j] > max_val)
                    max_val = mat.at(i)[j];
        return max_val;
    }

    //....................Matrix + calculation type value
    Matrix operator + (CalcType add_val)
    {
        Matrix sum(mat.size(), mat.at(0).size());

        for (int i=0; i<sum.Rows(); i++)
            for (int j=0; j<sum.Cols(); j++)
                sum[i][j] = mat.at(i)[j] + add_val;

        return sum;
    }

    //....................Matrix - calculation type value
    Matrix operator - (CalcType sub_val)
    {
        Matrix sub_result(mat.size(), mat.at(0).size());
        
        for (int i=0; i<sub_result.Rows(); i++)
            for (int j=0; j<sub_result.Cols(); j++)
                sub_result[i][j] = mat.at(i)[j] - sub_val;

        return sub_result;
    }
    
    //....................-Matrix
    Matrix operator - ()
    {
        Matrix minus_mat(mat.size(), mat.at(0).size());
        
        for (int i=0; i<minus_mat.Rows(); i++)
            for (int j=0; j<minus_mat.Cols(); j++)
                minus_mat[i][j] = -1 * mat.at(i)[j];
        
        return minus_mat;
    }

    //....................Matrix * calculation type value
    Matrix operator * (CalcType com_val)
    {
        Matrix composition(mat.size(), mat.at(0).size());
        
        for (int i=0; i<composition.Rows(); i++)
            for (int j=0; j<composition.Cols(); j++)
                composition[i][j] = mat.at(i)[j] * com_val;

        return composition;
    }

    //....................Matrix / calculation type value
    Matrix operator / (CalcType div_val)
    {
        Matrix quotient(mat.size(), mat.at(0).size());
        
        for (int i=0; i<quotient.Rows(); i++)
            for (int j=0; j<quotient.Cols(); j++)
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
    VectorOfMatrix(int n = 0) { vecofmat.resize(n); }

    //....................Vector of Matrix indexing operator
    Matrix& operator () (const int index) { return vecofmat.at(index); }
};


//____________________________MATRIX OF MATRIX OF CALCULATION TYPE_____________________________
class MatrixOfMatrix
{
    vector<vector<Matrix>> matofmat;

public:
    MatrixOfMatrix(const int r = 0, const int c = 0)
    {
        matofmat.resize(r);
        for (int i=0; i<r; i++)
            matofmat.at(i).resize(c);
    }

    //....................Matrix of Matrix indexing operator
    Matrix& operator () (const int row_index, const int col_index) { return matofmat.at(row_index).at(col_index); }

    //....................Returns numbers of rows in Matrix of Matrix
    int Rows() { return matofmat.size(); }

    //....................Returns numbers of columns in Matrix of Matrix
    int Cols() { return matofmat.at(0).size(); }
};