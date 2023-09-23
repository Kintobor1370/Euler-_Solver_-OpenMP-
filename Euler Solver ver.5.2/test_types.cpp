#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <utility> 
#include <cassert>
#include <tuple>
#include <algorithm>

#include "types.h"

using namespace std;

int main()
{
    //....................VECTOR TEST
    vector<CalcType> test_vec {1, 2, 3, 4, 5};
    Vector TestVec(test_vec);
    assert (TestVec.get_vector() == test_vec);
    
    // Indexing operator
    for (int i=0; i<TestVec.size(); i++)
        assert (TestVec[i] == i+1);

    // Vector is not empty
    assert (!TestVec.empty());

    // Vector's size
    assert (TestVec.size() == 5);

    // Resize to bigger size
    TestVec.resize(8);
    for (int i=0; i<5; i++)
        assert (TestVec[i] == i+1);
    for (int i=5; i<TestVec.size(); i++)
        assert (TestVec[i] == 0);

    // Resize to smaller size
    TestVec.resize(3);
    for (int i=0; i<TestVec.size(); i++)
        assert (TestVec[i] == i+1);

    // Clear Vector
    TestVec.clear();
    assert (TestVec.empty());
    assert (TestVec.size() == 0);
    

    //....................SQUARE MATRIX TEST
    Vector test_mat(vector<CalcType> { 
        1, 2, 3, 4, 5,
        6, 7, 8, 9, 10,
        11, 12, 13, 14, 15,
        16, 17, 18, 19, 20,
        21, 22, 23, 24, 25
    });
    Matrix TestMat = VectorToMatrix(test_mat, 5);
    
    CalcType arg = 2.54;

    Matrix TestOperation = TestMat + arg;
    for (int i=0; i<TestOperation.size(); i++)
        for (int j=0; j<TestOperation.size(); j++)
            assert (TestOperation[i][j] == TestMat[i][j] + arg);

    TestOperation = TestMat - arg;
    for (int i=0; i<TestOperation.size(); i++)
        for (int j=0; j<TestOperation.size(); j++)
            assert (TestOperation[i][j] == TestMat[i][j] - arg);

    TestOperation = -TestMat;
    for (int i=0; i<TestOperation.size(); i++)
        for (int j=0; j<TestOperation.size(); j++)
            assert (TestOperation[i][j] == TestMat[i][j] * -1);

    TestOperation = TestMat * arg;
    for (int i=0; i<TestOperation.size(); i++)
        for (int j=0; j<TestOperation.size(); j++)
            assert (TestOperation[i][j] == TestMat[i][j] * arg);

    TestOperation = TestMat / arg;
    for (int i=0; i<TestOperation.size(); i++)
        for (int j=0; j<TestOperation.size(); j++)
            assert (TestOperation[i][j] == TestMat[i][j] / arg);
    
    cout << "\n=========================ALL TESTS COMPLETE!=========================\n\n";
    return 0;
}