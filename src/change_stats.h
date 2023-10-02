#pragma once
#include <algorithm>
#include <ctime>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include "graph_class.h"
using namespace std;




double cs_order2(int i, int j, mlnet& mlnetwork, int k, int l) {
    if (k != l) {
        if (mlnetwork.is_edge(i, j, l)) {
            return 1;
        }
        else return 0;
    }
    else return 1;
    
}

double cs_order3(int i, int j, mlnet& mlnetwork, int k, int l) {
    if (k != l) {
        if (mlnetwork.is_edge(i, j, l)) {
            return 1;
        }
        else return 0;
    }
    else return 1;

}
