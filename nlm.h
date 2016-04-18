#ifndef NLM_H
#define NLM_H

#include "vtkImageData.h"
#include <queue>
#include <cmath>
#include <map>
#include <string>
#include <vector>

using namespace std;

const int SIZE = 550;
const int MSIZE = 150;
const int KSIZE = 20;
#define E 2.71828182845904523536
#define PI 3.141592653

#define ACCERLERATE 0
#define ROTATE 0
#define UNBIAS 0

class NLM {
public:
    NLM(vtkImageData *imageData);
    ~NLM();
    void setPatch(int);
    void setWindow(int);
    void setDelta(double s, double i, double m, double sm, double j, double sm2);
    vtkImageData *getImageData();

    /*Store the imageData in data(vecotr) and extend the edge as 0*/
    void extendEdge();

    /*Calculate the ROAD and Euclidean distance*/
    void pretreatment(int, double &);

    /*main function, Calculate the value of w*/
    void calcuW();

    /************************/
    void setAccelerate(double jun, double zhong, double fangcha, double er);
    
    int noise; //0:none, 1:gau, 2:imp, 3:mix

    void unbias(double level);

    void cutIllegal();
private:
    static const int upBoard = 256;
    static const int INF = 1 << 30;
    vtkImageData *imageData, *imageData_new;
    int dim[2];
    int patch, window;
    
    /*The ROAD value*/
    //double ROAD[SIZE][SIZE][MSIZE];
    vector<vector<int> > ROAD;

    /*Calculate ROAD*/
    int calcuROAD(int x, int y, int block_nr);

    /*main function, Calculate the value of w_S, w_I, w_M*/
    double calcuW_SIM(int, int);

    /*Store the imagedata in data(vector)*/
    void initData(int x, int y);
    
    /*store imagedata in vector*/
    vector<vector<int> > data;

    /*equal to window / 2*/
    int extend;

    /*equal to patch / 2*/
    int extendPatch;

    double delta_S, delta_I, delta_M, delta_SM, delta_J, delta_SM_2;

    double calcu_W_SM(int, int, int, int);
    double calcu_W_SM_rotate(int, int, int, int);
    double calcu_W_S(int, int, int, int);
    double calcu_W_I(int, int);
    double calcu_J_I(int, int, int, int);
    double calcu_W_M(double, double);

    /*distance * distance*/
    double dis2(int, int, int, int);

    /******************************/
    
    int mov[8][2];
    int outmov[16][2];

    static const int mov_step = 8;
    static const int outmov_step = 8;

    bool judge(int x, int y, int a, int b);

    struct Node {
        double jun, zhong, fangcha, er;
        double m;
        int sita, angle;
        int win_fc[24];
        Node() {
            memset(win_fc, 0, sizeof(win_fc));
        }
    };
    vector<vector<Node> > node;

    Node calcuNode(int x, int y);

    double JUN, ZHONG, FANGCHA, ER;

    void calcuSitaAndM(int, int);
    int calcuAngle(int, int);
    void calcuWinFc(int, int); //calculate variance
    int judge_window(int, int);

    map<int, int> mp_noise[4];
};

#endif //NLM_H
