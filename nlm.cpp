#include "nlm.h"


NLM::NLM(vtkImageData *imageData) {/*{{{*/
    this->imageData = imageData;
    imageData->GetDimensions(dim);

    mov[0][0] = -1, mov[0][1] = -1;
    mov[1][0] = -1, mov[1][1] = 0;
    mov[2][0] = -1, mov[2][1] = 1;
    mov[3][0] = 0, mov[3][1] = 1;
    mov[4][0] = 1, mov[4][1] = 1;
    mov[5][0] = 1, mov[5][1] = 0;
    mov[6][0] = 1, mov[6][1] = -1;
    mov[7][0] = 0, mov[7][1] = -1;
    
    outmov[0][0] = -2, outmov[0][1] = -2;
    outmov[1][0] = -2, outmov[1][1] = -1;
    outmov[2][0] = -2, outmov[2][1] = 0;
    outmov[3][0] = -2, outmov[3][1] = 1;
    outmov[4][0] = -2, outmov[4][1] = 2;
    outmov[5][0] = -1, outmov[5][1] = 2;
    outmov[6][0] = 0, outmov[6][1] = 2;
    outmov[7][0] = 1, outmov[7][1] = 2;
    outmov[8][0] = 2, outmov[8][1] = 2;
    outmov[9][0] = 2, outmov[9][1] = 1;
    outmov[10][0] = 2, outmov[10][1] = 0;
    outmov[11][0] = 2, outmov[11][1] = -1;
    outmov[12][0] = 2, outmov[12][1] = -2;
    outmov[13][0] = 1, outmov[13][1] = -2;
    outmov[14][0] = 0, outmov[14][1] = -2;
    outmov[15][0] = -1, outmov[15][1] = -2;

    /*imp noise*/
    mp_noise[2][0] = 7;
    mp_noise[2][1] = 9;
    mp_noise[2][2] = 11;
    for(int i = 3; i <= 10; i ++)
        mp_noise[2][3] = 13;

    /*gau noise*/
    //UNDO
}/*}}}*/

NLM::~NLM() {
}

void NLM::setPatch(int patch) {
    this->patch = patch;
    extendPatch = patch / 2;           
}

void NLM::setWindow(int window) {
    this->window = window;
    extend = window / 2;
}

void NLM::setDelta(double s, double i, double m, double sm, double j, double sm2) {
    delta_S = s;
    delta_I = i;
    delta_M = m;
    delta_SM = sm;
    delta_J = j;
    delta_SM_2 = sm;
}

void NLM::extendEdge() {/*{{{*/
    initData(dim[0] + extend * 2, dim[1] + extend * 2);

    /*copy imageData to data*/
    for(int i = 0; i < dim[0]; i ++)
        for(int j = 0; j < dim[1]; j ++) {
            int *tmp0 =  static_cast<int *>(imageData->GetScalarPointer(i, j, 0));
            data[extend + i][extend + j] = *tmp0;
        }
    
    /*extend the edge for the y*/
    for(int i = 0; i < dim[0]; i ++)
        for(int j = 1; j <= extend; j ++) {
            data[i][extend - j] = data[i][extend + j];
            data[i][extend + dim[1] + j - 1] = data[i][extend + dim[1] - j - 1];
        }

    /*extend the edge for the x*/
    for(int i = 1; i <= extend; i ++)
        for(int j = 0; j < dim[1] + extend * 2; j ++) {
           data[extend - i][j] = data[extend + i][j];
           data[extend + dim[0] + i - 1][j] = data[extend + dim[0] - i - 1][j];
       }
}/*}}}*/

/*remember the arrary's dimension*/
void NLM::initData(int x, int y) {/*{{{*/
    NLM::Node now;
    vector<int> tmp(y, 0);
    vector<NLM::Node> node_tmp(y, now);
    for(int i = 0; i < x; i ++) {
        data.push_back(tmp);
        ROAD.push_back(tmp);
        node.push_back(node_tmp);
    }
}/*}}}*/

void NLM::pretreatment(int block_nr, double &final_fangcha) {/*{{{*/
    int patch_board_beg = extendPatch;
    int patch_board_end = 2 * extend + dim[0] - extendPatch; //dim[0] == dim[1]
    /*Calculate ROAD*/
    for(int i = patch_board_beg; i < patch_board_end; i ++) {
        for(int j = patch_board_beg; j < patch_board_end; j ++) {
            ROAD[i][j] = calcuROAD(i, j, block_nr);
            node[i][j] = calcuNode(i, j);
            final_fangcha += node[i][j].fangcha;
            //calcuWinFc(i, j);
        }
        //printf("pretreatment:%d\n", i);
    }
    final_fangcha /= pow(double(patch_board_end - patch_board_beg), 2);
    puts("After calculate the ROAD");

    if(ROTATE) {
        for(int i = 0; i < extend * 2 + dim[0]; i ++)
            for(int j = 0; j < extend * 2 + dim[1]; j ++) {
                calcuSitaAndM(i, j);
            }
        puts("After calculate the Node");

        for(int i = patch_board_beg; i < patch_board_end; i ++) {
            for(int j = patch_board_beg; j < patch_board_end; j ++) {
                node[i][j].angle = calcuAngle(i, j);
            }
        }
        puts("After calculate the angle");
    }
}/*}}}*/

int NLM::calcuROAD(int x, int y, int block_nr) {/*{{{*/
    /*sort from less to greater*/
    priority_queue<int, vector<int>, greater<int> > q;
    while(!q.empty()) q.pop();
    //int sz = 2;
    int sz = 1;
    for(int i = x - sz; i <= x + sz; i ++)
        for(int j = y - sz; j <= y + sz; j ++) {
            if(i == x && j == y)
                continue;
            int tmp = abs(data[i][j] - data[x][y]);
            q.push(tmp);
        }

    int sum = 0;
    //int nr = (sz * 2 + 1) * (sz * 2 + 1) / 2;
    int nr = block_nr;
    for(int i = 0; i < nr; i ++) {
        sum += q.top();
        q.pop();
    }
    return sum;
}/*}}}*/

NLM::Node NLM::calcuNode(int x, int y) {/*{{{*/
    Node now;
    double sum = 0, sum2 = 0;
    vector<int> vec;
    for(int i = -extendPatch; i <= extendPatch; i ++)
        for(int j = -extendPatch; j <= extendPatch; j ++) {
            int val = data[x + i][y + j];
            sum += val;
            sum2 += val * val;
            vec.push_back(val);
        }
    sort(vec.begin(), vec.end());

    double ave = sum / vec.size();
    double fc = 0;
    for(int i = 0; i < vec.size(); i ++)
        fc += (vec[i] - ave) * (vec[i] - ave);

    now.jun = ave;
    now.zhong = vec[vec.size() / 2];
    now.fangcha = fc / vec.size();
    now.er = sum2 / vec.size();
    
    //cout << now.fangcha << endl;
    return now;
}/*}}}*/

void NLM::calcuWinFc(int x, int y) {
    int step = 11;
    int right_board = 2 * extend + dim[0]; //dim[0] == dim[1]
    int tmp = step / 2;
    while(step <= 11) {
        int tmp = step / 2;
        if(x < tmp || x >= right_board - tmp || \
                y < tmp || y >= right_board - tmp)
            break;
        int sum = 0;
        int beg_i = x - tmp, beg_j = y - tmp;
        int end_i = x + tmp, end_j = y + tmp;
        for(int i = beg_i; i <= end_i; i ++)
            for(int j = beg_j; j <= end_j; j ++)
                sum += data[i][j];
        double ave = sum * 1.0 / (step * step);
        double fc = 0;
        for(int i = beg_i; i <= end_i; i ++)
            for(int j = beg_j; j <= end_j; j ++)
                fc += (data[i][j] - ave) * (data[i][j] - ave);
        node[x][y].win_fc[step] = (int)(fc / step / step);

        //cout << "step:" << node[x][y].win_fc[step] << endl;
        step += 2;
    }
}

void NLM::calcuSitaAndM(int x, int y) {/*{{{*/
    double cha1 = 0, cha2 = 0;
    if(x + 1 < 2 * extend + dim[0] && x > 0 && y > 0 && y + 1 < 2 * extend + dim[1]) {
        cha1 = data[x + 1][y] - data[x - 1][y];
        cha2 = data[x][y + 1] - data[x][y - 1];
    }
    double m = sqrt(cha1 * cha1 + cha2 * cha2);
    double sita = atan2(cha2, cha1) * 180 / PI + 180;
    sita = (sita >= 360 ? sita - 360 : sita);

    node[x][y].m = m;
    node[x][y].sita = sita / 45;

}/*}}}*/

int NLM::calcuAngle(int x, int y) {/*{{{*/
    vector<double> his(8, 0);
    for(int ii = -extendPatch; ii <= extendPatch; ii ++)
        for(int jj = -extendPatch; jj <= extendPatch; jj ++)
            his[node[x + ii][y + jj].sita] += calcu_W_SM_rotate(ii, jj, 0, 0) * node[x + ii][y + jj].m; 

    int Max = -1, mark = 0;
    for(int i = 0; i < his.size(); i ++) {
        if(his[i] > Max) {
            Max = his[i];
            mark = i;
        }
    }
    return mark;
}/*}}}*/

int NLM::judge_window(int x, int y) {
    int dft = 11;
    //cout << "m" << node[x][y].win_fc[11] << endl;
    switch(noise) {
        case 1: break; //UNDO
        case 2: return mp_noise[2][node[x][y].win_fc[dft] / 1000];
        case 3: break; //UNDO
        default: break;
    }
    return dft;
}

double NLM::calcuW_SIM(int x, int y) {
    /*
    window = judge_window(x, y);
    cout << "w:" << window << endl;
    int half = window / 2;
    int right_board = 2 * extend + dim[0]; //dim[0] == dim[1]
    */
    int disDiff = extend - extendPatch;
    /*
    if(x - half >= extendPatch && x + half < right_board - extendPatch - 1 && \
        y - half >= extendPatch && y + half < right_board - extendPatch - 1)
        disDiff = half - extendPatch; 
    */
    double bigDenominator = 0, bigNumerator = 0;
    
    for(int i = x - disDiff; i <= x + disDiff; i ++)
        for(int j = y - disDiff; j <= y + disDiff; j ++) {
            if(ACCERLERATE && judge(x, y, i, j))
                continue;

            double sumDenominator = 0, sumNumerator = 0;
            for(int ii = -extendPatch; ii <= extendPatch; ii ++)
                for(int jj = -extendPatch; jj <= extendPatch; jj ++) {
                    double val = 0;
                    double tmp2 = 1;
                    val = data[i + ii][j + jj] - data[x + ii][y + jj];
                    sumDenominator += tmp2;
                    sumNumerator += tmp2 * val * val;
                }

            double tmp = calcu_W_M(sumNumerator, sumDenominator);

            bigDenominator += tmp;
            if(UNBIAS)
                bigNumerator += tmp * data[i][j] * data[i][j];
            else
                bigNumerator += tmp * data[i][j];
        }

    return bigDenominator == 0 ? 0 : bigNumerator / bigDenominator;
}

void NLM::calcuW() {
    imageData_new = vtkImageData::New();
    imageData_new->DeepCopy(imageData);
    for(int i = extend; i < extend + dim[0]; i ++) {
        for(int j = extend; j < extend + dim[1]; j ++) {
            int tmp = calcuW_SIM(i, j);
            //tmp = (tmp < 0 ? 0 : tmp);
            //tmp = (tmp >= upBoard ? upBoard - 1: tmp);
            //printf("~%d %d %d\n", i, j, tmp);
            int *tt = static_cast<int *>(imageData_new->GetScalarPointer(i - extend, j - extend, 0));
            *tt = (int)tmp;
        }
        //printf("row:%d\n", i);
    }
}

double NLM::calcu_W_SM_rotate(int x1, int y1, int x2, int y2) {/*{{{*/
    double tmp = dis2(x1, y1, x2, y2) / 2 / delta_SM_2 / delta_SM_2;
    return pow(E, -tmp);
}

double NLM::calcu_W_S(int x1, int y1, int x2, int y2) {
    double tmp = dis2(x1, y1, x2, y2) / 2 / delta_S / delta_S;
    //double tmp2 = max(abs(x1 - x2), abs(y1 - y2));
    //double tmp = tmp2 * tmp2 / 2 / delta_S / delta_S;
    return pow(E, -tmp);
}

double NLM::calcu_W_I(int x, int y) {
    double tmp = ROAD[x][y] * ROAD[x][y] / 2 / delta_I / delta_I;
    return pow(E, -tmp);
}

double NLM::calcu_J_I(int x1, int y1, int x2, int y2) {
    double tmp = (ROAD[x1][y1] + ROAD[x2][y2]) / 2;
    double tmp2 = tmp * tmp / 2 / delta_J / delta_J;
    return pow(E, -tmp2);
}

double NLM::calcu_W_M(double x, double y) {
    double tmp = (y == 0 ? INF : x / y);
    double tmp2 = tmp / 2 / delta_M / delta_M;
    return pow(E, -tmp2);
}

double NLM::dis2(int x1, int y1, int x2, int y2) {
    return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

vtkImageData *NLM::getImageData() {
    return imageData_new;
}/*}}}*/

void NLM::setAccelerate(double jun, double zhong, double fangcha, double er) {/*{{{*/
    this->JUN = jun;
    this->ZHONG = zhong;
    this->FANGCHA = fangcha;
    this->ER = er;
}/*}}}*/

bool NLM::judge(int x, int y, int a, int b) {/*{{{*/
    if (node[a][b].jun == 0)
        return 0;
    double cmp_jun = node[x][y].jun / node[a][b].jun;
    if(cmp_jun < JUN || cmp_jun > 2 - JUN)
        return 0;

    if (node[a][b].zhong == 0)
        return 0;
    double cmp_zhong = node[x][y].zhong / node[a][b].zhong;
    if(cmp_zhong < ZHONG || cmp_zhong > 2 - ZHONG)
        return 0;
    
    if (node[a][b].fangcha == 0)
        return 0;
    double cmp_fangcha = node[x][y].fangcha / node[a][b].fangcha;
    if(cmp_fangcha < FANGCHA || cmp_fangcha > 2 - FANGCHA)
        return 0;

    /*
    if (node[a][b].er == 0)
        return 0;
    double cmp_er = node[x][y].er / node[a][b].er;
    if(cmp_er < ER || cmp_er > 2 - ER)
        return 0;
    */
    return 1;
}/*}}}*/

void NLM::unbias(double level_p) {
    level_p = 2 * level_p * level_p;
    for(int i = 0; i < dim[0]; i ++)
        for(int j = 0; j < dim[1]; j ++) {
            int *val = static_cast<int *>(imageData_new->GetScalarPointer(i, j, 0));
            *val = (int)sqrt(max(0.0, *val - level_p));
        }
}

void NLM::cutIllegal() {
    for(int i = 0; i < dim[0]; i ++)
        for(int j = 0; j < dim[1]; j ++) {
            int *val = static_cast<int *>(imageData_new->GetScalarPointer(i, j, 0));
            int tmp = *val;
            tmp = (tmp < 0 ? 0 : tmp);
            tmp = (tmp >= upBoard ? upBoard - 1: tmp);
            *val = tmp;
            //printf("~%d %d %d\n", i, j, tmp);
        }
}
