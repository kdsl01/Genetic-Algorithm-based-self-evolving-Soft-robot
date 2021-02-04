#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <stdarg.h>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <array>
#include <vector>

#define _USE_MATH_DEFINES

#include <cmath>
#include<chrono>

using namespace std;
using namespace std::chrono;
random_device rd;
mt19937 e(rd());
static const double Gravity = 9.81;
static const double dt = 0.0001;
static const double u_s = 0.61;
static const double u_k = 0.47;
static const double damping = 0.9999;
struct mass {
    double m;
    double p[3];
    double v[3];
    double a[3];
};

struct spring {
    double k;
    double L0;
    double L0x;
    int m1;
    int m2;
};

struct sets {
    int pn;
    int sn;
    vector<double> p0;
    vector<double> p1;
    vector<double> p2;
    vector<double> b;
    vector<double> c;
    vector<int> sck;
//    double p[7][3];
//    double b[28];
//    double c[28];
//    int sck[28];
    double speed;
};

bool compare(spring a, spring b) {
//for descending order replace with a.distancce >b.distancce
    if (a.L0 < b.L0)
        return 1;
    else
        return 0;
}

bool comparevar(sets a, sets b) {
//for descending order replace with a.distancce >b.distancce
    if (a.speed > b.speed)
        return 1;
    else
        return 0;
}

void copy_var(sets origin[], sets tops[]) {
    for (int i = 0; i < 10; i++) {
        int sn = 0;
        int pn = 0;
        sn = origin[i].sn;
        pn = origin[i].pn - 1;
        tops[i].sn = sn;
        tops[i].pn = pn + 1;
        tops[i].p0.resize(pn);
        tops[i].p1.resize(pn);
        tops[i].p2.resize(pn);
        tops[i].b.resize(sn);
        tops[i].c.resize(sn);
        tops[i].sck.resize(sn);
        for (int k = 0; k < pn; k++) {
            tops[i].p0[k] = origin[i].p0[k];
            tops[i].p1[k] = origin[i].p1[k];
            tops[i].p2[k] = origin[i].p2[k];
        }
        for (int j = 0; j < sn; j++) {
            tops[i].b[j] = origin[i].b[j];
            tops[i].c[j] = origin[i].c[j];
            tops[i].sck[j] = origin[i].sck[j];
        }
        tops[i].speed = 0;
    }
}


void random_gv3(sets var[]) {
    uniform_real_distribution<double> first(70, 100);
    uniform_real_distribution<double> cxt(0, 19);
    uniform_int_distribution<int> spc1(5, 20);
    uniform_int_distribution<int> pointnumber(6, 15);
    uniform_real_distribution<double> posxy(-0.1, 0.1);
    uniform_real_distribution<double> posz(0, 0.1);
    for (int count = 0; count < 20; count++) {
        int pn = 0;
        int sn = 0;
        pn = pointnumber(e) + 1;
        sn = (pn * (pn - 1)) / 2;
        var[count].pn = pn;
        var[count].sn = sn;
        var[count].p0.resize(pn - 1);
        var[count].p1.resize(pn - 1);
        var[count].p2.resize(pn - 1);
        var[count].b.resize(sn);
        var[count].c.resize(sn);
        var[count].sck.resize(sn);
        for (int i = 0; i < pn - 1; i++) {
            var[count].p0[i] = posxy(e);
            var[count].p1[i] = posxy(e);
            var[count].p2[i] = posz(e);
        }
        for (int i = 0; i < sn; i++) {
            var[count].b[i] = first(e);
            var[count].c[i] = cxt(e);
            var[count].sck[i] = spc1(e) * 1000;
        }
        var[count].speed = 0;
    }
}

void mass_initv3(vector<mass> &point, double weight, double sp[], sets var) {
    for (int i = 0; i < var.pn; i++) {
        point[i].m = weight;
        point[i].v[0] = 0;
        point[i].v[1] = 0;
        point[i].v[2] = 0;
        point[i].a[0] = 0;
        point[i].a[1] = 0;
        point[i].a[2] = 0;
    }
    point[0].p[0] = sp[0];
    point[0].p[1] = sp[1];
    point[0].p[2] = sp[2];
    for (int i = 1; i < var.pn; i++) {
        point[i].p[0] = sp[0] + var.p0[i - 1];
        point[i].p[1] = sp[1] + var.p1[i - 1];
        point[i].p[2] = sp[2] + var.p2[i - 1];
    }
}

void spring_initv3(vector<spring> &edge, vector<mass> &point, vector<int> &spc, int sn, int pn) {
    for (int count = 0; count < sn; count++) {
        edge[count].k = spc[count];
    }
    int k = 0;
    for (int i = 0; i < pn; i++) {
        for (int j = i + 1; j < pn; j++) {
            double distance = 0;
            edge[k].m1 = i;
            edge[k].m2 = j;
            distance = sqrt(pow((point[i].p[0] - point[j].p[0]), 2) + pow((point[i].p[1] - point[j].p[1]), 2) +
                            pow((point[i].p[2] - point[j].p[2]), 2));
            edge[k].L0 = distance;
            edge[k].L0x = distance;
            k++;
        }
    }
}

void breathv3(vector<spring> &edge, double T, vector<double> &b, vector<double> &c, int sn) {
    for (int i = 0; i < sn; i++) {
        edge[i].L0 = edge[i].L0x + sin(23 * T + c[i]) / b[i];
    }
}

double update_cube_firctionv3(vector<spring> &edge, vector<mass> &point, double cT, int pn, int sn) {
    vector<double> fx;
    vector<double> fy;
    vector<double> fz;
    vector<bool> groundflag;
    fx.resize(pn);
    fy.resize(pn);
    fz.resize(pn);
    groundflag.resize(pn);
    for (int count = 0; count < pn; count++) {
        fx[count] = 0;
        fy[count] = 0;
        fz[count] = 0;
        groundflag[count] = false;
    }
    for (int i = 0; i < pn; i++) {
        for (int j = 0; j < sn; j++) {
            if (edge[j].m1 == i || edge[j].m2 == i) {
                int self = 0;
                int opp = 0;
                double X = 0;
                double Y = 0;
                double Z = 0;
                double L = 0;
                double F = 0;
                if (edge[j].m1 == i) {
                    self = edge[j].m1;
                    opp = edge[j].m2;
                } else if (edge[j].m2 == i) {
                    self = edge[j].m2;
                    opp = edge[j].m1;
                }
                X = pow((point[opp].p[0] - point[self].p[0]), 2);
                Y = pow((point[opp].p[1] - point[self].p[1]), 2);
                Z = pow((point[opp].p[2] - point[self].p[2]), 2);
                L = sqrt(X + Y + Z);
                F = edge[j].k * (L - edge[j].L0);
                fx[i] = fx[i] - F * ((point[self].p[0] - point[opp].p[0]) / L);
                fy[i] = fy[i] - F * ((point[self].p[1] - point[opp].p[1]) / L);
                fz[i] = fz[i] - F * ((point[self].p[2] - point[opp].p[2]) / L);
            }
        }
        fz[i] = fz[i] - point[i].m * Gravity;
        if (point[i].p[2] <= 0) {
            double fh = sqrt(pow(fx[i], 2) + pow(fy[i], 2));
            if (fh < abs(fz[i] * u_s)) {
                groundflag[i] = true;
            } else {
                fx[i] = fx[i] - (fx[i] / fh) * abs(fz[i] * u_k);
                fy[i] = fy[i] - (fy[i] / fh) * abs(fz[i] * u_k);
            }
            fz[i] = fz[i] + (0 - point[i].p[2]) * 100000;
        }
    }
    for (int count = 0; count < pn; count++) {
        point[count].a[0] = fx[count] / point[count].m;
        point[count].a[1] = fy[count] / point[count].m;
        point[count].a[2] = fz[count] / point[count].m;
        point[count].v[0] = point[count].v[0] + point[count].a[0] * dt;
        point[count].v[1] = point[count].v[1] + point[count].a[1] * dt;
        point[count].v[2] = point[count].v[2] + point[count].a[2] * dt;
        point[count].v[0] = damping * point[count].v[0];
        point[count].v[1] = damping * point[count].v[1];
        point[count].v[2] = damping * point[count].v[2];
        if (groundflag[count] == true) {
            point[count].v[0] = 0;
            point[count].v[1] = 0;
            groundflag[count] = false;
        }
        point[count].p[0] = point[count].p[0] + point[count].v[0] * dt;
        point[count].p[1] = point[count].p[1] + point[count].v[1] * dt;
        point[count].p[2] = point[count].p[2] + point[count].v[2] * dt;
    }
    return cT + dt;
}


void mutv3(sets &Var) {
    int indexx1 = 0;
    int indexx2 = 0;
    int sindexx = 0;
    int pn = 0;
    int sn = 0;
    int temppn = 0;
    int tempsn = 0;
    uniform_int_distribution<int> massindex(0, Var.pn - 2);
    uniform_int_distribution<int> springindex(0, Var.sn - 1);
    uniform_real_distribution<double> first(70, 100);
    uniform_real_distribution<double> cxt(0, 19);
    uniform_int_distribution<int> spc1(5, 20);
    uniform_real_distribution<double> posxy(-0.1, 0.1);
    uniform_real_distribution<double> posz(0, 0.1);
    uniform_int_distribution<int> pointnumber(6, 15);
    indexx1 = massindex(e);
    indexx2 = massindex(e);
    sindexx = springindex(e);
    while (indexx1 == indexx2) {
        indexx2 = massindex(e);
    }
    Var.p0[indexx1] = posxy(e);
    Var.p1[indexx1] = posxy(e);
    Var.p2[indexx1] = posz(e);
    Var.p0[indexx2] = posxy(e);
    Var.p1[indexx2] = posxy(e);
    Var.p2[indexx2] = posz(e);
    Var.b[sindexx] = first(e);
    Var.c[sindexx] = cxt(e);
    Var.sck[sindexx] = spc1(e) * 1000;
    pn = pointnumber(e)+1;
    while(pn==Var.pn){
        pn = pointnumber(e)+1;
    }
    sn = (pn*(pn-1))/2;
    Var.p0.resize(pn-1);
    Var.p1.resize(pn-1);
    Var.p2.resize(pn-1);
    Var.b.resize(sn);
    Var.c.resize(sn);
    Var.sck.resize(sn);
    if(Var.pn<pn){
        temppn = Var.pn-1;
        tempsn = Var.sn;
        for(int i =temppn;i<pn-1;i++){
            Var.p0[i] = posxy(e);
            Var.p1[i] = posxy(e);
            Var.p2[i] = posz(e);
        }
        for (int j = tempsn; j < sn; j++) {
            Var.b[j] = first(e);
            Var.c[j] = cxt(e);
            Var.sck[j] = spc1(e) * 1000;
        }
        Var.pn = pn;
        Var.sn = sn;
    }else{
        Var.pn = pn;
        Var.sn = sn;
    }
}

void crossoverv3(sets tops[], sets children[]) {
    uniform_int_distribution<int> rate(0, 20);
    uniform_int_distribution<int> parents(0, 9);
    int cr = 0;
    int mr = 0;
    int cp1 = 0;
    int cp2 = 0;
    int father = 0;
    int mother = 0;
    int fpn;
    int mpn;
    int cpn;
    int lpn;
    int lspn;
    int avgpn;
    for (int i = 0; i < 20; i++) {
        cr = rate(e);
        if (cr == 0 || cr == 2 || cr == 4 || cr == 6 || cr == 8 || cr == 10 || cr == 12 || cr == 14 || cr == 16 ||
            cr == 18) {
            father = parents(e);
            mother = parents(e);
            while (father == mother) {
                mother = parents(e);
            }
            fpn = tops[father].pn;
            mpn = tops[mother].pn;
            if (fpn >= mpn) {
                cpn = mpn;
                lpn = fpn;
            } else {
                cpn = fpn;
                lpn = mpn;
            }
            lspn = tops[father].sn;
            uniform_int_distribution<int> crossposition(1, cpn - 2);
            children[i].sn = tops[father].sn;
            children[i].pn = tops[father].pn;
            children[i].p0.resize(tops[father].pn - 1);
            children[i].p1.resize(tops[father].pn - 1);
            children[i].p2.resize(tops[father].pn - 1);
            children[i].b.resize(lspn);
            children[i].c.resize(lspn);
            children[i].sck.resize(lspn);
            cp1 = crossposition(e);
            for (int ctx = 0; ctx < cp1; ctx++) {
                children[i].p0[ctx] = tops[mother].p0[ctx];
                children[i].p1[ctx] = tops[mother].p1[ctx];
                children[i].p2[ctx] = tops[mother].p2[ctx];
            }
            for (int ctx = cp1; ctx < fpn - 1; ctx++) {
                children[i].p0[ctx] = tops[father].p0[ctx];
                children[i].p1[ctx] = tops[father].p1[ctx];
                children[i].p2[ctx] = tops[father].p2[ctx];
            }
            for (int spn = 0; spn < lspn; spn++) {
                children[i].b[spn] = tops[father].b[spn];
                children[i].c[spn] = tops[father].c[spn];
                children[i].sck[spn] = tops[father].sck[spn];
            }
        } else {
            father = parents(e);
            children[i].pn = tops[father].pn;
            children[i].sn = tops[father].sn;
            children[i].p0.resize(tops[father].pn - 1);
            children[i].p1.resize(tops[father].pn - 1);
            children[i].p2.resize(tops[father].pn - 1);
            children[i].b.resize(tops[father].sn);
            children[i].c.resize(tops[father].sn);
            children[i].sck.resize(tops[father].sn);
            for (int ct = 0; ct < tops[father].pn - 1; ct++) {
                children[i].p0[ct] = tops[father].p0[ct];
                children[i].p1[ct] = tops[father].p1[ct];
                children[i].p2[ct] = tops[father].p2[ct];
            }
            for (int ctx = 0; ctx < tops[father].sn; ctx++) {
                children[i].b[ctx] = tops[father].b[ctx];
                children[i].c[ctx] = tops[father].c[ctx];
                children[i].sck[ctx] = tops[father].sck[ctx];
            }
        }
        children[i].speed = 0;
        mr = rate(e);
        if (mr == 1 || mr == 3 || mr == 5 || mr == 7 || mr == 9 || mr == 11 || mr == 13 || mr == 15 || mr == 17 ||
            mr == 19) {
            mutv3(children[i]);
        }
    }
}


int main() {
    vector<mass> point;
    vector<spring> edge;
    sets varx[20];
    sets tophalf[10];
    double sp[3];
    double bestspeed;
    sp[0] = 0.0;
    sp[1] = 0.0;
    sp[2] = 0.0;
    double ctm[3];
    ctm[0] = 0;
    ctm[1] = 0;
    ctm[2] = 0;
    ofstream fout2("GA_Robot_test1v5_2.csv", ios::out);
//    ofstream fout4("GA_robot_for_diversity_sp.csv",ios::out);
//    ofstream fout5("GA_robot_for_diversity_np.csv",ios::out);
    random_gv3(varx);
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for (int pop = 0; pop < 20; pop++) {
        int T = 0;
        double xtxxx = 0;
        point.resize(varx[pop].pn);
        edge.resize(varx[pop].sn);
        mass_initv3(point, 1, sp, varx[pop]);
        spring_initv3(edge, point, varx[pop].sck, varx[pop].sn, varx[pop].pn);
        ctm[0] = 0;
        ctm[1] = 0;
        ctm[2] = 0;
        for (int xtc = 0; xtc < varx[pop].pn; xtc++) {
            ctm[0] += point[xtc].p[0];
            ctm[1] += point[xtc].p[1];
            ctm[2] += point[xtc].p[2];
        }
        ctm[0] = ctm[0] / varx[pop].pn;
        ctm[1] = ctm[1] / varx[pop].pn;
        ctm[2] = ctm[2] / varx[pop].pn;
        while (T < 1000) {
            xtxxx = update_cube_firctionv3(edge, point, xtxxx, varx[pop].pn, varx[pop].sn);
            breathv3(edge, xtxxx, varx[pop].b, varx[pop].c, varx[pop].sn);
            T = T + 1;
        }
        double ctmf[3];
        ctmf[0] = 0;
        ctmf[1] = 0;
        ctmf[2] = 0;
        for (int xtc = 0; xtc < varx[pop].pn; xtc++) {
            ctmf[0] += point[xtc].p[0];
            ctmf[1] += point[xtc].p[1];
            ctmf[2] += point[xtc].p[2];
        }
        ctmf[0] = ctmf[0] / varx[pop].pn;
        ctmf[1] = ctmf[1] / varx[pop].pn;
        ctmf[2] = ctmf[2] / varx[pop].pn;
        double dist = 0;
        dist = sqrt(pow(ctmf[0] - ctm[0], 2) + pow(ctmf[1] - ctm[1], 2));
//        cout<<dist<<","<<xtxxx<<"\n";
//        cout<<ctmf[0]<<","<<ctmf[1]<<","<<ctmf[2]<<"\n";
        varx[pop].speed = dist / xtxxx;
    }
    int n = sizeof(varx) / sizeof(sets);
    sort(varx, varx + n, comparevar);
    bestspeed = varx[0].speed;
    for (int i = 0; i < 10000; i++) {
        copy_var(varx, tophalf);
        crossoverv3(tophalf, varx);
        for (int pop = 0; pop < 20; pop++) {
            int T = 0;
            double xtxxx = 0;
            point.resize(varx[pop].pn);
            edge.resize(varx[pop].sn);
            mass_initv3(point, 1, sp, varx[pop]);
            spring_initv3(edge, point, varx[pop].sck, varx[pop].sn, varx[pop].pn);
            ctm[0] = 0;
            ctm[1] = 0;
            ctm[2] = 0;
            for (int xtc = 0; xtc < varx[pop].pn; xtc++) {
                ctm[0] += point[xtc].p[0];
                ctm[1] += point[xtc].p[1];
                ctm[2] += point[xtc].p[2];
            }
            ctm[0] = ctm[0] / varx[pop].pn;
            ctm[1] = ctm[1] / varx[pop].pn;
            ctm[2] = ctm[2] / varx[pop].pn;
            while (T < 1000) {
                xtxxx = update_cube_firctionv3(edge, point, xtxxx, varx[pop].pn, varx[pop].sn);
                breathv3(edge, xtxxx, varx[pop].b, varx[pop].c, varx[pop].sn);
                T = T + 1;
            }
            double ctmf[3];
            ctmf[0] = 0;
            ctmf[1] = 0;
            ctmf[2] = 0;
            for (int xtc = 0; xtc < varx[pop].pn; xtc++) {
                ctmf[0] += point[xtc].p[0];
                ctmf[1] += point[xtc].p[1];
                ctmf[2] += point[xtc].p[2];
            }
            ctmf[0] = ctmf[0] / varx[pop].pn;
            ctmf[1] = ctmf[1] / varx[pop].pn;
            ctmf[2] = ctmf[2] / varx[pop].pn;
            double dist = 0;
            dist = sqrt(pow(ctmf[0] - ctm[0], 2) + pow(ctmf[1] - ctm[1], 2));
            varx[pop].speed = dist / xtxxx;
        }
        int n = sizeof(varx) / sizeof(sets);
        sort(varx, varx + n, comparevar);
        if (varx[0].speed > bestspeed) {
            bestspeed = varx[0].speed;
            cout << "found bestspeed: " << bestspeed << " at evo: " << i << "\n";
            ofstream fout1("GA_Robot_Variablev5_2.csv", ios::out);
            ofstream fout3("GA_Robot_positionv5_2.csv", ios::out);
            for (int fo = 0; fo < varx[0].sn; fo++) {
                fout1 << varx[0].b[fo] << "," << varx[0].c[fo] << "," << varx[0].sck[fo] << "\n";
            }
            for(int fox = 0; fox<varx[0].pn-1;fox++){
                fout3<< varx[0].p0[fox]<< "," <<varx[0].p1[fox]<< "," <<varx[0].p2[fox]<<"\n";
            }
            fout1.close();
            fout3.close();
        }
        fout2<<bestspeed<<"\n";
//        for(int ccx =0; ccx<20; ccx++){
//            fout4<<varx[ccx].speed<<",";
//            fout5<<varx[ccx].pn<<",";
//        }
//        fout4<<"\n";
//        fout5<<"\n";
        if(i%100==0){
            cout<<i<<"\n";
        }
    }
    fout2.close();
//    fout4.close();
//    fout5.close();
}