//
//  MOEA-D.cpp
//  MOEA/D
//
//  Created by Molin on 04/08/2017.
//  Copyright © 2017 merlin. All rights reserved.
//
//
#ifndef MOEA_D_MOEAD_H
#define MOEA_D_MOEAD_H
#include <queue>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <random>
#include <exception>
#include <omp.h>
#include "randG.h"
#include "loader.h"
#include <algorithm>
#include <float.h>
//Population's sizetrue
#define N 150
#define MM 100
//K nearest neighbors
#define K 11//Cannot assign 10. I don't know why.
bool mainprocess = false;
bool p_1_3_detail = false;
bool show_gene = false;
bool initial_detail = false;
bool genetic_state = false;
std::vector<asset>raw_asset;
void setting(bool *s, const std::vector<asset>assetList){
    mainprocess = s[0];
    p_1_3_detail = s[1];
    show_gene = s[2];
    initial_detail = s[3];
    for(auto asset_item:assetList){
        raw_asset.push_back(asset_item);
    }
}
struct solution{//i.e. x_i
    std::vector<int> gene;
    std::vector<asset> raw_data;
    std::vector<unsigned int> dominate;
    double objectives[2]={};
    //NSGA II Fitness
    int fitness;
    int solution_id;
    unsigned int domination_num_raw= 0;
    unsigned int domination_num = 0;
    double Crowding = 0.0;
    void init(){
        for(auto asset_item:raw_asset){
            raw_data.push_back(asset_item);
        }
    }
    solution(){
        init();
    }

    struct solution &operator = (const struct solution x){
        this->gene.clear();
        for(int i = 0; i<x.gene.size(); i++){
            this->gene.push_back(x.gene[i]);
        }
        this->objectives[0] = x.objectives[0];
        this->objectives[1] = x.objectives[1];
        this->fitness = x.fitness;
        this->domination_num_raw = x.domination_num_raw;
        this->domination_num = x.domination_num;
        this->solution_id = x.solution_id;
        this->Crowding = x.Crowding;
        for(int i = 0; i<x.dominate.size(); i++){
            unsigned int input_buffer = x.dominate[i];
            dominate.push_back(input_buffer);
        }
        for(int i = 0; i<x.raw_data.size(); i++){
            this->raw_data[i] = x.raw_data[i];
        }
        return *this;
    }

};
bool operator<(const solution &a, const solution &b){
    return a.domination_num<b.domination_num;
}
bool lessincome(const solution &x, const solution &y){
    return x.objectives[0]<y.objectives[0];
}
bool lessrank(const solution&x, const solution&y){
    return x.fitness<y.fitness;
}
bool lessdom(const solution&x, const solution&y){
    return x.domination_num<y.domination_num;
}
bool largedistance(const solution&x, const solution&y){
    return x.Crowding>y.Crowding;
}

struct population{
    std::vector<solution> xi;
    std::vector<std::vector<solution>>d_front;
    struct population &operator = (const struct population x){
        this->xi.assign(x.xi.begin(), x.xi.end());
        this->d_front.assign((x.d_front.begin()), x.d_front.end());
        return *this;
    }
};


struct lamb{
    //
    //Size of weight vector is determined by the number of objective of MOP.
    //To make it easy to implement, we only consider profit and risk, where k = 2.
    //According to the contribution of [2], the corresponding
    // H = 99, and N = 100;
    //
    double v[2] = {};
    int id = 0;
    std::vector<lamb>k_nearest;
    lamb &operator = (const lamb&y){
        for(int i = 0; i<2; i++){
            this->v[i] = y.v[i];
            i++;
        }
        this->v[0] = y.v[0];
        this->v[1] = y.v[1];
        for(int i = 0; i<this->k_nearest.size(); i++){
            this->k_nearest[i] = y.k_nearest[i];
        }
        this->id = y.id;
        return *this;
    }

};
bool operator<(lamb a, lamb b){
    return a.v[0] < b.v[0];
}

void util_outx(const population&x, int gen){
    for(auto item:x.xi){
        std::cerr<<item.objectives[0]<<"\t"<<item.objectives[1]<<"\t"<<gen<<std::endl;
    }
}
void util_evalue( solution &xi, std::vector<asset>assetList){
    double income = 0;
    double risk = 0;
    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i]!=0){
            income += xi.gene[i] * assetList[i].mean_income;
            risk += xi.gene[i] * assetList[i].diviation_r;
        }
    }
    xi.objectives[0] = income;
    xi.objectives[1] = risk;
}
void util_print_gene( const solution&xi){
    int counter = 0;
    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i] != 0) counter++;
        std::cerr<<xi.gene[i]<<"\t";
    }
    std::cerr<<"Number: "<<counter<<"\tFitness: "<<xi.objectives[0]<<" "<<xi.objectives[1]<<std::endl;
    /*
    for(int i = 0; i<xi.raw_data.size(); i++){
        std::cerr<<xi.raw_data[i].id<<"\t";
    }
    std::cerr<<std::endl;
     */
    for(int i = 0; i<xi.raw_data.size(); i++){
        std::cerr<<xi.raw_data[i].max_buy<<"\t";
    }
    std::cerr<<std::endl;
}
bool util_dominate(const struct solution &x, const struct solution&y){
    return (x.objectives[0]>y.objectives[0]&&x.objectives[1]<y.objectives[1]);
}
void util_repair_gene( solution &xi, const std::vector<asset>&assetList ){
    //1. Check Cardinality Constraint
    //2. Check Fund Constraint
    //3. Check min&max Buy Limit Constraint
    //
    //Check cardinality constraint
    //
    int M = 10;
    bool cac = false;
    int counter = 0;
    std::vector<int>port;

    //Caculate Fund Constraint
    int fund = 0;
    for(int i = 0; i<assetList.size(); i++){
        fund += assetList[i].current_price*assetList[i].holding;
    }

    //Caculate Solution's Fund
    int local_fund = 0;
    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i]!=0){
            counter++;
            port.push_back(i);
        }
    }
    solution local_best;
    if(counter>M){
        int differ = counter - M;
        for(int i = 0; i<differ; i++){
            solution local_best_2;
            for(int j = 0; j<port.size(); j++){
                solution xi_buffer;
                if(i == 0)
                    xi_buffer= xi;
                else
                    xi_buffer = local_best;
                if(xi_buffer.gene[port[j]]!=0)
                    xi_buffer.gene[port[j]] = 0;
                else
                    continue;
                //std::cerr<<"\nx"<<j<<" buffer:\n";
                util_evalue(xi_buffer, assetList);
                //util_print_gene(xi_buffer);
                if(j==0){
                    local_best_2 = xi_buffer;
                    //util_print_gene(local_best_2);
                }else{
                    if(util_dominate(xi_buffer, local_best_2)){
                        local_best_2 = xi_buffer;
                    }
                }
            }
            local_best = local_best_2;
            /*
            std::cerr<<"\nLocal best 2\n";
            util_print_gene(local_best_2);
            std::cerr<<"\nLocal best:\n";
            util_print_gene(local_best);
            std::cerr<<std::endl;
             */
        }
        xi = local_best;
    }

    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i]!=0){
            if(xi.gene[i]>assetList[i].max_buy){
                xi.gene[i] = assetList[i].max_buy;
            }
            if(xi.gene[i]<assetList[i].min_buy){
                xi.gene[i] = assetList[i].min_buy;
            }
        }
    }

    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i]!=0){
            local_fund+=xi.gene[i]*assetList[i].current_price;
        }
    }
    if(local_fund>fund){
        int differ = local_fund-fund;
        for(int i = 0; i<xi.gene.size(); i++){
            int minus = xi.gene[i]/local_fund*differ;
            if(xi.gene[i] - minus >=assetList[i].min_buy)
                xi.gene[i]-=minus;
        }
    }
    util_evalue(xi, assetList);
}
void util_findK( const int i, int k, int upper, int*h){
    int u = 0;
    int d[2] = {};
    bool up = false;
    if(k%2!=0)
        up = true;
    if(i+k/2+(up?1:0)>=upper){
        d[1] = upper-(i+k/2)+(up?1:0);
        d[0] = k-d[1];
    }
    else{
        d[1] = k/2+(up?1:0);
        if(i-k/2-(up?1:0)>0)
            d[0] = k/2+(up?1:0);
        else{
            d[0] = i;
            d[1] = k-i;
        }
    }
    for (int m = 1; m<=d[1]; m++) {
        int buffer = i+m;
        h[u] = buffer;
        u++;
    }
    for(int m = 1; m<=d[0]; m++){
        if(i == 0 )
            std::cerr<<"error"<<std::endl;
        int buffer = i-m;
        h[u] = buffer;
        u++;
    }
}
bool util_isfeasible( const solution &xi, const std::vector<asset>&assetList){
    int fund = 0, local_fund = 0, counter = 0;
    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i]!=0)
            counter++;
        fund+=assetList[i].current_price*assetList[i].holding;
        local_fund +=xi.gene[i]*assetList[i].current_price;
        if(xi.gene[i]>assetList[i].max_buy||xi.gene[i]<assetList[i].min_buy)
            return false;
    }
    if(counter>10)  return false;
    if(local_fund>fund) return false;
    return true;
}
void util_display_nn(const std::vector<lamb> &lamblist){
    for(int i = 0; i<lamblist.size(); i++) {
        std::cout << lamblist[i].v[0] << ", " << lamblist[i].v[1] << ":\t\t";
        for (int j = 0; j < lamblist[i].k_nearest.size(); j++) {
            std::cout << lamblist[i].k_nearest[j].v[0] << "," <<
                      lamblist[i].k_nearest[j].v[1] << "\t";
        }
        std::cout << std::endl;
    }
}
void util_display_nn_id( const std::vector<lamb> &lamblist){
    for(int i = 0; i<lamblist.size(); i++){
        std::cout<<lamblist[i].id<<":\t\t";
        for(int j = 0; j<lamblist[i].k_nearest.size(); j++){
            std::cout<<lamblist[i].k_nearest[j].id<<",\t";
        }
        std::cout<<std::endl;
    }
}
void util_fundDistribute( std::vector<struct asset>&p, double&fund, int n, int mode){

    if(p[n].max_buy<p[n].min_buy){
        check:
        bool buyable = false;
        for (int i = 0; i < p.size(); i++) {
            if (p[i].max_buy > p[i].min_buy&&fund/p[i].current_price>p[i].min_buy) {
                buyable = true;
            }
        }
        if(buyable){
            for (int i = 0; i < p.size(); i++) {
                if (p[i].max_buy > p[i].min_buy&&fund/p[i].current_price>p[i].min_buy) {
                    util_fundDistribute(p, fund, i, 0);
                }
            }
        }
        if(buyable) goto check;
        return;
    }

    if(fund/p[n].current_price<p[n].min_buy){
        for(int i = 0; i<p.size(); i++){
            if(fund/p[i].current_price>=p[i].min_buy){
                util_fundDistribute(p, fund, i, 0);
            }
        }
        return;
    }

    double dice = randG();
    double raw_buy = (dice * (p[n].max_buy-p[n].min_buy));
    if(raw_buy<0){
        std::cerr<<"!!!!Error\n";
    }
    int buy_number = static_cast<int>(dice * (p[n].max_buy-p[n].min_buy) + p[n].min_buy);

    double local_fund = buy_number*p[n].current_price;
    if(local_fund>fund){
        local_fund = fund;
        buy_number = local_fund/p[n].current_price;
        if(buy_number == 0){
            //std::cerr<<local_fund<<std::endl;
        }
        return;
    }
    p[n].fundpool+=local_fund;
    p[n].buy_asset_number += buy_number;
    p[n].max_buy -=buy_number;
    p[n].history.push_back(-buy_number);
    fund -=local_fund;
    if(n<=0||mode == 1){
        rn:
        double s = randG();
        util_fundDistribute(p, fund, int(s*p.size()),1);
    }
    else {
        util_fundDistribute(p, fund, n - 1, 0);
    }
}
void util_preprocess(std::vector<struct asset> &assetArray){
    for(int i = 0; i<assetArray.size(); i++){
        assetArray[i].max_buy = assetArray[i].max_buy+assetArray[i].holding;
        //assetArray[i].holding = 0;
    }
}
void util_ffunction(solution &xi,
                    const std::vector<asset> &assetArray,
                    double &income,
                    double &risk){
    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i]!=0){
            income+=xi.gene[i]*assetArray[i].mean_income;
            risk+=xi.gene[i]*assetArray[i].diviation_r;
        }
    }
    xi.objectives[0] = income;
    xi.objectives[1] = risk;
    if(p_1_3_detail) {
        std::cerr << "[1.3]: Income:\t" << income
                  << "\tRisk:\t" << risk << std::endl;
    }
}
double max(const double a, const double b){
    if(a>b)
        return a;
    return b;
}
double util_tchebycheff(const struct solution &x,
                        const struct lamb lb,
                        const double *z_star){
    double result;
    double temp[2];
    for(int i = 0; i<2; i++){
        temp[i] = std::abs((x.objectives[i]-z_star[i])*lb.v[i]);
    }
    result = max(temp[0], temp[1]);
    return result;
}

void process_updateZ( const population &x, double z[]){
    for(int i = 0; i<x.xi.size(); i++){
        if(x.xi[i].objectives[0]>z[0]){
            z[0] = x.xi[i].objectives[0];
        }
        if(x.xi[i].objectives[1]<z[1]){
            z[1] = x.xi[i].objectives[1];
        }
    }
    if(mainprocess){
        std::cerr<<"[2.3]:\tUpdate z:\t"
            <<"\n\t\tMax income:\t"<<z[0]<<"\t"<<"\tMin risk:\t"
            <<z[1]<<std::endl;
    }
}
void process_genetic( solution &m, solution &n, const std::vector<asset>&assetList, solution&trail, solution&trailx){
    //
    //Crossover - one point
    //
    int M = 10;
    solution trail_y;
    solution trail_x;
    int dot = randG()*m.gene.size();
    for(int i = 0; i<m.gene.size(); i++) {
        int trail_y_buffer = 0;
        int trail_x_buffer = 0;
        if (i < dot) {
            trail_y_buffer = m.gene[i];
            trail_x_buffer = n.gene[i];
            trail_x.raw_data[i] = n.raw_data[i];
            trail_y.raw_data[i] = m.raw_data[i];
        }
        else{
            trail_y_buffer = n.gene[i];
            trail_x_buffer = m.gene[i];
            trail_x.raw_data[i] = m.raw_data[i];
            trail_y.raw_data[i] = n.raw_data[i];
        }
        trail_y.gene.push_back(trail_y_buffer);
        trail_x.gene.push_back(trail_x_buffer);
    }
    if(genetic_state) {
        std::cerr << "M:\n";
        util_print_gene(m);
        std::cerr << "N:\n";
        util_print_gene(n);
        std::cerr << "Trail before repair\n";
        util_print_gene(trail_y);

        std::cerr << "Trail after repair\n";
        util_print_gene(trail_y);
    }
    //util_repair_gene(trail_y, assetList);
    //util_repair_gene(trail_x, assetList);

    double roll_dice = randG();
    double mu = 1;
    if(roll_dice<mu){
        //
        //Mutation
        //
        int countery = 0;
        int counterx = 0;
        std::vector<asset>trail_asset;
        std::vector<asset>trail_assetx;
        std::vector<int>availible_listx;
        std::vector<int>availible_listy;
        double trail_fund = 0;
        double trail_fundx= 0;
        double fund = 0;
        for(int i = 0; i<trail_y.gene.size(); i++){
            fund += assetList[i].current_price*assetList[i].holding;
            if(trail_y.gene[i]!=0){
                countery++;
                trail_fund += trail_y.gene[i] * assetList[i].current_price;
            }
            if(trail_y.gene[i] == 0){
                availible_listy.push_back(i);
            }
            if(trail_x.gene[i] == 0){
                availible_listx.push_back(i);
            }
            if(trail_x.gene[i]!=0){
                counterx++;
                trail_fundx +=trail_x.gene[i]*assetList[i].current_price;
            }
        }
        std::vector<int>index_random;
        std::vector<int>index_randomx;
#pragma omp parallel num_threads(8)
#pragma omp section
        {
            for(int i = 0; i<M-counterx;i++) {
                int temp_index;
                roll1:
                temp_index = randG() * availible_listx.size();
                for (auto item:index_randomx) {
                    if (temp_index == item)
                        goto roll1;
                }
                index_randomx.push_back(temp_index);
                asset temp_trail_meta = trail_x.raw_data[availible_listx[temp_index]];
                //std::cerr<<"Empty index:\t"<<availible_listy[temp_index]<<"\t";
                trail_assetx.push_back(temp_trail_meta);
            }
        }
#pragma omp section
        {
            for(int i = 0; i<M-countery; i++){
                int temp_index;
                roll:
                temp_index = randG()*availible_listy.size();
                for(auto item:index_random){
                    if(temp_index == item)
                        goto roll;
                }
                index_random.push_back(temp_index);
                asset temp_trail_meta = trail_y.raw_data[availible_listy[temp_index]];
                //std::cerr<<"Empty index:\t"<<availible_listy[temp_index]<<"\t";
                trail_asset.push_back(temp_trail_meta);
            }
        };

#pragma omp barrier
#pragma omp section
        {
            if(countery<M&&trail_fund<fund){
                int differ = countery-M;
                double buffer_fund = fund - trail_fund;
                //std::cerr<<"Size of trail asset\t"<<trail_asset.size()<<std::endl;
                util_fundDistribute(trail_asset, buffer_fund, trail_asset.size()-1, 0);
            }
        };
#pragma omp section
        {
            if(counterx<M&&trail_fundx<fund){
                int differ = counterx-M;
                double buffer_fund = fund - trail_fundx;
                //std::cerr<<"Size of trail asset\t"<<trail_asset.size()<<std::endl;
                util_fundDistribute(trail_assetx, buffer_fund, trail_assetx.size()-1, 0);
            }
        };

        /*
        for(int i = 0; i<trail_asset.size(); i++){
            std::cerr<<trail_asset[i].buy_asset_number<<"\t";
        }
         */
        for(int i = 0; i<trail_y.gene.size(); i++){
            for(int j = 0; j<trail_asset.size(); j++){
                if(i == trail_asset[j].id){
                    if(trail_y.gene[i] != 0)
                        std::cerr<<"Error\n";
                    trail_y.gene[i] = trail_asset[j].buy_asset_number;
                    trail_y.raw_data[i] = trail_asset[j];
                }
            }
        }
#pragma omp section
        {
            for(int i = 0; i<trail_x.gene.size(); i++){
                for(int j = 0; j<trail_assetx.size(); j++){
                    if(i == trail_assetx[j].id){
                        if(trail_x.gene[i] != 0)
                            std::cerr<<"Error\n";
                        trail_x.gene[i] = trail_assetx[j].buy_asset_number;
                        trail_x.raw_data[i] = trail_assetx[j];
                    }
                }
            }
        };
    }
    util_evalue(trail_y, assetList);
    util_evalue(trail_x, assetList);
    //util_print_gene(trail_y);
    trail = trail_y;
    trailx = trail_x;
}
void process_updateP(const std::vector<lamb>&lamblist,
             double z[2],
             const population &x,
             population &new_x,
             const std::vector<asset>&assetList ){
    if(mainprocess){
        std::cerr<<"[2]:\tUpdate\n";
    }
    //
    //  [2.1] Reproduction
    //
    if(mainprocess){
        std::cerr<<"[2.1]:\tStart\tReporduction\n";
    }
    omp_set_num_threads(32);
    int chunksize = 10;
#pragma omp parallel shared(chunksize) private(i, t_id)
    {
#pragma omp for schedule(dynamic, chunksize)
        for (int i = 0; i < N; i++) {
            solution trail_buffer;
            int m = lamblist[i].k_nearest[static_cast<int>(randG() * lamblist[i].k_nearest.size())].id;
            int n = lamblist[i].k_nearest[static_cast<int>(randG() * lamblist[i].k_nearest.size())].id;

            solution m_buffer = x.xi[m];
            solution n_buffer = x.xi[n];
            //if(mainprocess) std::cerr<<"[2.1]:\tm:"<<m<<"\tn:"<<n<<std::endl;
            //process_genetic(m_buffer, n_buffer, assetList, trail_buffer);
            new_x.xi.push_back(trail_buffer);
        }
    }
}
void process_updateN( population&new_x,
                      population&x,
                      const std::vector<asset>&assetList,
                      const std::vector<lamb>&lambList,
                      const double z[]){
    if(mainprocess){
        std::cerr<<"[2.4]: Start\tUpdate Neighborhood Population"<<std::endl;
    }
    for(int i = 0; i<N; i++){
        for(int j = 0; j<lambList[i].k_nearest.size(); j++){
            double y = util_tchebycheff(new_x.xi[i], lambList[i].k_nearest[j], z);
            double x_j = util_tchebycheff(x.xi[lambList[i].k_nearest[j].id], lambList[i].k_nearest[j], z);
            if (y<=x_j){
                x.xi[lambList[i].k_nearest[j].id] = new_x.xi[i];
            }
        }
    }
    if(mainprocess)   std::cerr<<"[1.3]: *End*\tUpdate Neighborhood Population"<<std::endl;
}
void process_updateP( population&pool,
                      const std::vector<asset>&assetList){
    clock_t start, end;
    start = clock();
    if(mainprocess){
        std::cerr<<"[1.3]: START\tUpdate Pool"<<std::endl;
    }
    std::vector<solution> sub_population;
    for(int i = 0; i<pool.xi.size(); i++){
        int index1 = randG()*pool.xi.size();
        int index2 = randG()*pool.xi.size();
        solution parent1 = pool.xi[index1];
        solution parent2 = pool.xi[index2];
        solution trailx, traily;
        process_genetic(parent1, parent2, assetList, traily, trailx);
        sub_population.push_back(traily);
        sub_population.push_back(trailx);
    }
    pool.xi.insert(pool.xi.end(), sub_population.begin(), sub_population.end());
    end = clock();
    if(mainprocess){
        size_t runtime = end-start;
        std::cerr<<"[1.3]: Ranks:\t"<<pool.xi.size()<<std::endl;
        std::cerr<<"[1.3]: *END\tRuntime:\t"<<runtime<<std::endl;
    }
}
void process_updateEP( const population&x, population&ep){
    if(mainprocess){
        std::cerr<<"[2.5]: Start\tUpdate EP\n";
    }
    for(int i = 0; i<x.xi.size(); i++){
        for(int j = 0; j<x.xi.size(); j++){
            if(util_dominate(x.xi[i], ep.xi[j])){
                ep.xi[j] = x.xi[i];
            }
        }
    }
    if(mainprocess){
        std::cerr<<"[2.5]: *End*\tUpdate EP\n";
    }
}
void process_selection( population &x){
    population pool;
    clock_t start, end;
    start = clock();
    if(mainprocess){
        std::cerr<<"[1.3]: START\tSelection"<<std::endl;
    }

    int pool_size =125;
    int tour_size = 2;
    while(pool.xi.size()<pool_size){
        double dice1 = randG();
        double dice2 = randG();
        while(dice2 == dice1)
            dice2 = randG();

        int index1 = x.xi.size()*dice1;
        int index2 = x.xi.size()*dice2;

        if(x.xi[index1].fitness<x.xi[index2].fitness){
            solution buffer1 = x.xi[index1];
            pool.xi.push_back(buffer1);
        }
        else
        {
            if(x.xi[index1].fitness==x.xi[index2].fitness){
                solution buffer2= (x.xi[index1].Crowding>x.xi[index2].Crowding? x.xi[index1]:x.xi[index2]);
                pool.xi.push_back(buffer2);
            }
            else
            {
                solution buffer3 = x.xi[index2];
                pool.xi.push_back(buffer3);
            }
        }
    }
    x = pool;
    end = clock();
    if(mainprocess){
        size_t runtime = end-start;
        std::cerr<<"[1.3]: Pool rank:\t"<<pool.xi.size()<<std::endl;
        std::cerr<<"[1.3]: *END\tRuntime:\t"<<runtime<<std::endl;
    }
}

void process_Fast_Nondominated( population&x){
    clock_t start, end;
    start = clock();
    if(mainprocess){
        std::cerr<<"[1.x]: START\tFast Nondominated Sort"<<std::endl;
    }

    unsigned int n_size = x.xi.size();
    std::priority_queue<solution> handle_list;
    std::vector<solution>buffer_front;
#pragma omp parallel num_threads(16)
    {
#pragma omp for schedual(runtime)
        for(int i = 0; i<n_size; i++){
            for(int j = 0; j<n_size; j++){
                if(util_dominate(x.xi[j], x.xi[i])){
#pragma omp critical
                    {
                        x.xi[i].domination_num_raw++;
                        x.xi[i].domination_num++;
                    };
                }
            }
            /*
            if(x.xi[i].domination_num_raw == 0){
                solution buffer = x.xi[i];
                x.xi[i].fitness = 1;
                buffer_front.push_back(buffer);
            }
             */
        }
    }
    //x.d_front.push_back( buffer_front);
    //Rank
    sort(x.xi.begin(), x.xi.end(), lessdom);
    //Reset solution id;
    for(int i = 0; i<x.xi.size(); i++){
        x.xi[i].solution_id = i;
    }
    int rank = 1;
    int temper_dn = x.xi[0].domination_num;
    std::vector<solution> front_local;
    std::vector<solution>emp;
    x.d_front.push_back(front_local);
    for(int i = 0; i<n_size; i++){
        if(x.xi[i].domination_num == temper_dn){
            x.xi[i].fitness = rank;
            solution buffer = x.xi[i];
            if(buffer.solution_id!=x.xi[i].solution_id){
                std::cerr<<"RAW ERROR\n";
            }
            x.d_front[rank-1].push_back(buffer);
        }
        else{
            rank++;
            std::vector<solution> buffer;
            x.d_front.push_back(buffer);
            //buffer.assign(front_local.begin(), front_local.end());
            temper_dn = x.xi[i].domination_num;
            front_local.clear();
            x.xi[i].fitness = rank;
            solution buffer2 = x.xi[i];
            x.d_front[rank-1].push_back(buffer2);
        }
    }
    for(int i = 0; i<x.d_front.size(); i++){
        for(int j = 0; j<x.d_front[i].size(); j++){
            if(x.d_front[i][j].solution_id>x.xi.size()){
                std::cerr<<"i, j"<<i<<j<<std::endl;
            }
        }
    }

    if(0){

        for(int i = 0; i<x.d_front.size(); i++){
            std::cerr<<"[1.4]: Local rank of "<<i+1<<"\t"<<x.d_front[i].size()<<std::endl;
        }
    }

    //
    //[1.2.2] Calculate Crowding Distance
#pragma omp for schedual(dynamic)
    for(int i = 0; i<x.d_front.size(); i++){
        //Less income first
        sort(x.d_front[i].begin(), x.d_front[i].end(), lessincome);
        sort(x.d_front[i].begin(), x.d_front[i].end(), lessincome);
        if(x.d_front[i].size()>2){
            x.d_front[i][0].Crowding = x.d_front[i][x.d_front[i].size()-1].Crowding = DBL_MAX;
            for(int j = 1; j<x.d_front[i].size()-1; j++){
#pragma omp critical
                {
                    double neibor = (x.d_front[i][j+1].objectives[0]- x.d_front[i][j-1].objectives[0]);
                    double max_diff = (x.d_front[i][x.d_front[i].size()-1].objectives[0]-x.d_front[i][0].objectives[0]);
                    x.d_front[i][j].Crowding += neibor/max_diff;
                };
                //std::cerr<<"[1.2]: "<<i<<" "<<j<<" "<<"neibor\t"<<neibor<<"\tdiffer\t"<<max_diff<<"\t"<<x.d_front[i][j].Crowding<<std::endl;
            }
        }
    }

    for(int i = 0; i<x.d_front.size(); i++){
        for(int j = 0; j<x.d_front[i].size(); j++){
            double buffer = x.d_front[i][j].Crowding;
            int local_id = x.d_front[i][j].solution_id;
            x.xi[x.d_front[i][j].solution_id].Crowding = x.d_front[i][j].Crowding;
            x.xi[x.d_front[i][j].solution_id].fitness = x.d_front[i][j].fitness;
        }
    }

    end = clock();
    if(mainprocess){
        size_t runtime = end-start;
        std::cerr<<"[1.x]: Ranks:\t"<<x.d_front.size()<<std::endl;
        std::cerr<<"[1.x]: *END\tRuntime:\t"<<runtime<<std::endl;
    }
}

void process_fast_nondominated_lite( population&x){
    clock_t start, end;
    start = clock();
    if(mainprocess){
        std::cerr<<"[1.4]: START\tShrink"<<std::endl;
    }
    if(x.xi.size()>MM){
        process_Fast_Nondominated(x);
        sort(x.xi.begin(), x.xi.end(), lessrank);

        int local_rank =100;
        if(x.xi[MM-1].fitness == x.xi[MM].fitness){
            local_rank = x.xi[MM-1].fitness;
            sort(x.d_front[local_rank].begin(), x.d_front[local_rank].end(), largedistance);
            x.xi.resize(MM);
            while(x.xi.back().fitness == local_rank){
                x.xi.pop_back();
            }
            int iter = 0;
            while(x.xi.size()<MM){
                solution buffer = x.d_front[local_rank][iter];
                iter++;
                x.xi.push_back(buffer);
                if(iter>=x.d_front[local_rank].size())
                    break;
            }
        }
        else{
            x.xi.resize(MM);
        }

    }
    end = clock();
    if(mainprocess){
        size_t runtime = end-start;
        std::cerr<<"[1.4]: *END\tRuntime:\t"<<runtime<<std::endl;
    }
}
bool init_lamb( const int &k, const int &H, std::vector<lamb>&lamb_list ){
    for(int i = 0; i<N; i++){
        double source_number = randG();
        struct lamb buffer;
        buffer.id = i;
        buffer.v[0] = static_cast<double>(static_cast<int>(source_number*H))/ static_cast<double>(H);
        //std::cerr<<i<<":"<<buffer.v[0]<<std::endl;
        buffer.v[1] = 1 - buffer.v[0];
        lamb_list.push_back(buffer);
    }
    return true;
}
void init_distance( std::vector<lamb> &lamb_list ){
    if(mainprocess)   std::cerr<<"[1.2]: Start\tInit Distance"<<std::endl;
    std::priority_queue<lamb>nearest_list;
    int length = lamb_list.size();
    for(int i = 0; i<length; i++){
        lamb buffer = lamb_list.back();
        //std::cout<<buffer.v[0]<<"\t"<<buffer.v[1]<<std::endl;
        lamb_list.pop_back();
        nearest_list.push(buffer);
    }
    std::vector<lamb> toHandle;
    while(!nearest_list.empty()){
        lamb buffer = nearest_list.top();
        toHandle.push_back(buffer);
        nearest_list.pop();
    }

    for(int i = 0; i<toHandle.size(); i++){
        lamb buffer = toHandle[i];
        int h[K] = {};
        //util_findK(i, K, toHandle.size(), h);
        for(int j = 0; j<K; j++){
            buffer.k_nearest.push_back(toHandle[h[j]]);
        }
        lamb_list.push_back(buffer);
    }
    if(mainprocess)   std::cerr<<"[1.2]: *End*\tInit Distance"<<std::endl;
}


void init_solution( std::vector<struct asset> &tobuy,
                    double&fund){
    util_fundDistribute(tobuy, fund, tobuy.size()-1, 0);
    if(initial_detail){
        for(int i = 0; i<tobuy.size(); i++){
            std::cerr<<"ID: "<<tobuy[i].id<<"\tMax: "<<tobuy[i].max_buy<<"\tMin: "
                     <<tobuy[i].min_buy<<"\tBuy: "<<tobuy[i].buy_asset_number<<"\t";
            for(auto item:tobuy[i].history){
                std::cerr<<item<<"\t";
            }
            std::cerr<<std::endl;
        }
    }
}
void init_population(struct population &x,
                     const std::vector <struct asset> &asset,
                     const struct Constraint constraint,
                     const double (&correlation)[31][31]){
    clock_t start, end;
    start = clock();
    if(mainprocess)   std::cerr<<"[1.1]: Start\tInit Population"<<std::endl;
    //Cardinality Constraint M:
    int M = constraint.max_assets;
    int num_assets = constraint.num_assets;
    //std::cerr<<"Number of Assets:\t"<<num_assets<<std::endl;
    //
    double fundpool = 0;
    for(int i = 0; i<asset.size(); i++) {
        fundpool += asset[i].current_price * asset[i].holding;
    }
    for(int i = 0; i<N; i++){
        solution xi_buffer;
        struct Set candidate;
        //Handling cardinality constraint by using pointers.
        //std::cerr<<"[1.3]: \t\tStart\tassign assets\n";
        std::vector<struct asset>tobuy;
        for(int i = 0; i<M; i++){
            int buffer = static_cast<int>(randG()*(num_assets+1));
            if(buffer>=num_assets)
                buffer = -1;//Hold cash
            //std::cerr<<buffer<<"\t";
            if(buffer>=0) {
                if(!candidate.isin(buffer)){
                    candidate.data.push_back(buffer);
                }
            }
        }
        for(int i = 0; i<candidate.data.size(); i++){
            int buffer_id = candidate.data[i];
            struct asset asset_buffer = asset[buffer_id];
            asset_buffer.id = buffer_id;
            tobuy.push_back(asset_buffer);
        }
        std::priority_queue<int, std::vector<int>, std::greater<int>>tobuy_ase;
        for(int i = 0; i<tobuy.size(); i++){
            //std::cerr<<tobuy[i].id<<"\t";
            int buffer = tobuy[i].id;
            tobuy_ase.push(buffer);
        }
        if(p_1_3_detail) {
            std::cerr << "[1.3]:\tPort " << i << ":\t";
            while (!tobuy_ase.empty()) {
                std::cerr << tobuy_ase.top() << "\t";
                tobuy_ase.pop();
            }
            std::cerr << std::endl;
        }
        //if(mainprocess)   std::cerr<<"[1.3]: \t\tEnd\t\tassign assets\n";
        //if(mainprocess)   std::cerr<<"[1.3]: \t\tStart\tInitialize Solution\n";
        double fund_buffer = fundpool;
        init_solution(tobuy, fund_buffer);
        if(p_1_3_detail)   std::cerr<<"[1.3]: Fund remain:"<<fund_buffer<<std::endl;
        //if(mainprocess)   std::cerr<<"[1.3]: \t\tEnd\t\tInitialize Solution\n";
        for(int j = 0; j<asset.size(); j++){
            xi_buffer.gene.push_back(0);
        }
        for(auto item:tobuy){
            xi_buffer.gene[item.id] = item.buy_asset_number;
            xi_buffer.raw_data[item.id] = item;
        }
        double income_buffer = 0;
        double risk_buffer;
        util_ffunction(xi_buffer, asset, income_buffer, risk_buffer);
        xi_buffer.solution_id = i;
        x.xi.push_back(xi_buffer);
        if(show_gene) {
            for (int i = 0; i < xi_buffer.gene.size(); i++) {
                std::cout << xi_buffer.gene[i] << " ";
            }
            std::cout << std::endl;
        }
    }
    end = clock();
    if(mainprocess){
        size_t runtime = end-start;
        std::cerr<<"[1.1]: *END\tRuntime:\t"<<runtime<<std::endl;
    }
}
size_t covariance(const std::vector<int>&x,
                  const std::vector<struct asset> &asset,
                  const double (&correlation)[31][31]){
    size_t cov = 0;
    for(int i = 0; i<31; i++){
        for(int j = 0; j<31; j++){
            cov +=asset[i].current_price*(x[i]+asset[i].holding)*asset[j].current_price*(x[j]+asset[j].holding)*correlation[i][j];
        }
    }
    return cov;
}

#endif //MOEA_D_MOEAD_H
