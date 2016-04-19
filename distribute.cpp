
#include <bits/stdc++.h>
#include <map>
#include <iostream>
#include <queue> 
#include <string>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#define GENE_NUM 22212
#define EXIST_GENE_NUM 9496
#define FEATURE_NUM 21234
using namespace std;
using namespace boost;
typedef long long ll;
typedef double val_type;

struct gene{
    int id;
    val_type val;
    gene(int _id=0, val_type _val=0):id(_id),
           val(_val){}
};

struct top_1d{
    int fid; // feature id
    int correct; // number of correctly classified genes
    top_1d(int _fid=0, int _correct=0):fid(_fid),
           correct(_correct){}
    bool operator<(const top_1d& rhs) const
    {
        return correct < rhs.correct;
    }
};

vector<vector<val_type > > unsorted_f; // unsorted_f[f][g]
vector<vector<gene > > sorted_f;
vector<vector<val_type > > back_f;
vector<dynamic_bitset<unsigned long> > bitmaps;
map<string, int> gene_id;  // gene name -> gid
vector<int > pos_gid;
vector<int > neg_gid;
vector<top_1d > top_f;


bool compareByvalue(const gene &a, const gene &b)
{
    return a.val < b.val;
}

bool compareByvalueDe(const gene &a, const gene &b)
{
    return a.val > b.val;
}

void load_matrix(FILE * fin){  // matrix fille
    clock_t t = clock();
    unsorted_f.resize(FEATURE_NUM);
    
    char name[200];
    for(int i=0;i<FEATURE_NUM+GENE_NUM;i++) {
        int x=fscanf(fin,"%s ", name);  // feature name 
        // can store feature name if want
    }

    for(int gid=0; gid<EXIST_GENE_NUM; gid++){
        int x = fscanf(fin,"%s ", name);
        string feature_name(name);
        gene_id[feature_name] = gid;  // gene id
        
        val_type entry;
        for(int j=0; j<FEATURE_NUM;j++){
            x=fscanf(fin,"%lf ", &entry);
            //fprintf(stderr,"=====%lf======\n",entry);
            unsorted_f[j].push_back(entry);
        }
    }
    fprintf(stderr,"=====total number of genes: %ld======\n",gene_id.size());
    t = clock() - t;
    fprintf (stderr, "load matrix: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

int Load_exp(FILE *fin_pos, FILE *fin_neg){
    clock_t t = clock();
    char name[200];
    while(fscanf(fin_pos,"%s ", name)!=EOF){
        string gene_name(name);
        map<string, int>::iterator it = gene_id.find(gene_name);
        if(it != gene_id.end()){
            pos_gid.push_back(it->second);
            //fprintf(stderr,"%d ",it->second);
        }
    }
    //fprintf(stderr,"\n");
    while(fscanf(fin_neg,"%s ", name)!=EOF){
        string gene_name(name);
        map<string, int>::iterator it = gene_id.find(gene_name);
        if(it != gene_id.end()){
            neg_gid.push_back(it->second);
            //fprintf(stderr,"%d ",it->second);
        }
    }
    //fprintf(stderr,"\n");
    fprintf(stderr,"=====number of positive genes:%ld; number of negtive genes. %ld======\n",pos_gid.size(),neg_gid.size());
    
    int total_gene_num = pos_gid.size()+neg_gid.size();
    sorted_f.resize(FEATURE_NUM);
    //bitmaps.resize(FEATURE_NUM);
    for(int f=0;f<FEATURE_NUM;f++){
        bitmaps.push_back(dynamic_bitset<unsigned long >(total_gene_num));   // create bitmaps  total_gene_num
        //fprintf(stderr, "%d %d total_gene_num====\n",f, total_gene_num);
        //cout<<bitmaps[f].size()<<endl;
        //bitmaps[f].set(0);
        for(int i=0; i<pos_gid.size();i++){
            int cur_gid = pos_gid[i];
            gene cur_gene(i, unsorted_f[f][cur_gid]); //(cur_gid, unsorted_f[f][cur_gid]);
            sorted_f[f].push_back(cur_gene);    
            // if(f==48)
            //     fprintf(stderr, "%lf ",unsorted_f[f][cur_gid]);
        }
        // if(f==48)
        //     fprintf(stderr, "\n");
        int pos_num =pos_gid.size();
        for(int i=0; i<neg_gid.size();i++){
            int cur_gid = neg_gid[i];
            gene cur_gene(i+pos_num, unsorted_f[f][cur_gid]);//(cur_gid, unsorted_f[f][cur_gid]);
            sorted_f[f].push_back(cur_gene);    
            // if(f==48)
            //     fprintf(stderr, "%lf ",unsorted_f[f][cur_gid]);
        }
        // if(f==48)
        //     fprintf(stderr, "\n");
    }

    fprintf(stderr,"------done with extracting submatrix for this experiment-------\n");
    t = clock() - t;
    fprintf (stderr, "load exp: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    return pos_gid.size();
}


/*transformation and get top features*/
void transformation(int pos_num){
    clock_t t = clock();
    priority_queue<top_1d> top_1d_features;
    back_f.resize(FEATURE_NUM);

    for(int f=0; f<FEATURE_NUM; f++){
        vector<gene > pos_genes(sorted_f[f].begin(), sorted_f[f].begin()+pos_num);
        vector<gene > neg_genes(sorted_f[f].begin()+pos_num, sorted_f[f].end());
        //median of positive 
        sort(pos_genes.begin(), pos_genes.end(), compareByvalue);
        val_type median_p =pos_genes[pos_genes.size()/2+1].val;  //x+
        sort(neg_genes.begin(), neg_genes.end(), compareByvalue);
        val_type median_n =neg_genes[neg_genes.size()/2+1].val; //x-
        
                
        val_type w = median_p - median_n;  //(x+ - x-)
        val_type intercept = -(median_p*median_p-median_n*median_n)/2;  //-(x+^2- x-^2)/2
        
        // if(f==48)
        //     fprintf(stderr, "%f %f %f %f\n",median_p, median_n, w, intercept);
        
        int num_correct=0;
        for(int i=0; i<pos_num;i++){
            val_type tmp=sorted_f[f][i].val * w + intercept;
            sorted_f[f][i].val = tmp; 
            back_f[f].push_back(tmp);
            // if(f==48)
            //     fprintf(stderr, "%lf ",sorted_f[f][i].val);
            if(sorted_f[f][i].val>0){
                bitmaps[f][i]=1;
                num_correct++;
            }
            
        }
        for(int i=pos_num; i<sorted_f[f].size();i++){
            val_type tmp= -(sorted_f[f][i].val * w + intercept);
            sorted_f[f][i].val = tmp;
            back_f[f].push_back(tmp);
            // if(f==48)
            //     fprintf(stderr, "%lf ",sorted_f[f][i].val);
            if(sorted_f[f][i].val>0){
                bitmaps[f][i]=1;
                num_correct++;
            }
        }
        top_1d cur_f(f,num_correct);
        top_1d_features.push(cur_f);
    }
    fprintf(stderr,"------done with transformation-------\n");
    t = clock() - t;
    fprintf (stderr, "transformation: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    
    // for(int i=0; i<FEATURE_NUM; i++){
    //     fprintf(stderr, "top-%d feature %d : with accuracy %d \n",i+1,top_1d_features.top().fid, top_1d_features.top().correct);
    //     top_1d_features.pop();
    // }
    for(int i=0; i<FEATURE_NUM; i++){
        int f = top_1d_features.top().fid;
        if(i<10)
        fprintf(stderr,"% ld set bit in top-%d feature %d",bitmaps[f].count(),i,f);
        top_f.push_back(top_1d_features.top());
        top_1d_features.pop();
    }
}


void distribution(){
    clock_t t = clock();
    FILE *fout = fopen("distribution.txt","a");
    int num_genes = back_f[0].size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    
    /*each feature pair*/
    int best = num_genes;
    int worst = 0;
    int best_f1,best_f2,worst_f1,worst_f2;
    int feature_considered =  FEATURE_NUM;
    long total_checked=0;

    for(int i=0; i<feature_considered; i++){
        //fprintf(stderr,"i=======%d\n",i);
        for(int j=i; j<feature_considered; j++){
            int f1 = top_f[i].fid;
            int f2 = top_f[j].fid;
            int missed = 0, g=0;
            for(; g<num_genes; g++){
                if(back_f[f1][g]+back_f[f2][g] <= 0)
                    missed++;
            }
            best = min(best,missed);

            worst = max(worst,missed);
            //total_checked+=g;
        }
    }
    fprintf(stderr,"===================\n");
    //fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.*****\n",best_f1,best_f2,best, feature_considered);
    t = clock() - t;
    fprintf (stderr, "==========best:%d===========\n",best);
    fprintf (stderr, "find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

    int total_range = best-worst;
    int rangeY[20];
    double dis = ((double)(best-worst))/20;
    memset(rangeY,0,sizeof(rangeY));

    for(int i=0; i<feature_considered; i++){
        for(int j=i; j<feature_considered; j++){
            int f1 = top_f[i].fid;
            int f2 = top_f[j].fid;
            int missed = 0, g=0;
            for(; g<num_genes; g++){
                if(back_f[f1][g]+back_f[f2][g] <= 0)
                    missed++;
            }
            int bucket = (int)((missed-worst)/dis);
            int correct_bucket = (bucket==20? 19: bucket);
            rangeY[correct_bucket]++; 

        }
    }
    double tmp = worst;
    fprintf(stderr,"============================D=I=S=T=R=I=B=U=T=I=O=N=============================================\n");
    fprintf(fout,"=========worst: %d======D=I=S=T=R=I=B=U=T=I=O=N======best: %d=========\n", worst, best);
    for(int i=0;i<20;i++){
        fprintf(stderr,"*******Missed range from [%lf to %lf) : %lld\n",tmp,tmp+dis,rangeY[i]);
        fprintf(fout,"%lf\t%lf\t%lld\n",tmp, tmp+dis, rangeY[i]);
        tmp+=dis;
    }
    t = clock() - t;
    fprintf (stderr, "find best feature pair again and write out: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

    fclose(fout);
}


void distribution_1d(){
    clock_t t = clock();
    FILE *fout = fopen("distribution_1d.txt","a");
    int num_genes = back_f[0].size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    
    /*each feature pair*/
    val_type best, worst;
    int feature_considered =  FEATURE_NUM;
    srand (time(NULL));

    for(int i=0; i<100; i++){
        int f = rand()%FEATURE_NUM ;
        sort(back_f[f].begin(),back_f[f].end());
        best = back_f[f][num_genes-1];
        worst = back_f[f][0];
        double dis = (best-worst)/20;
        int pre_gid=0,ptr=0;
        int rangeY[20];
        for(int g=0; g<num_genes; g++){
            if(back_f[f][g]>=worst + (ptr+1)*dis)
            {
               rangeY[ptr++] = g-pre_gid; 
               pre_gid=g;
            }
        }
        if(ptr!=20)
        {
            fprintf(stderr, "ptr:%d, f:%d; pos_gid:%d ======\n", ptr,f,pre_gid);
            rangeY[ptr]=num_genes-pre_gid;
        }
        val_type tmp = worst;
        fprintf(stderr,"=========================1D: D=I=S=T=R=I=B=U=T=I=O=N=============================================\n");
        fprintf(fout,"=========worst: %lf======D=I=S=T=R=I=B=U=T=I=O=N======best: %lf=========\n", worst, best);
        for(int i=0;i<20;i++){
            fprintf(stderr,"*******Missed range from [%lf to %lf) : %lld\n",tmp,tmp+dis,rangeY[i]);
            fprintf(fout,"%lf\t%lf\t%lld\n",tmp, tmp+dis, rangeY[i]);
            tmp+=dis;
        }
    }
    t = clock() - t;
    fprintf (stderr, "find 1d distribution write out: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    fclose(fout);
}


void distribution_equi_depth_1d(){
    clock_t t = clock();
    FILE *fout = fopen("distribution_depth_1d.txt","a");
    int num_genes = back_f[0].size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    
    /*each feature pair*/
    val_type best, worst;
    int feature_considered =  FEATURE_NUM;
    srand (time(NULL));

    for(int i=0; i<100; i++){
        int f = rand()%FEATURE_NUM ;
        sort(back_f[f].begin(),back_f[f].end());
        val_type bucket[21];
        for(int ptr=0; ptr<20; ptr++){
            bucket[ptr] = back_f[f][ptr*num_genes/20];
        }
        bucket[20] = back_f[f][num_genes-1];

        fprintf(stderr,"=========================1D: D=I=S=T=R=I=B=U=T=I=O=N=============================================\n");
        fprintf(fout,"===============D=I=S=T=R=I=B=U=T=I=O=N===============\n");
        for(int i=0;i<20;i++){
            fprintf(stderr,"*******Missed range from [%lf to %lf) : %d\n",bucket[i],bucket[i+1],num_genes/20);
            fprintf(fout,"%lf\t%lf\t%d\n",bucket[i],bucket[i+1],num_genes/20);
        }
    }
    t = clock() - t;
    fprintf (stderr, "find 1d distribution write out: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    fclose(fout);
}

int main(int argc, char* argv[]){
    clock_t t = clock();
    
    FILE * fin = fopen(argv[1],"r");
    FILE * fin_pos = fopen(argv[2],"r");
    FILE * fin_neg = fopen(argv[3],"r");

    //test();
    load_matrix(fin);
    int num_pos =Load_exp(fin_pos, fin_neg);
    transformation(num_pos);
    distribution_equi_depth_1d();
    //baseline_withbitmap();
    //sorted_list();
    //sorted_list_bi();
    //sorted_list_desc();
    t = clock() - t;
    fprintf (stderr, "It took me %d clicks (%f seconds) for the whole program.\n",t,((float)t)/CLOCKS_PER_SEC);
    return 0;
}