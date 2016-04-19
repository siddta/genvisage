

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
    int fcorrect;
    gene(int _id=0, val_type _val=0, int _fcorrect=0):id(_id),
           val(_val),fcorrect(_fcorrect){}
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

struct gene_fnum{
    int gid; // feature id
    long fnum; // number of correctly classified genes
    gene_fnum(int _gid=0, long _fnum=0):gid(_gid),
           fnum(_fnum){}
    bool operator<(const gene_fnum& rhs) const
    {
        return fnum < rhs.fnum;
    }
};

struct feature{
    vector<val_type > vals;
    int correct;
    int fid;
    dynamic_bitset<unsigned long> bitmap;
    
    feature(vector<val_type > _vals, int _correct, int _fid,dynamic_bitset<unsigned long> _bitmap):vals(_vals),correct(_correct),fid(_fid),bitmap(_bitmap){}
    feature(){}
    bool operator<(const feature& a) const
    {
        return correct > a.correct;
    }
};

struct top_2d{
    int f1;
    int f2;
    int miss;
    top_2d(int _f1=0, int _f2=0, int _miss=0):f1(_f1),f2(_f2),miss(_miss){}
    bool operator<(const top_2d& a) const
    {
        return miss < a.miss;
    }
};


vector<vector<val_type > > unsorted_f; // unsorted_f[f][g]
vector<feature > sorted_f;
vector<feature > back_f;
vector<vector<val_type > > hist;  //store the histogram for each feature 
vector<vector<val_type > > hist2; 
//vector<dynamic_bitset<unsigned long> > bitmaps;
map<string, int> gene_id;  // gene name -> gid
vector<int > pos_gid;
vector<int > neg_gid;
vector<top_1d > top_f;


bool compareByCorrectD(const feature &a, const feature &b)
{
    return a.correct > b.correct;
}


bool compareByvalue(const gene &a, const gene &b)
{
    return a.val < b.val;
}

bool compareByfeature(const gene &a, const gene &b)
{
    return a.fcorrect < b.fcorrect;
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
        string gene_name(name);
        gene_id[gene_name] = gid;  // gene id
        
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
    back_f.resize(FEATURE_NUM);

    for(int f=0;f<FEATURE_NUM;f++){
        //bitmaps.push_back(dynamic_bitset<unsigned long >(total_gene_num));   // create bitmaps  total_gene_num
        back_f[f].bitmap = dynamic_bitset<unsigned long >(total_gene_num);
        //fprintf(stderr, "%d %d total_gene_num====\n",f, total_gene_num);
        //cout<<bitmaps[f].size()<<endl;
        //bitmaps[f].set(0);
        for(int i=0; i<pos_gid.size();i++){
            int cur_gid = pos_gid[i];
            back_f[f].vals.push_back(unsorted_f[f][cur_gid]);
        }
        int pos_num =pos_gid.size();
        for(int i=0; i<neg_gid.size();i++){
            int cur_gid = neg_gid[i];
            back_f[f].vals.push_back(unsorted_f[f][cur_gid]);    
        }
    }

    fprintf(stderr,"------done with extracting submatrix for this experiment-------\n");
    t = clock() - t;
    fprintf (stderr, "load exp: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    return pos_gid.size();
}


/*transformation and get top features*/
void transformation(int pos_num){
    clock_t t = clock();
    //priority_queue<top_1d> top_1d_features;
    int num_genes = back_f[0].vals.size();
    
    for(int f=0; f<FEATURE_NUM; f++){
        vector<val_type > pos_genes(back_f[f].vals.begin(), back_f[f].vals.begin()+pos_num);
        vector<val_type > neg_genes(back_f[f].vals.begin()+pos_num, back_f[f].vals.end());
        //median of positive 
        sort(pos_genes.begin(), pos_genes.end());
        val_type median_p =pos_genes[pos_genes.size()/2+1];  //x+
        sort(neg_genes.begin(), neg_genes.end());
        val_type median_n =neg_genes[neg_genes.size()/2+1]; //x-
        
                
        val_type w = median_p - median_n;  //(x+ - x-)
        val_type intercept = -(median_p*median_p-median_n*median_n)/2;  //-(x+^2- x-^2)/2
        
        int num_correct=0;
        for(int i=0; i<pos_num;i++){
            val_type tmp=back_f[f].vals[i] * w + intercept; //CHANGE: int TO val_type
            back_f[f].vals[i] = tmp; 
           
            if(back_f[f].vals[i]>0){
                //bitmaps[f][i]=1;
                back_f[f].bitmap[i]=1;
                num_correct++;
            }
        }
        for(int i=pos_num; i<num_genes; i++){
            val_type tmp= -(back_f[f].vals[i] * w + intercept); ////CHANGE: int TO val_type
            back_f[f].vals[i] = tmp;
            if(back_f[f].vals[i]>0){
                //bitmaps[f][i]=1;
                back_f[f].bitmap[i]=1;
                num_correct++;
            }
        }
        back_f[f].correct = num_correct;
        back_f[f].fid = f;
        //top_1d cur_f(f,num_correct);  // can be omitted 
        //top_1d_features.push(cur_f);  // can be omitted
    }
    fprintf(stderr,"------done with transformation-------\n");
    t = clock() - t;
    fprintf (stderr, "transformation: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    
    /*can be omitted*/
    /*
    for(int i=0; i<FEATURE_NUM; i++){
        int f = top_1d_features.top().fid;
        if(i<10)
            fprintf(stderr,"%ld set bit in top-%d feature %d \n",back_f[f].bitmap.count(),i,f);
        top_f.push_back(top_1d_features.top());
        top_1d_features.pop();
    }*/
    
}

void calculate_50percent(){
    int total_bad_f=0;
    int num_genes = back_f[0].vals.size();
    for(int f=0; f<FEATURE_NUM;f++){
        int miss = 0;
        for(int g=0;g<num_genes;g++){
            if(back_f[f].bitmap[g]==0)
                miss++;
        }
        if(miss>=num_genes/2)
            total_bad_f++;
    }
    fprintf(stderr,"bad features:%d; good feature: %d \n",total_bad_f, FEATURE_NUM-total_bad_f);
}

void sort_genes_fpair(){
    clock_t t = clock();
    fprintf(stderr,"here-============1\n");
    int num_genes = back_f[0].vals.size();
    vector<gene_fnum > gene_correct;
    fprintf(stderr,"here-============2\n");
    for(int g=0; g<num_genes;g++){
        vector<val_type > cur_gene;  //genes[0-genenumber]
        for(int f=0; f<FEATURE_NUM;f++){
            cur_gene.push_back(back_f[f].vals[g]);
        }
        sort(cur_gene.begin(),cur_gene.end());  //sort features
        int ptr1=0, ptr2=FEATURE_NUM-1;
        long total_fpair=0;
        while(ptr1<ptr2){
            if(cur_gene[ptr1]+cur_gene[ptr2]<=0){
                ptr1++;
                //total_fpair+=(FEATURE_NUM-ptr2-1);
            }
            else{
                total_fpair+=(ptr2-ptr1);
                ptr2--;
            }
        }
        gene_fnum tmp(g,total_fpair);
        gene_correct.push_back(tmp);
    }
    fprintf(stderr,"here-============3\n");
    sort(gene_correct.begin(),gene_correct.end());
    for(int f=0; f<FEATURE_NUM;f++){
        sorted_f[f].correct = back_f[f].correct;
        sorted_f[f].fid = back_f[f].fid;
        sorted_f[f].bitmap = dynamic_bitset<unsigned long >(num_genes);
        for(int g=0; g<num_genes;g++){
            sorted_f[f].vals.push_back(back_f[f].vals[gene_correct[g].gid]);  //gene_correct[g].fid: gid of gth gene
            sorted_f[f].bitmap[g] = back_f[f].bitmap[gene_correct[g].gid];  // update bitmap
        }
    }
    t = clock() - t;
    fprintf (stderr, "sort_genes FPAIR: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

void sort_genes(){
    clock_t t = clock();
    sorted_f.resize(FEATURE_NUM);
    vector<top_1d > gene_correct;  // use top_1d (id feature)-- fid->gid
    int num_genes = back_f[0].vals.size();
    for(int g=0;g<num_genes;g++){
        int g_corrects = 0;
        for(int f=0;f<FEATURE_NUM;f++){
            if(back_f[f].bitmap[g]==1)
            //if(bitmaps[f][g]==1)
                g_corrects++;
        }
        top_1d tmp(g,g_corrects);
        gene_correct.push_back(tmp);
    }
    sort(gene_correct.begin(),gene_correct.end());
    
    for(int f=0; f<FEATURE_NUM;f++){
        sorted_f[f].correct = back_f[f].correct;
        sorted_f[f].fid = back_f[f].fid;
        sorted_f[f].bitmap = dynamic_bitset<unsigned long >(num_genes);
        for(int g=0; g<num_genes;g++){
            sorted_f[f].vals.push_back(back_f[f].vals[gene_correct[g].fid]);  //gene_correct[g].fid: gid of gth gene
            sorted_f[f].bitmap[g] = back_f[f].bitmap[gene_correct[g].fid];  // update bitmap
        }
    }
    t = clock() - t;
    fprintf (stderr, "sort_genes : %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}



void baseline_unsorted_noprune(){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    
    /*each feature pair*/
    int best = num_genes;
    int best_f1,best_f2;
    int feature_considered =  FEATURE_NUM;
    long total_checked=0;
    for(int f1=0; f1<feature_considered; f1++){
        for(int f2=f1+1; f2<feature_considered; f2++){

            int missed = 0,g=0;
            for(; g<num_genes; g++){
                if(back_f[f1].vals[g]+back_f[f2].vals[g] <= 0)
                    missed++;
            }
            if(missed < best){
                best = missed;
                best_f1 = back_f[f1].fid;
                best_f2 = back_f[f2].fid;
            }
            total_checked+=g;
        }
    }
    
    fprintf(stderr,"=========%ld=======\n",total_checked);
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.*****\n",best_f1,best_f2,best, feature_considered);
    t = clock() - t;
    fprintf (stderr, "baseline_unsorted_noprune find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}


void baseline_unsorted(){  //unsorted features
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    
    FILE* fout = fopen("geneNum_unsorted.txt","w");
    /*each feature pair*/
    int best = num_genes;
    int best_f1,best_f2;
    int feature_considered =  FEATURE_NUM;
    long total_checked=0;
    for(int f1=0; f1<feature_considered; f1++){
        int g_f1=0;
        for(int f2=f1+1; f2<feature_considered; f2++){
            int missed = 0, g=0;
            for(; g<num_genes; g++){
                if(back_f[f1].vals[g]+back_f[f2].vals[g] <= 0)
                    missed++;
                if(missed >= best)
                    break;
            }
            if(missed < best){
                best = missed;
                best_f1 = back_f[f1].fid;
                best_f2 = back_f[f2].fid;
            }
            total_checked+=g;
            g_f1+=g;
        }
        if(f1!=feature_considered-1)
            fprintf(fout,"%d\t%d\t%d\n",f1,g_f1/(feature_considered-f1-1),best);
    }
     fclose(fout);
    fprintf(stderr,"=========%ld=======\n",total_checked);
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.*****\n",best_f1,best_f2,best, feature_considered);
    t = clock() - t;
    fprintf (stderr, "baseline_unsorted find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

void histogram(int num_bucket){  //B: num_bucket
    clock_t t = clock();
    hist.resize(FEATURE_NUM);
    hist2.resize(FEATURE_NUM);
    int num_genes = back_f[0].vals.size();
    for(int f=0; f<FEATURE_NUM; f++){
        vector<val_type > cur_f = back_f[f].vals; 
        // vector<val_type > cur_f(back_f[f].vals.begin(), back_f[f].vals.end()); 
        sort(cur_f.begin(),cur_f.end());
        for(int ptr=0; ptr<num_bucket; ptr++){
            hist[f].push_back(cur_f[ptr*num_genes/num_bucket]);
        }
        hist[f].push_back(cur_f[num_genes-1]);
        
        for(int ptr=0; ptr<num_genes; ptr++){
            hist2[f].push_back(cur_f[ptr]);
        }
        //hist2[f].push_back(cur_f[num_genes-1]);
        //fprintf(stderr,"%lf %lf\n",hist[f][0],hist[f][num_bucket]);
    }
    t = clock() - t;
    fprintf (stderr, "histogram: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

void pruning_histogram(int num_bucket){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    
    /*each feature pair*/
    int best = num_genes;
    int best_f1,best_f2;
    long pruned=0;
    for(int f1=0; f1<FEATURE_NUM; f1++){  //FEATURE_NUM
        for(int f2=f1+1; f2<FEATURE_NUM; f2++){
            //estimation
            double best_mis=0;
            
            for(int i=num_bucket; i>=0; i--){  // or linear scan??
                //int j = lower_bound(hist[f2].begin(), hist[f2].end(), -hist[f1][i])-hist[f2].begin();
                int j= upper_bound(hist[f2].begin(), hist[f2].end(), -hist[f1][i])-hist[f2].begin();
                //double est = (-hist[f1][i]-*(hist[f2].begin()+j-1))/(*(hist[f2].begin()+j)-*(hist[f2].begin()+j-1));
                
                if(j>=1){
                    
                    double est = (j==11)? 0: (-hist[f1][i]-hist[f2][j-1])/(hist[f2][j]-hist[f2][j-1]);
                    double cur_mis = (j+i-1+est)/(double)(num_bucket)-1;
                    //double cur_mis = (j+i)/(double)(num_bucket)-1;
                    if(best_mis < cur_mis)
                        best_mis = cur_mis;
                    if(best_mis*num_genes >= best)
                        break;
                    
                    //if(est >=1 ||est<0){
                     //   fprintf(stderr, "error======%lf %lf %lf %lf %d\n",est,*(hist[f2].begin()+j-1), -hist[f1][i], *(hist[f2].begin()+j), j);
                     //   break;
                    //}
                }
                //if(j1-j<1)
                    //fprintf(stderr, "error======%d %d\n",j,j1);
                
            }
            
           /* for(int i=num_bucket; i>=num_bucket/2; i--){  // or linear scan??
                //int j = lower_bound(hist[f2].begin(), hist[f2].end(), -hist[f1][i])-hist[f2].begin();
                int j= upper_bound(hist2[f2].begin(), hist2[f2].end(), -hist[f1][i])-hist2[f2].begin();
                //double est = (-hist[f1][i]-*(hist[f2].begin()+j-1))/(*(hist[f2].begin()+j)-*(hist[f2].begin()+j-1));
                //if(j1-j<1)
                    //fprintf(stderr, "error======%d %d\n",j,j1);
                double cur_mis = (i)/(double)(num_bucket)+(j-1)/(double)(num_genes)-1;
                if(best_mis < cur_mis)
                    best_mis = cur_mis;
                if(best_mis*num_genes >= best)
                    break;
            }
            for(int i=num_bucket; i>=num_bucket/2; i--){  // or linear scan??
                //int j = lower_bound(hist[f2].begin(), hist[f2].end(), -hist[f1][i])-hist[f2].begin();
                int j= upper_bound(hist2[f1].begin(), hist2[f1].end(), -hist[f2][i])-hist2[f1].begin();
                //double est = (-hist[f1][i]-*(hist[f2].begin()+j-1))/(*(hist[f2].begin()+j)-*(hist[f2].begin()+j-1));
                //if(j1-j<1)
                    //fprintf(stderr, "error======%d %d\n",j,j1);
                double cur_mis = (i)/(double)(num_bucket)+(j-1)/(double)(num_genes)-1;
                if(best_mis < cur_mis)
                    best_mis = cur_mis;
                if(best_mis*num_genes >= best)
                    break;
            }*/
            
            
            if(best_mis*num_genes < best)
            {
                int missed = 0;
                for(int g=0; g<num_genes; g++){
                    if(back_f[f1].vals[g]+back_f[f2].vals[g] <= 0)
                        missed++;
                }
                if(missed < best){
                    best = missed;
                    best_f1 = back_f[f1].fid;
                    best_f2 = back_f[f2].fid;
                }
                //if(missed < best_mis*num_genes)
                  //  fprintf(stderr,"error: %lf %d %lf\n",best_mis*num_genes, missed, (best_mis-1/(double)(num_bucket))*num_genes);
            }
            else
                pruned++;
            
            //fprintf(fout,"%d\t%lf\t%lf\t%d\t%d\n",missed,est*num_genes,est,out_f1,out_f2);
        }
            
    }
    fprintf (stderr, "******best feature pair (%d,%d:%d):%ld pruned.*****\n",best_f1,best_f2,best, pruned);
    t = clock() - t;
    fprintf (stderr, "histogram profile: find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

}



void pruning_histogram_trick(int num_bucket){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    
    /*each feature pair*/
    int best = num_genes;
    int best_f1,best_f2;
    long pruned=0;
    for(int f1=0; f1<FEATURE_NUM; f1++){  //FEATURE_NUM
        for(int f2=f1+1; f2<FEATURE_NUM; f2++){
            //estimation
            double best_mis=0;
            
            for(int i=num_bucket; i>=0; i--){  // or linear scan??
                //int j = lower_bound(hist[f2].begin(), hist[f2].end(), -hist[f1][i])-hist[f2].begin();
                int j= upper_bound(hist[f2].begin(), hist[f2].end(), -hist[f1][i])-hist[f2].begin();
                //if(j1-j<1)
                    //fprintf(stderr, "error======%d %d\n",j,j1);
                double cur_mis = (j+i-1)/(double)(num_bucket)-1;
                if(best_mis < cur_mis)
                    best_mis = cur_mis;
                if(best_mis*num_genes >= best)
                    break;
            }
            
            if(best_mis*num_genes < best)
            {
                int missed = 0;
                for(int g=0; g<num_genes; g++){
                    if(back_f[f1].vals[g]+back_f[f2].vals[g] <= 0)
                        missed++;
                    if(missed>=best)
                        break;
                }
                if(missed < best){
                    best = missed;
                    best_f1 = back_f[f1].fid;
                    best_f2 = back_f[f2].fid;
                }
                //if(missed < best_mis*num_genes)
                  //  fprintf(stderr,"error: %lf %d %lf\n",best_mis*num_genes, missed, (best_mis-1/(double)(num_bucket))*num_genes);
            }
            else
                pruned++;
            
            //fprintf(fout,"%d\t%lf\t%lf\t%d\t%d\n",missed,est*num_genes,est,out_f1,out_f2);
        }
            
    }
    fprintf (stderr, "******best feature pair (%d,%d:%d):%ld pruned.*****\n",best_f1,best_f2,best, pruned);
    t = clock() - t;
    fprintf (stderr, "histogram profile: find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

}

void pruning_histogram_sortG(int num_bucket){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    
    /*each feature pair*/
    int best = num_genes;
    int best_f1,best_f2;
    int feature_considered =  FEATURE_NUM;
    long total_checked=0, pruned=0;
    for(int f1=0; f1<feature_considered; f1++){
        for(int f2=f1+1; f2<feature_considered; f2++){
            double best_mis=0;
            
            for(int i=num_bucket; i>=0; i--){  // or linear scan??
                //int j = lower_bound(hist[f2].begin(), hist[f2].end(), -hist[f1][i])-hist[f2].begin();
                int j= upper_bound(hist[f2].begin(), hist[f2].end(), -hist[f1][i])-hist[f2].begin();
                //if(j1-j<1)
                    //fprintf(stderr, "error======%d %d\n",j,j1);
                double cur_mis = (j+i-1)/(double)(num_bucket)-1;
                if(best_mis < cur_mis)
                    best_mis = cur_mis;
                if(best_mis*num_genes >= best)
                    break;
            }
            
            if(best_mis*num_genes < best){
                int missed = 0,g=0;
                for(; g<num_genes; g++){
                    if(sorted_f[f1].vals[g]+sorted_f[f2].vals[g] <= 0)
                        missed++;
                    if(missed >= best)
                        break;
                }
                if(missed < best){
                    best = missed;
                    best_f1 = sorted_f[f1].fid;
                    best_f2 = sorted_f[f2].fid;
                }
                total_checked+=g;
            }
            else
                pruned++;
        }
    }
    fprintf(stderr,"=========%ld %d=======\n",total_checked, pruned);
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.*****\n",best_f1,best_f2,best, feature_considered);
    t = clock() - t;
    fprintf (stderr, "baseline_unsorted_sortG find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}


void estimation_profile(int num_bucket){  //profile the estimation of each feature pair score
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    vector<vector<int> > heatmap;
    
    heatmap.resize(num_bucket*num_bucket);
    for(int i=0; i<num_bucket*num_bucket; i++){
        for(int j=0; j<num_bucket*num_bucket; j++)
            heatmap[i].push_back(0);
    }

    FILE* fout = fopen("stat_estimation.txt","w");
    /*each feature pair*/
    int best = num_genes;
    int best_f1,best_f2;
    for(int f1=0; f1<FEATURE_NUM; f1++){
        for(int f2=f1+1; f2<FEATURE_NUM; f2++){
            //estimation
            vector<val_type>::iterator ptr_f1, ptr_f2; 
            ptr_f2 =upper_bound (hist[f2].begin(), hist[f2].end(), -hist[f1][0]);
            ptr_f1 = upper_bound (hist[f1].begin(), hist[f1].end(), -hist[f2][0]);
            
            int out_f1 = (ptr_f1 ==hist[f1].end())?1:0;;
            int out_f2 = (ptr_f2 ==hist[f2].end())?1:0;
            int h=ptr_f2-hist[f2].begin();
            int b=ptr_f1-hist[f1].begin();
            
            double est = (h-1)*(b-1)/(double)(2*num_bucket*num_bucket);
            
            int missed = 0;
            for(int g=0; g<num_genes; g++){
                if(back_f[f1].vals[g]+back_f[f2].vals[g] <= 0)
                    missed++;
            }
            
            int ptr1= ceil((double)((h-1)*(b-1))/2)-1;  // estimate
            int ptr2 = ceil((double)(missed*(num_bucket*num_bucket))/num_genes)-1;
            heatmap[ptr1][ptr2]++;
            
            if(missed < best){
                best = missed;
                best_f1 = back_f[f1].fid;
                best_f2 = back_f[f2].fid;
            }
            //if(rand()%100000<100)
            //    fprintf(fout,"%d\t%d\n",missed,(int)(est*num_genes));
            //fprintf(fout,"%d\t%lf\t%lf\t%d\t%d\n",missed,est*num_genes,est,out_f1,out_f2);
        }
    }
    
    fprintf(fout, "est\t"); //( ]
    for(int i=0; i< num_bucket*num_bucket; i++){
        fprintf(fout, "%dreal\t", (i+1)*num_genes/(num_bucket*num_bucket)); //( ]
    }
    fprintf(fout,"\n");
    
    for(int i=0; i<num_bucket*num_bucket; i++){
        fprintf(fout, "%dest", (i+1)*num_genes/(num_bucket*num_bucket)); //( ]
        for(int j=0; j<num_bucket*num_bucket; j++){
            fprintf(fout, "\t%d", heatmap[i][j]);
        }
        fprintf(fout,"\n");
    }
    fclose(fout);
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.*****\n",best_f1,best_f2,best, FEATURE_NUM);
    t = clock() - t;
    fprintf (stderr, "histogram profile: find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}


void baseline_unsorted_withbitmap(){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
     FILE* fout = fopen("geneNum_unsorted_bitmap.txt","w");
    //each feature pair
    int best = num_genes;
    int best_f1,best_f2;
    int feature_considered =  FEATURE_NUM, pruned = 0;
    long total_checked=0;

    for(int f1=0; f1<feature_considered; f1++){
        int g_f1=0,no_pruned=0;
        for(int f2=f1+1; f2<feature_considered; f2++){
            int upper_num = (back_f[f1].bitmap | back_f[f2].bitmap).count(); //(bitmaps[f1] | bitmaps[f2]).count();
            if(upper_num > (num_genes- best)){
                no_pruned++;
                int missed = 0,g=0;
                for(; g<num_genes; g++){

                    if(back_f[f1].vals[g]+back_f[f2].vals[g] <= 0)
                        missed++;
                    if(missed >= best)
                        break;
                }
                if(missed < best){
                    best = missed;
                    best_f1 = back_f[f1].fid;
                    best_f2 = back_f[f2].fid;
                }
                g_f1 +=g;
                total_checked +=g;
            }
            else
                pruned++;
        }
        if(no_pruned!=0)
            fprintf(fout,"%d\t%d\t%d\t%d\n",f1,g_f1/no_pruned,no_pruned,best);
        else
            fprintf(fout,"%d\t%d\t%d\t%d\n",f1,g_f1,no_pruned,best);
    }
     fclose(fout);
    fprintf(stderr,"=========%ld=======\n",total_checked);
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.prunned %d*****\n",best_f1,best_f2,best, feature_considered, pruned);
    t = clock() - t;
    fprintf (stderr, "baseline_nosort_withbitmap find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

void baseline_unsorted_bitmapoverhead(){
    clock_t t = clock();
    int t_bitmap = 0; 
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    
    int feature_considered =  FEATURE_NUM;

    for(int f1=0; f1<feature_considered; f1++){
        for(int f2=f1+1; f2<feature_considered; f2++){
            int test = (back_f[f1].bitmap | back_f[f2].bitmap).count(); //(bitmaps[f1] | bitmaps[f2]).count();
        }
    }
 
    t = clock() - t;
    fprintf (stderr, "baseline_nosort_bitmapoverhead find best feature pair: %d clicks (%f seconds).\n", t,((float)t)/CLOCKS_PER_SEC);
}


/*sorted features*/
void baseline_horizontal(int num_bucket){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    sort(back_f.begin(),back_f.end(),compareByCorrectD);
    fprintf(stderr,"%d correct classified genes in top1 feature %d\n",back_f[0].correct,back_f[0].fid);
    FILE* fout = fopen("geneNum_horizontal.txt","w");
    /*each feature pair*/
    int best = num_genes;
    int best_f1,best_f2;
    int feature_considered =  FEATURE_NUM;
    long total_checked=0, pruned=0;
    for(int f1=0; f1<feature_considered; f1++){
        int g_f1=0;
        for(int f2=f1+1; f2<feature_considered; f2++){
            double best_mis=0;
            int f1id=back_f[f1].fid, f2id=back_f[f2].fid;
            for(int i=num_bucket; i>=0; i--){  // or linear scan??
                int j= upper_bound(hist[f2id].begin(), hist[f2id].end(), -hist[f1id][i])-hist[f2id].begin();
                double cur_mis = (j+i-1)/(double)(num_bucket)-1;
                if(best_mis < cur_mis)
                    best_mis = cur_mis;
                if(best_mis*num_genes >= best)
                    break;
            }
            
            if(best_mis*num_genes < best){
                int missed = 0, g=0;
                for(; g<num_genes; g++){
                    if(back_f[f1].vals[g]+back_f[f2].vals[g] <= 0)
                        missed++;
                    if(missed >= best)
                        break;
                }
                if(missed < best){
                    best = missed;
                    best_f1 = back_f[f1].fid;
                    best_f2 = back_f[f2].fid;
                }
                total_checked+=g;
                g_f1+=g;
            }
            else
                pruned++;
        }
       // if(f1!=feature_considered-1)
         //   fprintf(fout,"%d\t%d\t%d\t%d\n",f1,back_f[f1].fid,g_f1/(feature_considered-f1-1),best);
    }
     fclose(fout);
    fprintf(stderr,"=========%ld %ld=======\n",total_checked,pruned);
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.*****\n",best_f1,best_f2,best, feature_considered);
    t = clock() - t;
    fprintf (stderr, "baseline find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}


void baseline_vertical(){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    sort(back_f.begin(),back_f.end(),compareByCorrectD);
    FILE* fout = fopen("geneNum_vertical.txt","w");
    //each feature pair
    int best = num_genes;
    int best_f1,best_f2;
    int feature_considered =  FEATURE_NUM, pruned = 0;
    long total_checked=0;
    for(int f2=0; f2<feature_considered; f2++){
        int g_f2=0;
        for(int f1=0; f1<f2; f1++){
            int missed = 0,g=0;
            for(; g<num_genes; g++){
                if(back_f[f1].vals[g]+back_f[f2].vals[g] <= 0)
                    missed++;
                if(missed >= best)
                    break;
            }
            if(missed < best){
                best = missed;
                best_f1 = back_f[f1].fid;
                best_f2 = back_f[f2].fid;
            }
            g_f2+=g;
            total_checked+=g;
        }
        if(f2!=0)
            fprintf(fout,"%d\t%d\t%d\t%d\n",f2,back_f[f2].fid,g_f2/f2,best);
    }
     fclose(fout);
    fprintf(stderr,"=========%ld=======\n",total_checked);
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.*****\n",best_f1,best_f2,best, feature_considered);
    t = clock() - t;
    fprintf (stderr, "baseline_vertical find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

void baseline_horizontal_bitmap(){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    sort(back_f.begin(),back_f.end(),compareByCorrectD);
    FILE* fout = fopen("geneNum_horizontal_bitmap.txt","w");
    //each feature pair
    int best = num_genes;
    int best_f1,best_f2;
    int feature_considered =  FEATURE_NUM, pruned = 0;
    long total_checked=0;
    for(int f1=0; f1<feature_considered; f1++){
        int g_f1=0,no_pruned=0;
        for(int f2=f1+1; f2<feature_considered; f2++){
            //int i=back_f[f1].fid;
            //int j=back_f[f2].fid;
            
            if((back_f[f1].bitmap | back_f[f2].bitmap).count()> (num_genes- best)){
            //if((bitmaps[i] | bitmaps[j]).count()> (num_genes- best)){
                no_pruned++;
                int missed = 0,g=0;
                for(; g<num_genes; g++){
                    if(back_f[f1].vals[g]+back_f[f2].vals[g] <= 0)
                        missed++;
                    if(missed >= best)
                        break;
                }
                if(missed < best){
                    best = missed;
                    best_f1 = back_f[f1].fid;
                    best_f2 = back_f[f2].fid;
                }
                g_f1+=g;
                total_checked+=g;
            }
            else
                pruned++;
            
        }
        if(no_pruned!=0)
            fprintf(fout,"%d\t%d\t%d\t%d\n",f1,g_f1/no_pruned,no_pruned,best);
        else
            fprintf(fout,"%d\t%d\t%d\t%d\n",f1,g_f1,no_pruned,best);
    }
    fclose(fout);
    fprintf(stderr,"=========%ld=======\n",total_checked);
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.prunned %d*****\n",best_f1,best_f2,best, feature_considered, pruned);
    t = clock() - t;
    fprintf (stderr, "baseline_horizontal_bitmap find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

void baseline_horizontal_bitmapoverhead(){
    clock_t t = clock();
    
    int feature_considered =  FEATURE_NUM, pruned = 0;
    long total_checked=0;
    for(int f1=0; f1<feature_considered; f1++){
        for(int f2=f1+1; f2<feature_considered; f2++){

            int test = (back_f[f1].bitmap | back_f[f2].bitmap).count(); //(bitmaps[i] | bitmaps[j]).count();
        }
    }
    t = clock() - t;
    fprintf (stderr, "baseline_horizontal_bitmapoverhead find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}




/*sort genes according to global order*/
void baseline_unsorted_sortG(){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    
    FILE* fout = fopen("geneNum_unsorted_sortG.txt","w");
    /*each feature pair*/
    int best = num_genes;
    int best_f1,best_f2;
    int feature_considered =  FEATURE_NUM;
    long total_checked=0;
    for(int f1=0; f1<feature_considered; f1++){
        int g_f1=0;
        for(int f2=f1+1; f2<feature_considered; f2++){
            int missed = 0,g=0;
            for(; g<num_genes; g++){
                if(sorted_f[f1].vals[g]+sorted_f[f2].vals[g] <= 0)
                    missed++;
                if(missed >= best)
                    break;
            }
            if(missed < best){
                best = missed;
                best_f1 = sorted_f[f1].fid;
                best_f2 = sorted_f[f2].fid;
            }
            total_checked+=g;
            g_f1+=g;
        }
        if(f1!=feature_considered-1)
            fprintf(fout,"%d\t%d\t%d\n",f1,g_f1/(feature_considered-f1-1),best);
    }
    fclose(fout);
    fprintf(stderr,"=========%ld=======\n",total_checked);
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.*****\n",best_f1,best_f2,best, feature_considered);
    t = clock() - t;
    fprintf (stderr, "baseline_unsorted_sortG find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

void baseline_unsorted_sortG_bitmap(){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    
    FILE* fout = fopen("geneNum_unsorted_sortG_bitmap.txt","w");
    /*each feature pair*/
    int best = num_genes;
    int best_f1,best_f2;
    int feature_considered =  FEATURE_NUM, pruned=0;
    long total_checked=0;
    for(int f1=0; f1<feature_considered; f1++){
        int g_f1=0,no_pruned=0;
        for(int f2=f1+1; f2<feature_considered; f2++){
            
            if((sorted_f[f1].bitmap | sorted_f[f2].bitmap).count()> (num_genes- best)){
            //if((bitmaps[i] | bitmaps[j]).count()> (num_genes- best)){
                no_pruned++;
               int missed = 0, g=0;
                for(; g<num_genes; g++){
                    if(sorted_f[f1].vals[g]+sorted_f[f2].vals[g] <= 0)
                        missed++;
                    if(missed >= best)
                        break;
                }
                if(missed < best){
                    best = missed;
                    best_f1 = sorted_f[f1].fid;
                    best_f2 = sorted_f[f2].fid;
                }
                total_checked+=g;
                g_f1+=g;
            }
            else
                pruned++;
        }
        if(no_pruned!=0)
            fprintf(fout,"%d\t%d\t%d\t%d\n",f1,g_f1/no_pruned,no_pruned,best);
        else
            fprintf(fout,"%d\t%d\t%d\t%d\n",f1,g_f1,no_pruned,best);
    }
    fclose(fout);
    fprintf(stderr,"=========%ld=======\n",total_checked);
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.%d prunned*****\n",best_f1,best_f2,best, feature_considered,pruned);
    t = clock() - t;
    fprintf (stderr, "baseline_unsorted_sortG_bitmap find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}


void baseline_unsorted_sortG_bitmap_overhead(){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    
    
    for(int f1=0; f1<FEATURE_NUM; f1++){
        for(int f2=f1+1; f2<FEATURE_NUM; f2++){
            int test = (sorted_f[f1].bitmap | sorted_f[f2].bitmap).count();
        }
    }
    t = clock() - t;
    fprintf (stderr, "baseline_unsorted_sortG_bitmap_overhead find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

void baseline_horizontal_sortG(){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    sort(sorted_f.begin(),sorted_f.end(),compareByCorrectD);
    fprintf(stderr,"%d correct classified genes in top1 feature %d\n",sorted_f[0].correct,sorted_f[0].fid);
    FILE* fout = fopen("geneNum_horizontal_sortG.txt","w");
    /*each feature pair*/
    int best = num_genes;
    int best_f1,best_f2;
    int feature_considered =  FEATURE_NUM;
    long total_checked=0;
    for(int f1=0; f1<feature_considered; f1++){
        int g_f1=0;
        for(int f2=f1+1; f2<feature_considered; f2++){
            int missed = 0, g=0;
            for(; g<num_genes; g++){
                if(sorted_f[f1].vals[g]+sorted_f[f2].vals[g] <= 0)
                    missed++;
                if(missed >= best)
                    break;
            }
            if(missed < best){
                best = missed;
                best_f1 = sorted_f[f1].fid;
                best_f2 = sorted_f[f2].fid;
            }
            total_checked+=g;
            g_f1+=g;
        }
        if(f1!=feature_considered-1)
            fprintf(fout,"%d\t%d\t%d\t%d\n",f1,sorted_f[f1].fid,g_f1/(feature_considered-f1-1),best);
    }
    fclose(fout);
    fprintf(stderr,"=========%ld=======\n",total_checked);
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.*****\n",best_f1,best_f2,best, feature_considered);
    t = clock() - t;
    fprintf (stderr, "baseline find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}




void baseline_vertical_sortG(){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    sort(sorted_f.begin(),sorted_f.end(),compareByCorrectD);
    FILE* fout = fopen("geneNum_vertical_sortG.txt","w");
    //each feature pair
    int best = num_genes;
    int best_f1,best_f2;
    int feature_considered =  FEATURE_NUM, pruned = 0;
    long total_checked=0;
    for(int f2=0; f2<feature_considered; f2++){
        int g_f2=0;
        for(int f1=0; f1<f2; f1++){
            int missed = 0,g=0;
            for(; g<num_genes; g++){
                if(sorted_f[f1].vals[g]+sorted_f[f2].vals[g] <= 0)
                    missed++;
                if(missed >= best)
                    break;
            }
            if(missed < best){
                best = missed;
                best_f1 = sorted_f[f1].fid;
                best_f2 = sorted_f[f2].fid;
            }
            g_f2+=g;
            total_checked+=g;
        }
        if(f2!=0)
            fprintf(fout,"%d\t%d\t%d\t%d\n",f2,sorted_f[f2].fid,g_f2/f2,best);
    }
    fclose(fout);
    fprintf(stderr,"=========%ld=======\n",total_checked);
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.*****\n",best_f1,best_f2,best, feature_considered);
    t = clock() - t;
    fprintf (stderr, "baseline_vertical find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}



void baseline_horizontal_sortG_withbitmap(){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    fprintf(stderr,"number of genes in this exp: %d \n", num_genes);
    sort(sorted_f.begin(),sorted_f.end(),compareByCorrectD);
    FILE* fout = fopen("geneNum_horizontal_sortG_bitmap.txt","w");
    
    //each feature pair
    int best = num_genes;
    int best_f1,best_f2;
    int feature_considered =  FEATURE_NUM, pruned = 0;
    long total_checked=0;
    for(int f1=0; f1<feature_considered; f1++){
        int g_f1 =0, no_pruned=0;;
        for(int f2=f1; f2<feature_considered; f2++){
            
            if((sorted_f[f1].bitmap | sorted_f[f2].bitmap).count()> (num_genes- best)){
            //if((bitmaps[i] | bitmaps[j]).count()> (num_genes- best)){
               int missed = 0, g=0;
               no_pruned++;
                for(; g<num_genes; g++){
                    if(sorted_f[f1].vals[g]+sorted_f[f2].vals[g] <= 0)
                        missed++;
                    if(missed >= best)
                        break;
                }
                if(missed < best){
                    best = missed;
                    best_f1 = sorted_f[f1].fid;
                    best_f2 = sorted_f[f2].fid;
                }
                total_checked+=g;
                g_f1+=g;
            }
            else
                pruned++;
            
        }
        if(no_pruned!=0)
            fprintf(fout,"%d\t%d\t%d\t%d\n",f1,g_f1/no_pruned,no_pruned,best);
        else
            fprintf(fout,"%d\t%d\t%d\t%d\n",f1,g_f1,no_pruned,best);
        
    }
    fprintf(stderr,"=========%ld=======\n",total_checked);
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered.prunned %d*****\n",best_f1,best_f2,best, feature_considered, pruned);
    t = clock() - t;
    fprintf (stderr, "baseline_horizontal_sortG_withbitmap find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}


void baseline_horizontal_sortG_withbitmap_overhead(){
    clock_t t = clock();
    int num_genes = back_f[0].vals.size();  // number of genes in this experiment
    
    for(int f1=0; f1<FEATURE_NUM; f1++){
        for(int f2=f1; f2<FEATURE_NUM; f2++){
            int test = (sorted_f[f1].bitmap | sorted_f[f2].bitmap).count();
        }
    }
    t = clock() - t;
    fprintf (stderr, "baseline_horizontal_sortG_withbitmap find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}




void test(){
    dynamic_bitset<> test1(50);
    fprintf(stderr,"=======%ld====\n",test1.size());
    vector<dynamic_bitset<> > test2;
    test2.push_back(dynamic_bitset<>(50));
    test2.push_back(dynamic_bitset<>(50));
    test2[0][1] = 1; test2[1][2] = 1;
    test1= test2[1]&test2[0];
    cout<< test1.count() << " " << test2[1].count() <<" "<< (test2[1]&test2[0]).count()<<endl;
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
   // calculate_50percent();
    histogram(10);
    //pruning_histogram(10);
    //pruning_histogram_trick(10);
    //sort_genes();
    //pruning_histogram_sortG(10);
    
    //baseline_horizontal(10);
    estimation_profile(10);
/*    baseline_unsorted_noprune();
    baseline_unsorted();
    baseline_unsorted_withbitmap();
    baseline_unsorted_bitmapoverhead();
    
    baseline_horizontal();
    baseline_vertical();
    baseline_horizontal_bitmap();
    baseline_horizontal_bitmapoverhead();*/
    
    
//    sort_genes_fpair();
    //
//    baseline_unsorted_sortG();
    //baseline_unsorted_sortG_bitmap();
    //baseline_unsorted_sortG_bitmap_overhead();
    /*baseline_horizontal_sortG();
    baseline_vertical_sortG();
    baseline_horizontal_sortG_withbitmap();
    baseline_horizontal_sortG_withbitmap_overhead();*/
  
    //sorted_list();
    //sorted_list_bi();
    //sorted_list_desc();
    t = clock() - t;
    fprintf (stderr, "It took me %d clicks (%f seconds) for the whole program.\n",t,((float)t)/CLOCKS_PER_SEC);
    return 0;
}


/*
//ascending
void sorted_list(){
    clock_t t = clock();
    int num_genes = sorted_f[0].size();  // number of genes in this experiment
    
    for(int i=0; i<FEATURE_NUM; i++){
        //sort each feature
        sort(sorted_f[i].begin(), sorted_f[i].end(), compareByvalue);  //ascending 
        //fprintf(stderr, "%lf %lf \n",sorted_f[f][5].val,sorted_f[f][6].val);
    }
    //each feature pair
    int best = num_genes;
    int best_f1,best_f2, pruned = 0, overall_check=0;
    int feature_considered =  FEATURE_NUM;
    for(int i=0; i<feature_considered; i++){
        int checked_total =0;
        for(int j=i; j<feature_considered; j++){
            int f1 = top_f[i].fid;
            int f2 = top_f[j].fid;
            if((bitmaps[f1] | bitmaps[f2]).count()> (num_genes- best))
            {
                int missed = 0;
                //bool * check = (bool *)malloc(sizeof(bool) * num_genes);
                //memset(check,false,sizeof(bool) * num_genes);
                double thres;
                for( int ptr2 =0;ptr2 < num_genes; ptr2++){
                    thres = sorted_f[f1][0].val + sorted_f[f2][ptr2].val;
                    if(thres >0){
                        pruned++;
                        break;
                    }
                    int cur_gid =  sorted_f[f2][ptr2].id;
                    // if(back_f[f2][cur_gid] != sorted_f[f2][ptr2].val)
                    //     fprintf (stderr, "error========================.\n");
                    if(back_f[f1][cur_gid] + sorted_f[f2][ptr2].val <=0)
                        missed++;
                    
                    if(missed >= best)
                        break;
                }
                //checked_total+=ptr2;
                if(missed < best){
                    best = missed;
                    best_f1 = f1;
                    best_f2 = f2;
                }
            }
            
        }
        //overall_check+=checked_total/feature_considered;
    }
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered. %d pruned, %d checked *****\n",best_f1,best_f2,best, feature_considered,pruned,overall_check/feature_considered);
    t = clock() - t;
    fprintf (stderr, "find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}


//descending
void sorted_list_desc(){
    clock_t t = clock();
    int num_genes = sorted_f[0].size();  // number of genes in this experiment
    
    for(int i=0; i<FEATURE_NUM; i++){
        //sort each feature
        sort(sorted_f[i].begin(), sorted_f[i].end(), compareByvalueDe);  //ascending ?descending?
        //fprintf(stderr, "%lf %lf \n",sorted_f[f][5].val,sorted_f[f][6].val);
    }
    //each feature pair
    int best = num_genes;
    int best_f1,best_f2, pruned = 0;
    int feature_considered =  FEATURE_NUM;
    for(int i=0; i<feature_considered; i++){
        for(int j=i; j<feature_considered; j++){
            int f1 = top_f[i].fid;
            int f2 = top_f[j].fid;
            int missed = 0;
            
            double thres;
            for(int ptr2 =0; ptr2 < num_genes; ptr2++){
                thres = sorted_f[f1][0].val + sorted_f[f2][ptr2].val;
                if(thres <= 0){
                    missed += (num_genes-ptr2);
                    pruned++;
                    break;
                }
                int cur_gid =  sorted_f[f2][ptr2].id;
                // if(back_f[f2][cur_gid] != sorted_f[f2][ptr2].val)
                //     fprintf (stderr, "error========================.\n");
                if(back_f[f1][cur_gid] + sorted_f[f2][ptr2].val <=0)
                    missed++;
                
                if(missed >= best)
                    break;
            }
            if(missed < best){
                best = missed;
                best_f1 = f1;
                best_f2 = f2;
            }
            
        }
    }
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered. %d pruned*****\n",best_f1,best_f2,best, feature_considered,pruned);
    t = clock() - t;
    fprintf (stderr, "find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}




//ascending bi direction
void sorted_list_bi(){
    clock_t t = clock();
    int num_genes = sorted_f[0].size();  // number of genes in this experiment
    
    for(int i=0; i<FEATURE_NUM; i++){
        //sort each feature
        sort(sorted_f[i].begin(), sorted_f[i].end(), compareByvalue);  //ascending ?descending?
    }
    //each feature pair
    int best = num_genes;
    int best_f1,best_f2, pruned = 0,overall_check=0;
    int feature_considered =  FEATURE_NUM;
    for(int i=0; i<feature_considered; i++){
        //fprintf (stderr, "%d th iteration\n",i);
        int checked_total =0;
        for(int j=i; j<feature_considered; j++){
            
            if((top_f[i].correct + top_f[j].correct) > (num_genes - best)){
                

                int f1 = top_f[i].fid;
                int f2 = top_f[j].fid;
                int missed = 0;
                bool * check = (bool *)malloc(sizeof(bool) * num_genes);
                memset(check,false,sizeof(bool) * num_genes);
                double thres;
                double delta1=100000001, delta2=100000000;
                int ptr1 = 0, ptr2 = 0, checked_num=0;
                while(ptr1 <  num_genes & ptr2 <  num_genes){  //ptr1 and ptr2 pointing to the entry not checked yet
                    if(delta1 > delta2){
                        int cur_gid =  sorted_f[f1][ptr1].id;
                        while(check[cur_gid] == true){
                            ptr1++;
                            cur_gid =  sorted_f[f1][ptr1].id;
                        }
                        if(sorted_f[f1][ptr1].val + back_f[f2][cur_gid] <=0){
                            missed++;
                        }
                        check[cur_gid] = true;
                        checked_num++;
                        ptr1++;
                        delta1 = sorted_f[f1][ptr1].val - sorted_f[f1][ptr1-1].val;
                    }
                    else{
                        int cur_gid =  sorted_f[f2][ptr2].id;
                        while(check[cur_gid] == true){
                            //fprintf(stderr,"%d %d %d %d\n",ptr1, ptr2, cur_gid, checked_num );
                            ptr2++;
                            cur_gid =  sorted_f[f2][ptr2].id;
                        }
                        if(sorted_f[f2][ptr2].val + back_f[f1][cur_gid] <=0){
                            missed++;
                        }
                        check[cur_gid] = true;
                        checked_num++;
                        ptr2++;
                        delta2 = sorted_f[f2][ptr2].val - sorted_f[f2][ptr2-1].val;
                    }
                    if(checked_num ==num_genes)
                    //if(missed >= best || checked_num ==num_genes)
                        break;
                    thres = sorted_f[f1][ptr1].val + sorted_f[f2][ptr2].val;
                    if(thres > 0){
                        pruned++;
                        break;
                    }
                    
                }
                checked_total +=checked_num;
                free(check);
                if(missed < best){
                    best = missed;
                    best_f1 = f1;
                    best_f2 = f2;
                }
            }
            
        }
        overall_check += checked_total/feature_considered;
        //fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered. %d pruned*****\n",best_f1,best_f2,best, feature_considered,pruned);
    }
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered. %d pruned, checked %d*****\n",best_f1,best_f2,best, feature_considered,pruned,overall_check/feature_considered);
    t = clock() - t;
    fprintf (stderr, "find best feature pair: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}*/
