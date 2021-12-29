#include"mtrand.h"
#include<iostream.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<fstream.h>
#include<time.h>
#include<assert.h>

#define SSIZE   50                          //Sample size.
#define    fr   1                           //Output all biallelic SNPs with minor allele count >=fr.
#define     L   500.0                       //Mean length of the geometrically distributed tract length in basepairs.
#define   LEN   100000                      //Length of the DNA sequence in basepairs = 1 + Number of potential recombination breakpoints.
#define THETA   100.0                       //Population mutation rate(i.e. 4Nu).
#define   RHO   100.0                       //Population crossing-over rate(i.e. 4Nr).
#define   GC    100.0                       //Population gene-conversion rate(i.e. 4Nc).
#define ghot    0                           //0 for uniform gene-conversion and 1 for gene-conversion hotspot.
#define gwidth  2000                        //Width of the conversion hotspot in basepairs.
#define gstart  24000                       //Start of the conversion hotspot in basepairs.
#define gfr     0.50                        //Fraction of conversion events occuring in hotspots(0 - 1.0).
#define rhot    0                           //0 for uniform crossing-over and 1 for crossing-over hotspot.
#define rwidth  2000                        //Width of the crossing-over hotspot in basepairs.
#define rstart  24000                       //Start of the crossing-over hotspot in basepairs.
#define rfrac   0.50                        //Fraction of crossing-over events occuring in hotspots(0 - 1.0).
#define RUNS    1                           //Total number of replicates.
#define ZERO    0                           //Constant. Do not change.    
#define MAXIMUM 100000 + int(1200*(RHO+GC)) //Maximum number of nodes in the graph. Increase for very high RHO and GAMMA.
#define RSIZE   (SSIZE-1)/31                //Reduced size in bitwise notation. Do not change.
#define MAXSNPS 1000 + int(25*THETA)        //Maximum number of SNPs in a chromosome.


unsigned long init[4] = {0x123, 0x234, 0x345, 0x456},length = 4;
MTRand drand;//double in [0, 1) generator, already init
MTRand_int32 irand(init, length);// 32-bit int generator
MTRand_int32 seed (init, length);// 32-bit int generator
MTRand mt;


//RANDOM FUNCTIONS
static inline double expon(double rate){double u;u = drand();return(-log(u)/rate);}
static inline int poisson(double mean){int count=0;double q,bound;
bound = exp(-mean);for(q=1;q>=bound;q*=drand()){++count;if(count>10000){return(int(mean));}}
return(count-1);}
static inline int choose(int i, int j){return(i + int((j-i)*drand()));}
static inline double  minimum(double i,double j){return(i>j) ? j:i;}
static inline double  maximum(double i,double j){return(i>j) ? i:j;}

//DEFINES A NODE IN THE GRAPH
class member
{   public:
    int size,(*p)[RSIZE+3];
    double times;
};
member tree[MAXIMUM];int current[500*SSIZE],matrix[MAXSNPS][RSIZE+2],polymorphism[LEN],repeats[LEN];

main(){ 
int i,j,k,l,m,v,w,x,y,mu,count,max,min,n,w1,w2,v1,v2,y1,mcount,rcount;
double a,b,c,d,e,f,q,r,times,times1,times2,sum,tot, gcount, summ, avgtime, theor;
int brk,stop,total,flag1,flag2,null[RSIZE+1];
 
double GAMMA;
GAMMA = GC*double(LEN + 50*L)/LEN;
double gfrac = gfr*GC/GAMMA;

//Initialize random seed
ifstream input("seed",ios::in);
while(!input.eof()){
for(i=0;i<length;i++){input>>init[i];cout<<init[i]<<" ";}cout<<"\n\n";break;}
input.close();
mt.seed(init,length);

ofstream output("seed",ios::out);
for(i=0;i<length;i++){output<<irand()<<" ";}
output.close();

memset(&null[0],'\0',4*(RSIZE+1));


for(n=RUNS; n > 0; --n){ 

//INITIALIZATION OF MEMBERS
memset(&polymorphism[0],'\0',4*LEN);
for(x=0;x<SSIZE;x++){
tree[x].p = new int[2][RSIZE+3];
for(i=0;i<RSIZE+3;i++){tree[x].p[0][i]=0;tree[x].p[1][i]=0;}
tree[x].size=2;tree[x].p[0][0]=0;tree[x].p[1][0]=LEN;
tree[x].p[0][1+x/31]|= (1<<(x%31));
tree[x].times=0;tree[x].p[0][RSIZE+2]=1;
}
for(x=0;x<500*SSIZE;x++){if(x<SSIZE){current[x]=x+1;}else{current[x]=0;}}

count=SSIZE;total=SSIZE;times=0;mcount=0;rcount=0;stop=0;gcount=0;    
//ORDER OF EVENTS TILL MRCA
while(count!=1){

//GENERATE TWO COMPETING EXPONENTIAL RANDOM VARIABLES
c  = expon((count*count-count)/2.0); 
r  = expon(count*(RHO+GAMMA)/2.0);

if(count>500*SSIZE){cout<<"Out of memory. Increase the size of the array named current in the top of the code.\n";return(0);}
if(total>MAXIMUM  ){cout<<"Out of memory. Increase the size of variable MAXIMUM in the top of the code.\n";return(0);}
if(mcount>MAXSNPS ){cout<<"Out of memory. Increase the size of variable MAXSNPS in the top of the code.\n";return(0);}
  
//IF THE NEXT EVENT IS COALESCENCE  
if(c<= r){
i = choose(1,count+1);j = choose(1,count+1);while(i==j){j = choose(1,count+1);}
if(i>j){max=i;min=j;}if(j>i){max=j;min=i;}
j = current[max-1] - 1;i = current[min-1] - 1;
--count;++total;times  = times + c;tree[total-1].times = times;
   
//UPDATE CURRENT   
current[min-1] = total;memmove(&current[max-1],&current[max],4*(count-max+1));

//MAKE A  NEW LINEAGE 
tree[total-1].p = new int[tree[i].size + tree[j].size + 1][RSIZE+3];w=0;w1=0;tree[total-1].size=0;
while(ZERO==0){
tree[total-1].p[tree[total-1].size][0]= minimum(tree[i].p[w][0],tree[j].p[w1][0]);
if(tree[total-1].p[tree[total-1].size][0]==tree[i].p[w][0]){++w;}
if(tree[total-1].p[tree[total-1].size][0]==tree[j].p[w1][0]){++w1;}
++tree[total-1].size;
if(w==tree[i].size){for(l=w1;l<tree[j].size;l++){tree[total-1].p[tree[total-1].size][0] = tree[j].p[l][0];++tree[total-1].size;}break;}
if(w1==tree[j].size){for(l=w;l<tree[i].size;l++){tree[total-1].p[tree[total-1].size][0] = tree[i].p[l][0];++tree[total-1].size;}break;}
}

q=0;
for(x=0;x<tree[i].size-1;x++){if(tree[i].p[x][RSIZE+2]==SSIZE){q += tree[i].p[x][0] - tree[i].p[x+1][0];}}
for(x=0;x<tree[j].size-1;x++){if(tree[j].p[x][RSIZE+2]==SSIZE){q += tree[j].p[x][0] - tree[j].p[x+1][0];}}

//COALESCE
w=0;w1=0;
for(x=0;x<tree[total-1].size-1;x++){
flag1=0;flag2=0; memcpy(&tree[total-1].p[x][1],&null[0],4*(RSIZE+1));  

for(y=w;y<tree[i].size-1;y++){
if(tree[i].p[y][0]<=tree[total-1].p[x][0] && tree[i].p[y+1][0]>=tree[total-1].p[x+1][0]){ 
flag1 = tree[i].p[y][RSIZE+2];if(flag1>0){memcpy(&tree[total-1].p[x][1],&tree[i].p[y][1],4*(RSIZE+1));}
w=y;break;}} 
   
for(y1=w1;y1<tree[j].size-1;y1++){
if(tree[j].p[y1][0]<=tree[total-1].p[x][0] && tree[j].p[y1+1][0]>=tree[total-1].p[x+1][0]){
flag2 = tree[j].p[y1][RSIZE+2];
if(flag1==0){if(flag2>0){memcpy(&tree[total-1].p[x][1],&tree[j].p[y1][1],4*(RSIZE+1));}}
else{if(flag2>0){for(l=1;l<RSIZE+2;l++){tree[total-1].p[x][l] = tree[total-1].p[x][l] | tree[j].p[y1][l];}}}
w1=y1;break;}}

tree[total-1].p[x][RSIZE+2]=flag1+flag2;
if(tree[total-1].p[x][RSIZE+2]==SSIZE){q += tree[total-1].p[x+1][0] - tree[total-1].p[x][0];}
}

//TRIM THE ENDS IF THEY HAVE REACHED MRCA
v2=tree[total-1].size-2;w2=0;
while(v2>=0){
l=tree[total-1].p[v2][RSIZE+2];
if(l>=1 && l<=SSIZE-1){break;}
if(l==0 || l==SSIZE){v2=v2-1;}}

while(w2<tree[total-1].size-1){
l=tree[total-1].p[w2][RSIZE+2];   
if(l>=1 && l<=SSIZE-1){break;} 
if(l==0 || l==SSIZE){w2=w2+1;}}
 
if(v2<w2){memmove(&current[min-1],&current[min],4*(count-min+1));delete tree[total-1].p;--count;}

//MERGE ALL-ANCESTRAL AND NON-ANCESTRAL INTERVALS 
if(v2>=w2){ 
w=0;x=w2; 
while(x<=v2){
l=tree[total-1].p[x][RSIZE+2];
if(l==0||l==SSIZE){
tree[total-1].p[w][0] = tree[total-1].p[x][0];
while(l==0||l==SSIZE){x=x+1;l=tree[total-1].p[x][RSIZE+2];}
for(y=1;y<RSIZE+3;y++){tree[total-1].p[w][y]=0;}
++w;}

if(l>=1 && l<=SSIZE-1){memcpy(&tree[total-1].p[w][0],&tree[total-1].p[x][0],4*(RSIZE+3));++w;++x;}
}
tree[total-1].p[w][0] = tree[total-1].p[v2+1][0];tree[total-1].size=w+1;}
   
//PUT MUTATIONS ACCORDING TO A POISSON PROCESS
times1 = tree[total-1].times - tree[i].times;
y = poisson(times1*THETA/2);
for(x=0;x<y;x++){ 
mu = int(LEN*drand())+1;while(polymorphism[mu-1]>0){mu = int(LEN*drand())+1;}
if(mu>tree[i].p[0][0] && mu<=tree[i].p[tree[i].size-1][0]){
for(k=1;k<tree[i].size;k++){if(mu<=tree[i].p[k][0]){break;}}
if(tree[i].p[k-1][RSIZE+2]>=fr && tree[i].p[k-1][RSIZE+2]<=SSIZE-fr){
if(polymorphism[mu-1]==0){
polymorphism[mu-1] = 1 + mcount;matrix[mcount][0]=mu;memcpy(&matrix[mcount][1],&tree[i].p[k-1][1],4*(RSIZE+1));
++mcount;}
}
}}

times2 = tree[total-1].times - tree[j].times;
y = poisson(times2*THETA/2);
for(x=0;x<y;x++){
mu = int(LEN*drand())+1;while(polymorphism[mu-1]>0){mu = int(LEN*drand())+1;}
if(mu>tree[j].p[0][0] && mu<=tree[j].p[tree[j].size-1][0]){
for(k=1;k<tree[j].size;k++){if(mu<=tree[j].p[k][0]){break;}}
if(tree[j].p[k-1][RSIZE+2]>=fr && tree[j].p[k-1][RSIZE+2]<=SSIZE-fr){
if(polymorphism[mu-1]==0){
polymorphism[mu-1] = 1 + mcount;matrix[mcount][0]=mu;memcpy(&matrix[mcount][1],&tree[j].p[k-1][1],4*(RSIZE+1));
++mcount;}
}
}}

//FREE MEMORY SPACE
delete tree[i].p;delete tree[j].p;

//CHECK IF ALL THE SEGMENTS HAVE REACHED THEIR MARGINAL MRCAs
if(q<0){"Fatal Error\n";return(0);}stop += q;
if(stop==LEN){for(x=0;x<count;x++){delete tree[current[x]-1].p;}break;} 
}
  
//IF THE NEXT EVENT IS RECOMBINATION
else{
//CHOOSE A LINEAGE AT RANDOM TO RECOMBINE AND UPDATE TIMES
j = choose(1,count+1);i = current[j-1] - 1;times += r;

double g,z;int start,end,tract;
g = drand();f = GAMMA/(RHO + GAMMA);brk = 2*LEN;

//IF NEXT EVENT IS GENE-CONVERSION
if(g < f){
tract = 1;z = drand(); while(z > 1.0/L){z = drand(); ++tract;}

//IF GENE-CONVERSION IS UNIFORM
if(ghot == 0){start = -50*L + int((LEN + 50*L - 1)*drand()) + 1;end = start + tract;}

//IF GENE-CONVERSION IS NON-UNIFORM
else{
z = drand();
if(z <= gfrac){start = gstart + int((gwidth - 1)*drand()) + 1; end = start + tract;}
else{
start = -50*L + int((LEN + 50*L - 1)*drand()) + 1;
while(start > gstart && start <= gstart + gwidth){start = -50*L + int((LEN + 50*L - 1)*drand()) + 1;}
end = start + tract;}
}

if(start >  tree[i].p[0][0] && start < tree[i].p[tree[i].size-1][0] && end >= tree[i].p[tree[i].size-1][0]){brk = start;}
if(start <= tree[i].p[0][0] && end   > tree[i].p[0][0]              && end <  tree[i].p[tree[i].size-1][0]){brk = end;}
}

//IF NEXT EVENT IS CROSSING-OVER
else{
//IF CROSSING-OVER IS UNIFORM
if(rhot == 0){brk = 1 + ((LEN-1)*drand());}

//IF CROSSING-OVER IS NON-UNIFORM
else{
z=drand();
if(z<=rfrac){ brk = rstart + int(rwidth*drand());}
else{brk = 1 + ((LEN-1)*drand());while(brk>=rstart && brk<rstart+rwidth){brk = 1 + ((LEN-1)*drand());}}
}}

//PERFORM CROSSING-OVER IF THE BREAKPOINT LIES WITHIN ANCESTRAL REGION
if(brk>tree[i].p[0][0] && brk<tree[i].p[tree[i].size-1][0]){

//PUT MUTATIONS ALONG THE BRANCHES
y = poisson((times - tree[i].times)*THETA/2);
for(x=0;x<y;x++){
mu = int(LEN*drand())+1;while(polymorphism[mu-1]>0){mu = int(LEN*drand())+1;}
if(mu>tree[i].p[0][0] && mu<=tree[i].p[tree[i].size-1][0]){
for(k=1;k<tree[i].size;k++){if(mu<=tree[i].p[k][0]){break;}}
if(tree[i].p[k-1][RSIZE+2]>=fr && tree[i].p[k-1][RSIZE+2]<=SSIZE-fr){
if(polymorphism[mu-1]==0){
polymorphism[mu-1] = 1 + mcount;matrix[mcount][0]=mu;memcpy(&matrix[mcount][1],&tree[i].p[k-1][1],4*(RSIZE+1));
++mcount;}
}
}}

++count;++total;tree[total-1].times = times;

//UPDATE CURRENT 
current[count-1] = total; 

//DETERMINE SIZE OF RECOMBINANTS
for(x=0;x<tree[i].size;x++){if(brk<=tree[i].p[x][0]){break;}}

if(brk<tree[i].p[x][0]){
w=x-1;v=x-1;
while(v>=0){
if(tree[i].p[v][RSIZE+2]>=1 && tree[i].p[v][RSIZE+2]<=SSIZE-1){break;}
else{--v;}}

while(w<tree[i].size-1){
if(tree[i].p[w][RSIZE+2]>=1 && tree[i].p[w][RSIZE+2]<=SSIZE-1){break;} 
else{++w;}}

if(tree[i].size-w<=0){cout<<"Error5\n";return(0);}if(v==-1){cout<<"Error6\n";return(0);}
tree[total-1].size = tree[i].size - w;tree[total-1].p = new int[tree[i].size - w][RSIZE+3];
 
//SPLIT THE LINEAGE TO 2 NEW ARRAYS
memcpy(&tree[total-1].p[0][0],&tree[i].p[w][0],4*(RSIZE+3)*tree[total-1].size);
if(w==x-1){tree[total-1].p[0][0] = brk;}
tree[i].size = v+2;tree[i].times = times;if(v==x-1){tree[i].p[x][0] = brk;}
}

else{
w=x;v=x-1;
while(v>=0){
if(tree[i].p[v][RSIZE+2]>=1 && tree[i].p[v][RSIZE+2]<=SSIZE-1){break;}
else{--v;}}

while(w<tree[i].size-1){
if(tree[i].p[w][RSIZE+2]>=1 && tree[i].p[w][RSIZE+2]<=SSIZE-1){break;}
else{++w;}}

if(tree[i].size-w==0){cout<<"Error5\n";return(0);}if(v==-1){cout<<"Fatal Error6\n";return(0);}
tree[total-1].size = tree[i].size - w;tree[total-1].p = new int[tree[i].size - w][RSIZE+3];
memcpy(&tree[total-1].p[0][0],&tree[i].p[w][0],4*(RSIZE+3)*tree[total-1].size);
tree[i].size = v+2;tree[i].times = times;}

if(g < f){++gcount;}else{++rcount;}
continue;}

//PERFORM GENE-CONVERSION IF BREAKPOINTS COVER ANCESTRAL REGION
if(g<f && start>tree[i].p[0][0] && end<tree[i].p[tree[i].size-1][0]){
int x1,x2;
//DETERMINE SIZE OF RECOMBINANTS
for(x1=0;x1<tree[i].size;x1++){if(start<=tree[i].p[x1][0]){break;}}
for(x2=x1;x2<tree[i].size;x2++){if(end<=tree[i].p[x2][0]){break;}}
 
v1=x1-1;if(tree[i].p[x1][0]==start){w1=x1;}else{w1=x1-1;}
while(v1>=0){
if(tree[i].p[v1][RSIZE+2]>=1 && tree[i].p[v1][RSIZE+2]<=SSIZE-1){break;}
else{--v1;}}

while(w1<tree[i].size-1){
if(tree[i].p[w1][RSIZE+2]>=1 && tree[i].p[w1][RSIZE+2]<=SSIZE-1){break;}
else{++w1;}}

v2=x2-1;if(tree[i].p[x2][0]==end){w2=x2;}else{w2=x2-1;}
while(v2>=0){
if(tree[i].p[v2][RSIZE+2]>=1 && tree[i].p[v2][RSIZE+2]<=SSIZE-1){break;}
else{--v2;}}

while(w2<tree[i].size-1){
if(tree[i].p[w2][RSIZE+2]>=1 && tree[i].p[w2][RSIZE+2]<=SSIZE-1){break;}
else{++w2;}}

//IF THE TRACT IS COMPLETELY NON-ANCESTRAL OR ALL-ANCESTRAL
if(v2<w1){continue;}

//IF THE TRACT CARRIES PARTLY ANCESTRAL MATERIAL
else{
++count;++total;++total;tree[total-1].times = times;tree[total-2].times=times;

//PUT MUTATIONS ALONG THE BRANCHES
y = poisson((times - tree[i].times)*THETA/2);
for(x=0;x<y;x++){
mu = int(LEN*drand())+1;while(polymorphism[mu-1]>0){mu = int(LEN*drand())+1;}
if(mu>tree[i].p[0][0] && mu<=tree[i].p[tree[i].size-1][0]){
for(k=1;k<tree[i].size;k++){if(mu<=tree[i].p[k][0]){break;}}
if(tree[i].p[k-1][RSIZE+2]>=fr && tree[i].p[k-1][RSIZE+2]<=SSIZE-fr){
if(polymorphism[mu-1]==0){
polymorphism[mu-1] = 1 + mcount;matrix[mcount][0]=mu;memcpy(&matrix[mcount][1],&tree[i].p[k-1][1],4*(RSIZE+1));
++mcount;}
}
}}
  
//UPDATE CURRENT
current[count-1] = total;current[j-1]=total-1;

//MAKE NEW LINEAGES
tree[total-2].size = v1 + 2 + tree[i].size - w2;tree[total-2].p = new int[tree[total-2].size][RSIZE+3];
memcpy(&tree[total-2].p[0][0],&tree[i].p[0][0],4*(RSIZE+3)*(v1+2));
if(v1==x1-1){tree[total-2].p[v1+1][0]= start;}tree[total-2].p[v1+1][RSIZE+2]=0;
memcpy(&tree[total-2].p[v1+2][0],&tree[i].p[w2][0],4*(RSIZE+3)*(tree[i].size-w2));
if(w2==x2-1){tree[total-2].p[v1+2][0]= end;}

tree[total-1].size = v2-w1+2;tree[total-1].p = new int[tree[total-1].size][RSIZE+3]; 
memcpy(&tree[total-1].p[0][0],&tree[i].p[w1][0],4*(RSIZE+3)*tree[total-1].size);
if(w1==x1-1){tree[total-1].p[0][0]=start;}
if(v2==x2-1){tree[total-1].p[tree[total-1].size-1][0]=end;}
delete tree[i].p;++gcount;continue;}
}}
}


int position[mcount], abc, counter;
cout<<"SNP count: "<<mcount<<"\n\n";
 
abc = 0;
for(l=0;l<LEN;l++){
if(polymorphism[l]>0){
position[abc] = l+1;

counter = 1;
if(l+1!=matrix[polymorphism[l]-1][0]){cout<<"ERR\n";return(0);}
cout<<l+1<<"\t";
     
for(x=1;x<RSIZE+1;x++){
w=matrix[polymorphism[l]-1][x];
for(m=0;m<31;m++){
cout<<(w%2);++counter;w=w/2;}}

w=matrix[polymorphism[l]-1][x];for(m=0;m<=(SSIZE-1)%31;m++){
cout<<(w%2);++counter;w=w/2;}
cout<<"\n";
abc++;assert(counter == SSIZE+1);}
}

assert(abc  ==  mcount);

}
}



 



 


 














