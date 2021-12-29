#include"mtrand.h"
#include<iostream.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<fstream.h>
#include<time.h>


#define SSIZE   50                                           //Sample size.
#define fr      2                                            //Minor allele count >=fr. 2 means only non-singleton SNPs will be output.
#define L       500.0                                        //Mean length of the geometrically distributed tract length in basepairs.
#define LEN     50000                                        //Length of the DNA sequence in basepairs.
#define ghot    0                                            //0 for uniform gene-conversion and 1 for gene-conversion hotspot.
#define gstart  24000                                        //Start of the conversion hotspot in basepairs.
#define gwidth  2000                                         //Width of the conversion hotspot in basepairs.
#define gfrac   0.50                                         //Fraction of conversion events in hotspots.(0-1.0)
#define rhot    0                                            //0 for uniform crossing-over and 1 for crossing-over hotspot.
#define rstart  24000                                        //Start of the crossing-over hotspot in basepairs.
#define rwidth  2000                                         //Width of the crossing-over hotspot in basepairs.
#define rfrac   0.50                                         //Fraction of crossing-over events in hotspots.(0-1.0)
#define RUNS    1                                            //Total number of replicates.
#define ZERO    0                                            //Constant. Do not change.    
#define RSIZE   (SSIZE-1)/31                                 //Reduced size in bitwise notation. Do not change.
double RHO;                                                  //Population crossing-over rate.(i.e. 4Nr)
#define rgrid 1                                              //Number of RHO values in the grid.
double rh[rgrid]  = {500.0};                                 //RHO values in ascending order.      
double GAMMA;                                                //Population gene-conversion rate.(i.e.4Nc)
#define ggrid 1                                              //Number of GAMMA values in the grid.
double gm[ggrid]  = {500.0};                                 //GAMMA values in ascending order.     
#define sites 100                                            //Number of biallelic SNPs in the dataset.
#define MAXIMUM 80000000                                     //Maximum number of nodes in the graph. Increase for very high RHO and GAMMA.
//Set these values at the start of main()
int position[sites];                                         //Positions of the SNPs in basepairs in ascending order.(1-LEN) 
int lfs[sites],ufs[sites];                                   //Simulate datasets with minor allele counts for SNPs between lfs and ufs.

unsigned long init[4] = {0x123, 0x234, 0x345, 0x456},length = 4;
MTRand drand;//double in [0, 1) generator, already init
MTRand_int32 irand(init, length);// 32-bit int generator
MTRand_int32 seed (init, length);// 32-bit int generator
MTRand mt;
double pairwise[sites][sites],prob[sites][2*SSIZE],prob1[sites][SSIZE][SSIZE],prob2[sites][SSIZE][SSIZE];
int current[500*SSIZE],matrix[sites][SSIZE],ccount[sites],all[sites][2*SSIZE][RSIZE+1];

//RANDOM FUNCTIONS
static inline double expon(double rate){double u;u = drand();return(-log(u)/rate);}
static inline int choose(int i, int j){return(i + int((j-i)*drand()));}
static inline double  minimum(double i,double j){return(i>j) ? j:i;}
static inline double  maximum(double i,double j){return(i>j) ? i:j;}

//DEFINES A NODE IN THE GRAPH
class member{   
public:
int size,(*p)[RSIZE+3];
double time;};member tree[MAXIMUM];

main(){ 
int i,j,k,l,m,q,v,w,x,y,mu,count,max,min,n,w1,w2,v1,v2,y1,maf1,maf2;
double a,b,c,d,e,f,r,time,time1,time2;int brk,stop,total,flag1,flag2,n1,n2;

//Specify positions and minor allele frequency ranges here. 
for(i=0;i<sites;i++){lfs[i]=fr;ufs[i]=SSIZE-fr;position[i]=(i+1)*LEN/sites;}

ifstream input("seed",ios::in);
while(!input.eof()){
for(i=0;i<length;i++){input>>init[i];cout<<init[i]<<" ";}cout<<"\n\n";break;}
input.close();
mt.seed(init,length);

ofstream output("seed",ios::out);
for(i=0;i<length;i++){output<<irand()<<" ";}
output.close();
 
for(n1=0;n1<rgrid;n1++){
for(n2=0;n2<ggrid;n2++){
RHO = rh[n1];GAMMA = gm[n2];n = 0;

while(n<RUNS){
for(i=0;i<sites;i++){
ccount[i]=0;
for(j=0;j<2*SSIZE;j++){prob[i][j]=0;
if(j<SSIZE){for(k=0;k<SSIZE;k++){prob1[i][j][k]=0;prob2[i][j][k]=0;}}
}}
  
//INITIALIZATION OF MEMBERS
for(x=0;x<SSIZE;x++){
tree[x].p = new int[2][RSIZE+3];
for(i=0;i<RSIZE+3;i++){tree[x].p[0][i]=0;tree[x].p[1][i]=0;}
tree[x].size=2;tree[x].p[0][0]=0;tree[x].p[1][0]=LEN;
tree[x].p[0][1+x/31]|= (1<<(x%31));
tree[x].time=0;tree[x].p[0][RSIZE+2]=1;
}

for(x=0;x<500*SSIZE;x++){if(x<SSIZE){current[x]=x+1;}else{current[x]=0;}}count=SSIZE;total=SSIZE;time=0;stop=0;
    
//ORDER OF EVENTS TILL MRCA
while(count!=1){
c  = expon((count*count-count)/2.0); 
r  = expon(count*(RHO+GAMMA)/2.0);

if(count>500*SSIZE){cout<<"Error1\n";return(0);}
if(total>MAXIMUM){cout<<"Error2\n";return(0);}
  
//IF THE NEXT EVENT IS COALESCENCE  
if(c<= r){

i = choose(1,count+1);j = choose(1,count+1);
while(i==j){j = choose(1,count+1);}
if(i>j){max=i;min=j;}if(j>i){max=j;min=i;}
j = current[max-1] - 1;i = current[min-1] - 1;
--count;++total;time  = time + c;tree[total-1].time = time;
   
//UPDATE CURRENT   
current[min-1] = total;memmove(&current[max-1],&current[max],4*(count-max+1));

//MAKE A  NEW LINEAGE 
tree[total-1].p = new int[tree[i].size + tree[j].size + 1][RSIZE+3];w=0;w1=0;tree[total-1].size=0;
while(ZERO==0){
tree[total-1].p[tree[total-1].size][0]= int(minimum(tree[i].p[w][0],tree[j].p[w1][0]));
if(tree[total-1].p[tree[total-1].size][0]==tree[i].p[w][0]){++w;}
if(tree[total-1].p[tree[total-1].size][0]==tree[j].p[w1][0]){++w1;}
tree[total-1].size++;
if(w==tree[i].size){for(l=w1;l<tree[j].size;l++){tree[total-1].p[tree[total-1].size][0] = tree[j].p[l][0];tree[total-1].size++;}break;}
if(w1==tree[j].size){for(l=w;l<tree[i].size;l++){tree[total-1].p[tree[total-1].size][0] = tree[i].p[l][0];tree[total-1].size++;}break;}
}

q=0;
for(x=0;x<tree[i].size-1;x++){if(tree[i].p[x][RSIZE+2]==SSIZE){q = q + tree[i].p[x][0] - tree[i].p[x+1][0];}}
for(x=0;x<tree[j].size-1;x++){if(tree[j].p[x][RSIZE+2]==SSIZE){q = q + tree[j].p[x][0] - tree[j].p[x+1][0];}}

//COALESCE
w=0;w1=0;
for(x=0;x<tree[total-1].size-1;x++){
for(l=1;l<RSIZE+2;l++){tree[total-1].p[x][l]=0;}
flag1=0;flag2=0;   
for(y=w;y<tree[i].size-1;y++){
if(tree[i].p[y][0]<=tree[total-1].p[x][0] && tree[i].p[y+1][0]>=tree[total-1].p[x+1][0]){ 
flag1 = tree[i].p[y][RSIZE+2];if(flag1>0){memcpy(&tree[total-1].p[x][1],&tree[i].p[y][1],4*(RSIZE+1));}
w=y;break;}}
   
for(y1=w1;y1<tree[j].size-1;y1++){
if(tree[j].p[y1][0]<=tree[total-1].p[x][0] && tree[j].p[y1+1][0]>=tree[total-1].p[x+1][0]){
flag2 = tree[j].p[y1][RSIZE+2];if(flag1==0){if(flag2>0){memcpy(&tree[total-1].p[x][1],&tree[j].p[y1][1],4*(RSIZE+1));}}
else{if(flag2>0){for(l=1;l<RSIZE+2;l++){tree[total-1].p[x][l] = tree[total-1].p[x][l] | tree[j].p[y1][l];}}}
w1=y1;break;}}
tree[total-1].p[x][RSIZE+2] = flag1 + flag2;
if(tree[total-1].p[x][RSIZE+2]==SSIZE){q = q + tree[total-1].p[x+1][0] - tree[total-1].p[x][0];}
}

//TRIM THE ENDS IF THEY ARE ALL-ANCESTRAL OR NON-ANCETSRAL
v2=tree[total-1].size-2;w2=0;
while(v2>=0){
if(tree[total-1].p[v2][RSIZE+2]>=1 && tree[total-1].p[v2][RSIZE+2]<=SSIZE-1){break;}
else{--v2;}}

while(w2<tree[total-1].size-1){
if(tree[total-1].p[w2][RSIZE+2]>=1 && tree[total-1].p[w2][RSIZE+2]<=SSIZE-1){break;} 
else{++w2;}}
 
if(v2<w2){memmove(&current[min-1],&current[min],4*(count-min+1));
delete tree[total-1].p;--count;}

//MERGE ALL-ANCESTRAL AND NON-ANCESTRAL INTERVALS 
else{ 
w=0;x=w2; 
while(x<=v2){
l=tree[total-1].p[x][RSIZE+2];
if(l==0||l==SSIZE){
tree[total-1].p[w][0] = tree[total-1].p[x][0];
while(l==0||l==SSIZE){x=x+1;l=tree[total-1].p[x][RSIZE+2];}
for(y=1;y<RSIZE+3;y++){tree[total-1].p[w][y]=0;}
++w;}

if(l>=1 && l<=SSIZE-1){
memmove(&tree[total-1].p[w][0],&tree[total-1].p[x][0],4*(RSIZE+3));
++w;++x;}
}
tree[total-1].p[w][0] = tree[total-1].p[v2+1][0];tree[total-1].size=w+1;
}


//PUT MUTATIONS IN THE SPECIFIED SITES
time1 = time - tree[i].time;k=1;
for(x=0;x<sites;x++){ 
if(position[x] > tree[i].p[0][0] && position[x]<=tree[i].p[tree[i].size-1][0]){
while(k<tree[i].size){if(position[x]<=tree[i].p[k][0]){break;}++k;}
maf1 = int(minimum(tree[i].p[k-1][RSIZE+2],SSIZE-tree[i].p[k-1][RSIZE+2]));
maf2 = tree[i].p[k-1][RSIZE+2];

if(maf1>=lfs[x] && maf1<=ufs[x]){
for(l=1;l<RSIZE+2;l++){
w = tree[i].p[k-1][l];
for(m=0;m<31;m++){if(w%2==1){break;}w = w/2;}
if(w%2==1){break;}}
m = 31*(l-1) + m;
if(prob1[x][m][maf2]==1){prob[x][int(prob2[x][m][maf2])] = prob[x][int(prob2[x][m][maf2])]+time1;}
if(prob1[x][m][maf2]==0){memcpy(&all[x][ccount[x]][0],&tree[i].p[k-1][1],4*(RSIZE+1));
prob1[x][m][maf2]=1;prob2[x][m][maf2]=ccount[x];
prob[x][ccount[x]]=time1;++ccount[x];}
}}}


time2 = time - tree[j].time;k=1;
for(x=0;x<sites;x++){
if(position[x] > tree[j].p[0][0] && position[x]<=tree[j].p[tree[j].size-1][0]){
while(k<tree[j].size){if(position[x]<=tree[j].p[k][0]){break;}++k;}
maf1 = int(minimum(tree[j].p[k-1][RSIZE+2],SSIZE-tree[j].p[k-1][RSIZE+2]));
maf2 = tree[j].p[k-1][RSIZE+2];

if(maf1>=lfs[x] && maf1<=ufs[x]){
for(l=1;l<RSIZE+2;l++){
w = tree[j].p[k-1][l];
for(m=0;m<31;m++){if(w%2==1){break;}w=w/2;}
if(w%2==1){break;}}
m = 31*(l-1) + m;

if(prob1[x][m][maf2]==1){prob[x][int(prob2[x][m][maf2])]= prob[x][int(prob2[x][m][maf2])]+time2;}
if(prob1[x][m][maf2]==0){memcpy(&all[x][ccount[x]][0],&tree[j].p[k-1][1],4*(RSIZE+1));
prob1[x][m][maf2]=1;prob2[x][m][maf2]=ccount[x];
prob[x][ccount[x]]=time2;++ccount[x];}
}}}


//FREE MEMORY SPACE
delete tree[i].p;delete tree[j].p;


//CHECK FOR MRCA
if(q<0){"Fatal Error\n";return(0);}
stop = stop + q;
if(stop==LEN){
for(x=0;x<count;x++){delete tree[current[x]-1].p;}
break;} 
}

  
//IF THE NEXT EVENT IS RECOMBINATION
else{

//CHOOSE A LINEAGE AT RANDOM TO RECOMBINE AND UPDATE TIME
j = choose(1,count+1);i = current[j-1] - 1;time = time + r;

double g,z;int start,end,tract;
g = drand();f = GAMMA/(RHO+GAMMA);brk=2*LEN;

//IF NEXT EVENT IS GENE-CONVERSION
if(g<=f){
tract=1;z=drand();while(z>1.00/L){z=drand();++tract;}

//IF GENE-CONVERSION IS UNIFORM
if(ghot==0){start = -tract + int((LEN+tract-1)*drand())+1;end = start + tract;}

//IF GENE-CONVERSION IS NON-UNIFORM
else{
z=drand();
if(z<=gfrac){start = gstart-tract + int((tract+gwidth-1)*drand()) + 1;end = start + tract;}
else{
start = -tract + int((LEN+tract-1)*drand()) + 1;
while(start > gstart-tract && start <=gstart+gwidth){start = -tract + int((LEN+tract-1)*drand()) + 1;}
end = start + tract;}
}

if(start>tree[i].p[0][0]  && start<tree[i].p[tree[i].size-1][0] && end>=tree[i].p[tree[i].size-1][0]){brk = start;}
if(start<=tree[i].p[0][0] && end>tree[i].p[0][0] && end<tree[i].p[tree[i].size-1][0]){brk = end;}
}


//IF NEXT EVENT IS CROSSING-OVER
else{

//IF CROSSING-OVER IS UNIFORM
if(rhot==0){brk = 1 + int((LEN-1)*drand());}

//IF CROSSING-OVER IS NON-UNIFORM
else{
z=drand();
if(z<=rfrac){ brk = rstart + int(rwidth*drand());}
else{brk = 1 + int((LEN-1)*drand());while(brk>=rstart && brk<rstart+rwidth){brk = 1 + int((LEN-1)*drand());}}
}}


//PERFORM CROSSING-OVER IF THE BREAKPOINT LIES WITHIN ANCESTRAL REGION
if(brk>tree[i].p[0][0] && brk<tree[i].p[tree[i].size-1][0]){

//PUT MUTATIONS IN THE SPECIFIED SITES
time1 = time - tree[i].time;k=1;
for(x=0;x<sites;x++){
if(position[x] > tree[i].p[tree[i].size-1][0]){break;}
if(position[x] > tree[i].p[0][0]){
while(k<tree[i].size){if(position[x]<=tree[i].p[k][0]){break;}++k;}
maf1 = int(minimum(tree[i].p[k-1][RSIZE+2],SSIZE-tree[i].p[k-1][RSIZE+2]));
maf2 = tree[i].p[k-1][RSIZE+2];

if(maf1>=lfs[x] && maf1<=ufs[x]){
for(l=1;l<RSIZE+2;l++){
w = tree[i].p[k-1][l];
for(m=0;m<31;m++){if(w%2==1){break;} w = w/2;}
if(w%2==1){break;}}
m=31*(l-1) +m;
if(prob1[x][m][maf2]==1){prob[x][int(prob2[x][m][maf2])]=prob[x][int(prob2[x][m][maf2])]+time1;}
if(prob1[x][m][maf2]==0){memcpy(&all[x][ccount[x]][0],&tree[i].p[k-1][1],4*(RSIZE+1));
prob1[x][m][maf2]=1;prob2[x][m][maf2]=ccount[x];
prob[x][ccount[x]]=time1;++ccount[x];}
}}
}

++count;++total;tree[total-1].time = time;

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


tree[total-1].size = tree[i].size - w;tree[total-1].p = new int[tree[i].size - w][RSIZE+3];
if(tree[i].size-w==0){cout<<"Error\n";return(0);}if(v==-1){cout<<"Error\n";return(0);}
 
//SPLIT THE LINEAGE TO 2 NEW ARRAYS
memcpy(&tree[total-1].p[0][0],&tree[i].p[w][0],4*(RSIZE+3)*tree[total-1].size);
if(w==x-1){tree[total-1].p[0][0] = brk;}
tree[i].size = v+2;tree[i].time = time;if(v==x-1){tree[i].p[x][0] = brk;}
}


else{
w=x;v=x-1;
while(v>=0){
if(tree[i].p[v][RSIZE+2]>=1 && tree[i].p[v][RSIZE+2]<=SSIZE-1){break;}
else{--v;}}

while(w<tree[i].size-1){
if(tree[i].p[w][RSIZE+2]>=1 && tree[i].p[w][RSIZE+2]<=SSIZE-1){break;}
else{++w;}}

tree[total-1].size = tree[i].size - w;tree[total-1].p = new int[tree[i].size - w][RSIZE+3];
if(tree[i].size-w==0){cout<<"Error\n";return(0);}if(v==-1){cout<<"Fatal Error6\n";return(0);}
memcpy(&tree[total-1].p[0][0],&tree[i].p[w][0],4*(RSIZE+3)*tree[total-1].size);
tree[i].size = v+2;tree[i].time = time;}
continue;}


//PERFORM GENE-CONVERSION IF BREAKPOINTS COVER ANCESTRAL REGION
if(g<=f && start>tree[i].p[0][0] && end<tree[i].p[tree[i].size-1][0]){

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
++count;++total;++total;
tree[total-1].time = time;tree[total-2].time=time;

//PUT MUTATIONS IN THE SPECIFIED SITES
time1 = time - tree[i].time;k=1;
for(x=0;x<sites;x++){
if(position[x] > tree[i].p[tree[i].size-1][0]){break;}
if(position[x] > tree[i].p[0][0]){
while(k<tree[i].size){if(position[x]<=tree[i].p[k][0]){break;}++k;}
maf1 = int(minimum(tree[i].p[k-1][RSIZE+2],SSIZE-tree[i].p[k-1][RSIZE+2]));
maf2 = tree[i].p[k-1][RSIZE+2];

if(maf1>=lfs[x] && maf1<=ufs[x]){
for(l=1;l<RSIZE+2;l++){
w = tree[i].p[k-1][l];
for(m=0;m<31;m++){if(w%2==1){break;} w = w/2;}
if(w%2==1){break;}}
m=31*(l-1) + m;
if(prob1[x][m][maf2]==1){prob[x][int(prob2[x][m][maf2])]=prob[x][int(prob2[x][m][maf2])]+time1;}
if(prob1[x][m][maf2]==0){memcpy(&all[x][ccount[x]][0],&tree[i].p[k-1][1],4*(RSIZE+1));
prob1[x][m][maf2]=1;prob2[x][m][maf2]=ccount[x];
prob[x][ccount[x]]=time1;++ccount[x];}
}}
}

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
delete tree[i].p;continue;
}
}}
}



int choice[sites],indic[SSIZE],flag,Hap,distance,cont;
double e1=0,e2=0,e3=0,e4=0,e5=0,e6=0,e7=0,e8=0,e9=0,e10=0;double e11=0,e12=0,e13=0,e14=0,e15=0,e16=0,e17=0,e18=0,e19=0,f1=0,f2=0,f3=0,f4=0,D,t1,t2,t3;
for(i=0;i<SSIZE;i++){indic[i]=0;}
cont=0;for(i=0;i<sites;i++){choice[i]=0;if(ccount[i]==0){cont=1;}}
if(cont==1){continue;}


for(i=0;i<sites;i++){
for(j=1;j<ccount[i];j++){prob[i][j]=prob[i][j]+prob[i][j-1];}
a = prob[i][ccount[i]-1]*drand();
for(j=0;j<ccount[i];j++){if(a<=prob[i][j]){break;}}
choice[i]=j;}


for(i=0;i<sites;i++){
k=0;
for(l=0;l<RSIZE;l++){w = all[i][choice[i]][l];for(m=0;m<31;m++){matrix[i][k]=w%2;++k;w = w/2;}}
w = all[i][choice[i]][l];for(m=0;m<=(SSIZE-1)%31;m++){matrix[i][k]=w%2;++k;w = w/2;} 
}
 

cout<<"Total SNPs: "<<sites<<"\n";
cout<<"Positions: ";for(i=0;i<sites;i++){cout<<position[i]<<" ";}cout<<"\n";
for(i=0;i<SSIZE;i++){for(j=0;j<sites;j++){cout<<matrix[j][i];}cout<<"\n";}
++n;}
}}



}






 
 



 



 


 














