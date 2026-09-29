// Exact two-strand targeted k-mer counts over every contig in an input FASTA.
// Two-bit uint128 encoding is collision-free for k<=63; hash collisions are
// resolved by unordered_map equality. Counts are genomic positions, not strands.
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cstdint>
using U=unsigned __int128;
struct Hash {size_t operator()(U x)const {uint64_t a=x,b=x>>64; return a^(b+0x9e3779b97f4a7c15ULL+(a<<6)+(a>>2));}};
struct Entry {std::string seq;uint64_t count=0;std::vector<std::string>loc;};
struct Table {int k;U f=0,r=0,mask;int valid=0;std::unordered_map<U,Entry,Hash> map;Table(int n):k(n),mask((U(1)<<(2*n))-1){} };
int val(char c){switch(c){case 'A':case 'a':return 0;case 'C':case 'c':return 1;case 'G':case 'g':return 2;case 'T':case 't':return 3;default:return -1;}}
U key(const std::string&s){U a=0,b=0;for(int i=0;i<(int)s.size();++i){a=(a<<2)|val(s[i]);b|=U(3-val(s[i]))<<(2*i);}return std::min(a,b);}
int main(int argc,char**argv){
 if(argc!=4){std::cerr<<"usage: scan probes.tsv reference.fa output.tsv\n";return 2;}
 std::vector<Table>ts;ts.emplace_back(31);ts.emplace_back(51);
 std::ifstream pin(argv[1]);std::string line;getline(pin,line);
 while(getline(pin,line)){if(line.empty())continue;bool ok=false;for(auto&t:ts)if(line.size()==(size_t)t.k){t.map.emplace(key(line),Entry{line,0,{}});ok=true;}if(!ok)return 3;}
 std::ifstream in(argv[2]);if(!in)return 4;std::string chrom;uint64_t pos=0,total=0,contigs=0;
 while(getline(in,line)){
  if(line.empty())continue;
  if(line[0]=='>'){chrom=line.substr(1,line.find_first_of(" \t")-1);pos=0;contigs++;for(auto&t:ts)t.valid=0;continue;}
  for(char c:line){if(c=='\r')continue;++pos;++total;int v=val(c);
   for(auto&t:ts){if(v<0){t.valid=0;t.f=0;t.r=0;continue;}
    t.f=((t.f<<2)|v)&t.mask;t.r=(t.r>>2)|(U(3-v)<<(2*(t.k-1)));
    if(t.valid<t.k)++t.valid;if(t.valid<t.k)continue;
    auto it=t.map.find(std::min(t.f,t.r));if(it!=t.map.end()){auto&e=it->second;++e.count;if(e.loc.size()<10)e.loc.push_back(chrom+":"+std::to_string(pos-t.k+1));}
   }
  }
 }
 std::ofstream out(argv[3]);out<<"canonical\tk\treference_occurrences\tfirst_positions\n";
 for(auto&t:ts)for(auto&kv:t.map){auto&e=kv.second;out<<e.seq<<'\t'<<t.k<<'\t'<<e.count<<'\t';for(size_t i=0;i<e.loc.size();++i){if(i)out<<',';out<<e.loc[i];}out<<'\n';}
 std::cerr<<"scanned_contigs="<<contigs<<" reference_bases="<<total<<" probes="<<ts[0].map.size()+ts[1].map.size()<<'\n';
}
