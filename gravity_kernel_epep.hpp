#include<pikg_vector.hpp>
#include<cmath>
#include<limits>
#include<chrono>

#include <pikg_avx512.hpp>
struct CalcForceLongEPEP{
PIKG::F32 eps2;
CalcForceLongEPEP(){}
CalcForceLongEPEP(PIKG::F32 eps2):eps2(eps2){}
void initialize(PIKG::F32 eps2_){
eps2 = eps2_;
}
int kernel_id = 0;
void operator()(const EPI_t* __restrict__ epi,const int ni,const EPJ_t* __restrict__ epj,const int nj,Force_t* __restrict__ force,const int kernel_select = 1){
static_assert(sizeof(EPI_t) == 48,"check consistency of EPI member variable definition between PIKG source and original source");
static_assert(sizeof(EPJ_t) == 112,"check consistency of EPJ member variable definition between PIKG source and original source");
static_assert(sizeof(Force_t) == 32,"check consistency of FORCE member variable definition between PIKG source and original source");
if(kernel_select>=0) kernel_id = kernel_select;
if(kernel_id == 0){
std::cout << "ni: " << ni << " nj:" << nj << std::endl;
Force_t* force_tmp = new Force_t[ni];
std::chrono::system_clock::time_point  start, end;
double min_time = std::numeric_limits<double>::max();
{ // test Kernel_I16_J1
for(int i=0;i<ni;i++) force_tmp[i] = force[i];
start = std::chrono::system_clock::now();
Kernel_I16_J1(epi,ni,epj,nj,force_tmp);
end = std::chrono::system_clock::now();
double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
std::cerr << "kerel 1: " << elapsed << " ns" << std::endl;
if(min_time > elapsed){
min_time = elapsed;
kernel_id = 1;
}
}
{ // test Kernel_I1_J16
for(int i=0;i<ni;i++) force_tmp[i] = force[i];
start = std::chrono::system_clock::now();
Kernel_I1_J16(epi,ni,epj,nj,force_tmp);
end = std::chrono::system_clock::now();
double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
std::cerr << "kerel 2: " << elapsed << " ns" << std::endl;
if(min_time > elapsed){
min_time = elapsed;
kernel_id = 2;
}
}
delete[] force_tmp;
} // if(kernel_id == 0)
if(kernel_id == 1) Kernel_I16_J1(epi,ni,epj,nj,force);
if(kernel_id == 2) Kernel_I1_J16(epi,ni,epj,nj,force);
} // operator() definition 
void Kernel_I16_J1(const EPI_t* __restrict__ epi,const PIKG::S32 ni,const EPJ_t* __restrict__ epj,const PIKG::S32 nj,Force_t* __restrict__ force){
PIKG::S32 i;
PIKG::S32 j;
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_x[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_y[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_z[ni];
PIKG::F32  __attribute__ ((aligned(64))) rsearchiloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) routiloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_x[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_y[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_z[nj];
PIKG::F32  __attribute__ ((aligned(64))) mjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) rsearchjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) routjloc_tmp[nj];
for(i = 0;i < ni;++i){
xiloc_tmp_x[i] = (epi[i].pos.x-epi[0].pos.x);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_y[i] = (epi[i].pos.y-epi[0].pos.y);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_z[i] = (epi[i].pos.z-epi[0].pos.z);
} // loop of i
for(i = 0;i < ni;++i){
rsearchiloc_tmp[i] = epi[i].r_search;
} // loop of i
for(i = 0;i < ni;++i){
routiloc_tmp[i] = epi[i].r_out;
} // loop of i
for(j = 0;j < nj;++j){
xjloc_tmp_x[j] = (epj[j].pos.x-epi[0].pos.x);
} // loop of j
for(j = 0;j < nj;++j){
xjloc_tmp_y[j] = (epj[j].pos.y-epi[0].pos.y);
} // loop of j
for(j = 0;j < nj;++j){
xjloc_tmp_z[j] = (epj[j].pos.z-epi[0].pos.z);
} // loop of j
for(j = 0;j < nj;++j){
mjloc_tmp[j] = epj[j].mass;
} // loop of j
for(j = 0;j < nj;++j){
rsearchjloc_tmp[j] = epj[j].r_search;
} // loop of j
for(j = 0;j < nj;++j){
routjloc_tmp[j] = epj[j].r_out;
} // loop of j
for(i = 0;i < (ni/16)*16;i += 16){
__m512i idi;

alignas(32) int32_t index_gather_load0[16] = {0,12,24,36,48,60,72,84,96,108,120,132,144,156,168,180};
__m512i vindex_gather_load0 = _mm512_load_epi32(index_gather_load0);
idi = _mm512_i32gather_epi32(vindex_gather_load0,((int*)&epi[i+0].id),4);
__m512 routiloc;

routiloc = _mm512_load_ps(((float*)&routiloc_tmp[i+0]));
__m512 rsearchiloc;

rsearchiloc = _mm512_load_ps(((float*)&rsearchiloc_tmp[i+0]));
__m512x3 xiloc;

xiloc.v0 = _mm512_load_ps(((float*)&xiloc_tmp_x[i+0]));
xiloc.v1 = _mm512_load_ps(((float*)&xiloc_tmp_y[i+0]));
xiloc.v2 = _mm512_load_ps(((float*)&xiloc_tmp_z[i+0]));
__m512x3 af;

af.v0 = _mm512_set1_ps(0.0f);
af.v1 = _mm512_set1_ps(0.0f);
af.v2 = _mm512_set1_ps(0.0f);
__m512i idngb;

idngb = _mm512_set1_epi32(std::numeric_limits<int32_t>::lowest());
__m512i nngb;

nngb = _mm512_set1_epi32(0);
__m512 pf;

pf = _mm512_set1_ps(0.0f);
for(j = 0;j < (nj/1)*1;++j){
__m512i idj;

idj = _mm512_set1_epi32(epj[j].id);

__m512 mjloc;

mjloc = _mm512_set1_ps(mjloc_tmp[j+0]);
__m512 routjloc;

routjloc = _mm512_set1_ps(routjloc_tmp[j+0]);
__m512 rsearchjloc;

rsearchjloc = _mm512_set1_ps(rsearchjloc_tmp[j+0]);
__m512x3 xjloc;

xjloc.v0 = _mm512_set1_ps(xjloc_tmp_x[j+0]);
xjloc.v1 = _mm512_set1_ps(xjloc_tmp_y[j+0]);
xjloc.v2 = _mm512_set1_ps(xjloc_tmp_z[j+0]);
__m512 rout;

__m512 rsearch;

__m512 rout2;

__m512 __fkg_tmp2;

__m512 rsearch2;

__m512x3 rij;

__m512 __fkg_tmp4;

__m512 __fkg_tmp3;

__m512 r2_real;

__m512 r2;

__m512i __fkg_tmp1;

__m512i __fkg_tmp0;

__m512 r_inv;

__m512 __fkg_tmp5;

__m512 tmp;

__m512 __fkg_tmp6;

__m512 r2_inv;

__m512 mr_inv;

__m512 mr3_inv;

rout = max(routjloc,routiloc);
rsearch = max(rsearchjloc,rsearchiloc);
rout2 = _mm512_mul_ps(rout,rout);
__fkg_tmp2 = _mm512_mul_ps(rsearch,rsearch);
rsearch2 = _mm512_mul_ps(__fkg_tmp2,_mm512_set1_ps(1.1025f));
rij.v0 = _mm512_sub_ps(xjloc.v0,xiloc.v0);
rij.v1 = _mm512_sub_ps(xjloc.v1,xiloc.v1);
rij.v2 = _mm512_sub_ps(xjloc.v2,xiloc.v2);
__fkg_tmp4 = _mm512_fmadd_ps(rij.v0,rij.v0,_mm512_set1_ps(eps2));
__fkg_tmp3 = _mm512_fmadd_ps(rij.v1,rij.v1,__fkg_tmp4);
r2_real = _mm512_fmadd_ps(rij.v2,rij.v2,__fkg_tmp3);
r2 = r2_real;
r2 = max(r2_real,rout2);
{
__mmask16 pg1;
__mmask16 pg0;
pg1 = _mm512_cmp_ps_mask(r2_real,rsearch2,_CMP_LT_OQ);
pg0 = pg1;

__fkg_tmp1 = _mm512_add_epi32(nngb,_mm512_set1_epi32(1));
{
__mmask16 pg3;
__mmask16 pg2;
pg3 = _mm512_cmp_epi32_mask(idi,idj,_MM_CMPINT_NE);
pg2 = pg3;
pg3 = _kand_mask16(pg3,pg1);

__fkg_tmp0 = max(idngb,idj);
idngb = _mm512_mask_blend_epi32(pg3,idngb,__fkg_tmp0);;
}

nngb = _mm512_mask_blend_epi32(pg1,nngb,__fkg_tmp1);;
}

r_inv = rsqrt(r2);
__fkg_tmp5 = _mm512_mul_ps(r_inv,r_inv);
tmp = _mm512_fnmadd_ps(r2,__fkg_tmp5,_mm512_set1_ps(3.0f));
__fkg_tmp6 = _mm512_mul_ps(tmp,_mm512_set1_ps(0.5f));
r_inv = _mm512_mul_ps(r_inv,__fkg_tmp6);
r2_inv = _mm512_mul_ps(r_inv,r_inv);
mr_inv = _mm512_mul_ps(mjloc,r_inv);
mr3_inv = _mm512_mul_ps(r2_inv,mr_inv);
af.v0 = _mm512_fmadd_ps(mr3_inv,rij.v0,af.v0);
af.v1 = _mm512_fmadd_ps(mr3_inv,rij.v1,af.v1);
af.v2 = _mm512_fmadd_ps(mr3_inv,rij.v2,af.v2);
pf = _mm512_sub_ps(pf,mr_inv);
} // loop of j

{
__m512 __fkg_tmp_accum;
alignas(32) int32_t index_gather_load1[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_gather_load1 = _mm512_load_epi32(index_gather_load1);
__fkg_tmp_accum = _mm512_i32gather_ps(vindex_gather_load1,((float*)&force[i+0].acc.x),4);
__fkg_tmp_accum = _mm512_add_ps(__fkg_tmp_accum,af.v0);
int32_t index_scatter_store0[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_scatter_store0 = _mm512_load_epi32(index_scatter_store0);
_mm512_i32scatter_ps(((float*)&force[i+0].acc.x),vindex_scatter_store0,__fkg_tmp_accum,4);
}

{
__m512 __fkg_tmp_accum;
alignas(32) int32_t index_gather_load2[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_gather_load2 = _mm512_load_epi32(index_gather_load2);
__fkg_tmp_accum = _mm512_i32gather_ps(vindex_gather_load2,((float*)&force[i+0].acc.y),4);
__fkg_tmp_accum = _mm512_add_ps(__fkg_tmp_accum,af.v1);
int32_t index_scatter_store1[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_scatter_store1 = _mm512_load_epi32(index_scatter_store1);
_mm512_i32scatter_ps(((float*)&force[i+0].acc.y),vindex_scatter_store1,__fkg_tmp_accum,4);
}

{
__m512 __fkg_tmp_accum;
alignas(32) int32_t index_gather_load3[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_gather_load3 = _mm512_load_epi32(index_gather_load3);
__fkg_tmp_accum = _mm512_i32gather_ps(vindex_gather_load3,((float*)&force[i+0].acc.z),4);
__fkg_tmp_accum = _mm512_add_ps(__fkg_tmp_accum,af.v2);
int32_t index_scatter_store2[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_scatter_store2 = _mm512_load_epi32(index_scatter_store2);
_mm512_i32scatter_ps(((float*)&force[i+0].acc.z),vindex_scatter_store2,__fkg_tmp_accum,4);
}

{
__m512i __fkg_tmp_accum;
alignas(32) int32_t index_gather_load4[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_gather_load4 = _mm512_load_epi32(index_gather_load4);
__fkg_tmp_accum = _mm512_i32gather_epi32(vindex_gather_load4,((int*)&force[i+0].id_neighbor),4);
__fkg_tmp_accum = max(__fkg_tmp_accum,idngb);
int32_t index_scatter_store3[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_scatter_store3 = _mm512_load_epi32(index_scatter_store3);
_mm512_i32scatter_epi32(((int*)&force[i+0].id_neighbor),vindex_scatter_store3,__fkg_tmp_accum,4);
}

{
__m512i __fkg_tmp_accum;
alignas(32) int32_t index_gather_load5[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_gather_load5 = _mm512_load_epi32(index_gather_load5);
__fkg_tmp_accum = _mm512_i32gather_epi32(vindex_gather_load5,((int*)&force[i+0].neighbor),4);
__fkg_tmp_accum = _mm512_add_epi32(__fkg_tmp_accum,nngb);
int32_t index_scatter_store4[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_scatter_store4 = _mm512_load_epi32(index_scatter_store4);
_mm512_i32scatter_epi32(((int*)&force[i+0].neighbor),vindex_scatter_store4,__fkg_tmp_accum,4);
}

{
__m512 __fkg_tmp_accum;
alignas(32) int32_t index_gather_load6[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_gather_load6 = _mm512_load_epi32(index_gather_load6);
__fkg_tmp_accum = _mm512_i32gather_ps(vindex_gather_load6,((float*)&force[i+0].phi),4);
__fkg_tmp_accum = _mm512_add_ps(__fkg_tmp_accum,pf);
int32_t index_scatter_store5[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_scatter_store5 = _mm512_load_epi32(index_scatter_store5);
_mm512_i32scatter_ps(((float*)&force[i+0].phi),vindex_scatter_store5,__fkg_tmp_accum,4);
}

} // loop of i
{ // tail loop of reference 
for(;i < ni;++i){
PIKG::S32 idi;

idi = epi[i+0].id;
PIKG::F32 routiloc;

routiloc = routiloc_tmp[i+0];
PIKG::F32 rsearchiloc;

rsearchiloc = rsearchiloc_tmp[i+0];
PIKG::F32vec xiloc;

xiloc.x = xiloc_tmp_x[i+0];
xiloc.y = xiloc_tmp_y[i+0];
xiloc.z = xiloc_tmp_z[i+0];
PIKG::F32vec af;

af.x = 0.0f;
af.y = 0.0f;
af.z = 0.0f;
PIKG::S32 idngb;

idngb = std::numeric_limits<int32_t>::lowest();
PIKG::S32 nngb;

nngb = 0;
PIKG::F32 pf;

pf = 0.0f;
for(j = 0;j < nj;++j){
PIKG::S32 idj;

idj = epj[j].id;
PIKG::F32 mjloc;

mjloc = mjloc_tmp[j+0];
PIKG::F32 routjloc;

routjloc = routjloc_tmp[j+0];
PIKG::F32 rsearchjloc;

rsearchjloc = rsearchjloc_tmp[j+0];
PIKG::F32vec xjloc;

xjloc.x = xjloc_tmp_x[j+0];
xjloc.y = xjloc_tmp_y[j+0];
xjloc.z = xjloc_tmp_z[j+0];
PIKG::F32 rout;

PIKG::F32 rsearch;

PIKG::F32 rout2;

PIKG::F32 __fkg_tmp2;

PIKG::F32 rsearch2;

PIKG::F32vec rij;

PIKG::F32 __fkg_tmp4;

PIKG::F32 __fkg_tmp3;

PIKG::F32 r2_real;

PIKG::F32 r2;

PIKG::S32 __fkg_tmp1;

PIKG::S32 __fkg_tmp0;

PIKG::F32 r_inv;

PIKG::F32 __fkg_tmp5;

PIKG::F32 tmp;

PIKG::F32 __fkg_tmp6;

PIKG::F32 r2_inv;

PIKG::F32 mr_inv;

PIKG::F32 mr3_inv;

rout = max(routjloc,routiloc);
rsearch = max(rsearchjloc,rsearchiloc);
rout2 = (rout*rout);
__fkg_tmp2 = (rsearch*rsearch);
rsearch2 = (__fkg_tmp2*1.1025f);
rij.x = (xjloc.x-xiloc.x);
rij.y = (xjloc.y-xiloc.y);
rij.z = (xjloc.z-xiloc.z);
__fkg_tmp4 = (rij.x*rij.x+eps2);
__fkg_tmp3 = (rij.y*rij.y+__fkg_tmp4);
r2_real = (rij.z*rij.z+__fkg_tmp3);
r2 = r2_real;
r2 = max(r2_real,rout2);
if((r2_real<rsearch2)){
__fkg_tmp1 = (nngb+1);
if((idi!=idj)){
__fkg_tmp0 = max(idngb,idj);
idngb = __fkg_tmp0;
}
nngb = __fkg_tmp1;
}
r_inv = rsqrt(r2);
__fkg_tmp5 = (r_inv*r_inv);
tmp = (3.0f - r2*__fkg_tmp5);
__fkg_tmp6 = (tmp*0.5f);
r_inv = (r_inv*__fkg_tmp6);
r2_inv = (r_inv*r_inv);
mr_inv = (mjloc*r_inv);
mr3_inv = (r2_inv*mr_inv);
af.x = (mr3_inv*rij.x+af.x);
af.y = (mr3_inv*rij.y+af.y);
af.z = (mr3_inv*rij.z+af.z);
pf = (pf-mr_inv);
} // loop of j

force[i+0].acc.x = (force[i+0].acc.x+af.x);
force[i+0].acc.y = (force[i+0].acc.y+af.y);
force[i+0].acc.z = (force[i+0].acc.z+af.z);
force[i+0].id_neighbor = max(idngb,force[i+0].id_neighbor);
force[i+0].neighbor = (force[i+0].neighbor+nngb);
force[i+0].phi = (force[i+0].phi+pf);
} // loop of i
} // end loop of reference 
} // Kernel_I16_J1 definition 
void Kernel_I1_J16(const EPI_t* __restrict__ epi,const PIKG::S32 ni,const EPJ_t* __restrict__ epj,const PIKG::S32 nj,Force_t* __restrict__ force){
PIKG::S32 i;
PIKG::S32 j;
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_x[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_y[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_z[ni];
PIKG::F32  __attribute__ ((aligned(64))) rsearchiloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) routiloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_x[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_y[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_z[nj];
PIKG::F32  __attribute__ ((aligned(64))) mjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) rsearchjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) routjloc_tmp[nj];
for(i = 0;i < ni;++i){
xiloc_tmp_x[i] = (epi[i].pos.x-epi[0].pos.x);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_y[i] = (epi[i].pos.y-epi[0].pos.y);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_z[i] = (epi[i].pos.z-epi[0].pos.z);
} // loop of i
for(i = 0;i < ni;++i){
rsearchiloc_tmp[i] = epi[i].r_search;
} // loop of i
for(i = 0;i < ni;++i){
routiloc_tmp[i] = epi[i].r_out;
} // loop of i
for(j = 0;j < nj;++j){
xjloc_tmp_x[j] = (epj[j].pos.x-epi[0].pos.x);
} // loop of j
for(j = 0;j < nj;++j){
xjloc_tmp_y[j] = (epj[j].pos.y-epi[0].pos.y);
} // loop of j
for(j = 0;j < nj;++j){
xjloc_tmp_z[j] = (epj[j].pos.z-epi[0].pos.z);
} // loop of j
for(j = 0;j < nj;++j){
mjloc_tmp[j] = epj[j].mass;
} // loop of j
for(j = 0;j < nj;++j){
rsearchjloc_tmp[j] = epj[j].r_search;
} // loop of j
for(j = 0;j < nj;++j){
routjloc_tmp[j] = epj[j].r_out;
} // loop of j
for(i = 0;i < (ni/1)*1;++i){
__m512i idi;

idi = _mm512_set1_epi32(epi[i+0].id);

__m512 routiloc;

routiloc = _mm512_set1_ps(routiloc_tmp[i+0]);
__m512 rsearchiloc;

rsearchiloc = _mm512_set1_ps(rsearchiloc_tmp[i+0]);
__m512x3 xiloc;

xiloc.v0 = _mm512_set1_ps(xiloc_tmp_x[i+0]);
xiloc.v1 = _mm512_set1_ps(xiloc_tmp_y[i+0]);
xiloc.v2 = _mm512_set1_ps(xiloc_tmp_z[i+0]);
__m512x3 af;

af.v0 = _mm512_set1_ps(0.0f);
af.v1 = _mm512_set1_ps(0.0f);
af.v2 = _mm512_set1_ps(0.0f);
__m512i idngb;

idngb = _mm512_set1_epi32(std::numeric_limits<int32_t>::lowest());
__m512i nngb;

nngb = _mm512_set1_epi32(0);
__m512 pf;

pf = _mm512_set1_ps(0.0f);
for(j = 0;j < (nj/16)*16;j += 16){
__m512i idj;

alignas(32) int32_t index_gather_load7[16] = {0,28,56,84,112,140,168,196,224,252,280,308,336,364,392,420};
__m512i vindex_gather_load7 = _mm512_load_epi32(index_gather_load7);
idj = _mm512_i32gather_epi32(vindex_gather_load7,((int*)&epj[j].id),4);
__m512 mjloc;

mjloc = _mm512_load_ps(((float*)&mjloc_tmp[j+0]));
__m512 routjloc;

routjloc = _mm512_load_ps(((float*)&routjloc_tmp[j+0]));
__m512 rsearchjloc;

rsearchjloc = _mm512_load_ps(((float*)&rsearchjloc_tmp[j+0]));
__m512x3 xjloc;

xjloc.v0 = _mm512_load_ps(((float*)&xjloc_tmp_x[j+0]));
xjloc.v1 = _mm512_load_ps(((float*)&xjloc_tmp_y[j+0]));
xjloc.v2 = _mm512_load_ps(((float*)&xjloc_tmp_z[j+0]));
__m512 rout;

__m512 rsearch;

__m512 rout2;

__m512 __fkg_tmp2;

__m512 rsearch2;

__m512x3 rij;

__m512 __fkg_tmp4;

__m512 __fkg_tmp3;

__m512 r2_real;

__m512 r2;

__m512i __fkg_tmp1;

__m512i __fkg_tmp0;

__m512 r_inv;

__m512 __fkg_tmp5;

__m512 tmp;

__m512 __fkg_tmp6;

__m512 r2_inv;

__m512 mr_inv;

__m512 mr3_inv;

rout = max(routjloc,routiloc);
rsearch = max(rsearchjloc,rsearchiloc);
rout2 = _mm512_mul_ps(rout,rout);
__fkg_tmp2 = _mm512_mul_ps(rsearch,rsearch);
rsearch2 = _mm512_mul_ps(__fkg_tmp2,_mm512_set1_ps(1.1025f));
rij.v0 = _mm512_sub_ps(xjloc.v0,xiloc.v0);
rij.v1 = _mm512_sub_ps(xjloc.v1,xiloc.v1);
rij.v2 = _mm512_sub_ps(xjloc.v2,xiloc.v2);
__fkg_tmp4 = _mm512_fmadd_ps(rij.v0,rij.v0,_mm512_set1_ps(eps2));
__fkg_tmp3 = _mm512_fmadd_ps(rij.v1,rij.v1,__fkg_tmp4);
r2_real = _mm512_fmadd_ps(rij.v2,rij.v2,__fkg_tmp3);
r2 = r2_real;
r2 = max(r2_real,rout2);
{
__mmask16 pg1;
__mmask16 pg0;
pg1 = _mm512_cmp_ps_mask(r2_real,rsearch2,_CMP_LT_OQ);
pg0 = pg1;

__fkg_tmp1 = _mm512_add_epi32(nngb,_mm512_set1_epi32(1));
{
__mmask16 pg3;
__mmask16 pg2;
pg3 = _mm512_cmp_epi32_mask(idi,idj,_MM_CMPINT_NE);
pg2 = pg3;
pg3 = _kand_mask16(pg3,pg1);

__fkg_tmp0 = max(idngb,idj);
idngb = _mm512_mask_blend_epi32(pg3,idngb,__fkg_tmp0);;
}

nngb = _mm512_mask_blend_epi32(pg1,nngb,__fkg_tmp1);;
}

r_inv = rsqrt(r2);
__fkg_tmp5 = _mm512_mul_ps(r_inv,r_inv);
tmp = _mm512_fnmadd_ps(r2,__fkg_tmp5,_mm512_set1_ps(3.0f));
__fkg_tmp6 = _mm512_mul_ps(tmp,_mm512_set1_ps(0.5f));
r_inv = _mm512_mul_ps(r_inv,__fkg_tmp6);
r2_inv = _mm512_mul_ps(r_inv,r_inv);
mr_inv = _mm512_mul_ps(mjloc,r_inv);
mr3_inv = _mm512_mul_ps(r2_inv,mr_inv);
af.v0 = _mm512_fmadd_ps(mr3_inv,rij.v0,af.v0);
af.v1 = _mm512_fmadd_ps(mr3_inv,rij.v1,af.v1);
af.v2 = _mm512_fmadd_ps(mr3_inv,rij.v2,af.v2);
pf = _mm512_sub_ps(pf,mr_inv);
} // loop of j

if(j<nj){ // tail j loop
__m512x3 __fkg_tmp7;

__fkg_tmp7.v0 = af.v0;
__fkg_tmp7.v1 = af.v1;
__fkg_tmp7.v2 = af.v2;
__m512i __fkg_tmp8;

__fkg_tmp8 = idngb;
__m512i __fkg_tmp9;

__fkg_tmp9 = nngb;
__m512 __fkg_tmp10;

__fkg_tmp10 = pf;
for(;j < nj;++j){
__m512i idj;

idj = _mm512_set1_epi32(epj[j].id);

__m512 mjloc;

mjloc = _mm512_set1_ps(mjloc_tmp[j+0]);
__m512 routjloc;

routjloc = _mm512_set1_ps(routjloc_tmp[j+0]);
__m512 rsearchjloc;

rsearchjloc = _mm512_set1_ps(rsearchjloc_tmp[j+0]);
__m512x3 xjloc;

xjloc.v0 = _mm512_set1_ps(xjloc_tmp_x[j+0]);
xjloc.v1 = _mm512_set1_ps(xjloc_tmp_y[j+0]);
xjloc.v2 = _mm512_set1_ps(xjloc_tmp_z[j+0]);
__m512 rout;

__m512 rsearch;

__m512 rout2;

__m512 __fkg_tmp2;

__m512 rsearch2;

__m512x3 rij;

__m512 __fkg_tmp4;

__m512 __fkg_tmp3;

__m512 r2_real;

__m512 r2;

__m512i __fkg_tmp1;

__m512i __fkg_tmp0;

__m512 r_inv;

__m512 __fkg_tmp5;

__m512 tmp;

__m512 __fkg_tmp6;

__m512 r2_inv;

__m512 mr_inv;

__m512 mr3_inv;

rout = max(routjloc,routiloc);
rsearch = max(rsearchjloc,rsearchiloc);
rout2 = _mm512_mul_ps(rout,rout);
__fkg_tmp2 = _mm512_mul_ps(rsearch,rsearch);
rsearch2 = _mm512_mul_ps(__fkg_tmp2,_mm512_set1_ps(1.1025f));
rij.v0 = _mm512_sub_ps(xjloc.v0,xiloc.v0);
rij.v1 = _mm512_sub_ps(xjloc.v1,xiloc.v1);
rij.v2 = _mm512_sub_ps(xjloc.v2,xiloc.v2);
__fkg_tmp4 = _mm512_fmadd_ps(rij.v0,rij.v0,_mm512_set1_ps(eps2));
__fkg_tmp3 = _mm512_fmadd_ps(rij.v1,rij.v1,__fkg_tmp4);
r2_real = _mm512_fmadd_ps(rij.v2,rij.v2,__fkg_tmp3);
r2 = r2_real;
r2 = max(r2_real,rout2);
{
__mmask16 pg1;
__mmask16 pg0;
pg1 = _mm512_cmp_ps_mask(r2_real,rsearch2,_CMP_LT_OQ);
pg0 = pg1;

__fkg_tmp1 = _mm512_add_epi32(nngb,_mm512_set1_epi32(1));
{
__mmask16 pg3;
__mmask16 pg2;
pg3 = _mm512_cmp_epi32_mask(idi,idj,_MM_CMPINT_NE);
pg2 = pg3;
pg3 = _kand_mask16(pg3,pg1);

__fkg_tmp0 = max(idngb,idj);
idngb = _mm512_mask_blend_epi32(pg3,idngb,__fkg_tmp0);;
}

nngb = _mm512_mask_blend_epi32(pg1,nngb,__fkg_tmp1);;
}

r_inv = rsqrt(r2);
__fkg_tmp5 = _mm512_mul_ps(r_inv,r_inv);
tmp = _mm512_fnmadd_ps(r2,__fkg_tmp5,_mm512_set1_ps(3.0f));
__fkg_tmp6 = _mm512_mul_ps(tmp,_mm512_set1_ps(0.5f));
r_inv = _mm512_mul_ps(r_inv,__fkg_tmp6);
r2_inv = _mm512_mul_ps(r_inv,r_inv);
mr_inv = _mm512_mul_ps(mjloc,r_inv);
mr3_inv = _mm512_mul_ps(r2_inv,mr_inv);
af.v0 = _mm512_fmadd_ps(mr3_inv,rij.v0,af.v0);
af.v1 = _mm512_fmadd_ps(mr3_inv,rij.v1,af.v1);
af.v2 = _mm512_fmadd_ps(mr3_inv,rij.v2,af.v2);
pf = _mm512_sub_ps(pf,mr_inv);
} // loop of j
af.v0 = _mm512_mask_blend_ps(_cvtu32_mask16(0b00000001),__fkg_tmp7.v0,af.v0);
af.v1 = _mm512_mask_blend_ps(_cvtu32_mask16(0b00000001),__fkg_tmp7.v1,af.v1);
af.v2 = _mm512_mask_blend_ps(_cvtu32_mask16(0b00000001),__fkg_tmp7.v2,af.v2);
idngb = _mm512_mask_blend_epi32(_cvtu32_mask16(0b00000001),__fkg_tmp8,idngb);
nngb = _mm512_mask_blend_epi32(_cvtu32_mask16(0b00000001),__fkg_tmp9,nngb);
pf = _mm512_mask_blend_ps(_cvtu32_mask16(0b00000001),__fkg_tmp10,pf);
} // if of j tail loop

((float*)&force[i+0].acc.x)[0] += _mm512_reduce_add_ps(af.v0);

((float*)&force[i+0].acc.y)[0] += _mm512_reduce_add_ps(af.v1);

((float*)&force[i+0].acc.z)[0] += _mm512_reduce_add_ps(af.v2);

((int*)&force[i+0].id_neighbor)[0] = _mm512_reduce_max_epi32(idngb);

((int*)&force[i+0].neighbor)[0] += _mm512_reduce_add_epi32(nngb);

((float*)&force[i+0].phi)[0] += _mm512_reduce_add_ps(pf);

} // loop of i
{ // tail loop of reference 
for(;i < ni;++i){
PIKG::S32 idi;

idi = epi[i+0].id;
PIKG::F32 routiloc;

routiloc = routiloc_tmp[i+0];
PIKG::F32 rsearchiloc;

rsearchiloc = rsearchiloc_tmp[i+0];
PIKG::F32vec xiloc;

xiloc.x = xiloc_tmp_x[i+0];
xiloc.y = xiloc_tmp_y[i+0];
xiloc.z = xiloc_tmp_z[i+0];
PIKG::F32vec af;

af.x = 0.0f;
af.y = 0.0f;
af.z = 0.0f;
PIKG::S32 idngb;

idngb = std::numeric_limits<int32_t>::lowest();
PIKG::S32 nngb;

nngb = 0;
PIKG::F32 pf;

pf = 0.0f;
for(j = 0;j < nj;++j){
PIKG::S32 idj;

idj = epj[j].id;
PIKG::F32 mjloc;

mjloc = mjloc_tmp[j+0];
PIKG::F32 routjloc;

routjloc = routjloc_tmp[j+0];
PIKG::F32 rsearchjloc;

rsearchjloc = rsearchjloc_tmp[j+0];
PIKG::F32vec xjloc;

xjloc.x = xjloc_tmp_x[j+0];
xjloc.y = xjloc_tmp_y[j+0];
xjloc.z = xjloc_tmp_z[j+0];
PIKG::F32 rout;

PIKG::F32 rsearch;

PIKG::F32 rout2;

PIKG::F32 __fkg_tmp2;

PIKG::F32 rsearch2;

PIKG::F32vec rij;

PIKG::F32 __fkg_tmp4;

PIKG::F32 __fkg_tmp3;

PIKG::F32 r2_real;

PIKG::F32 r2;

PIKG::S32 __fkg_tmp1;

PIKG::S32 __fkg_tmp0;

PIKG::F32 r_inv;

PIKG::F32 __fkg_tmp5;

PIKG::F32 tmp;

PIKG::F32 __fkg_tmp6;

PIKG::F32 r2_inv;

PIKG::F32 mr_inv;

PIKG::F32 mr3_inv;

rout = max(routjloc,routiloc);
rsearch = max(rsearchjloc,rsearchiloc);
rout2 = (rout*rout);
__fkg_tmp2 = (rsearch*rsearch);
rsearch2 = (__fkg_tmp2*1.1025f);
rij.x = (xjloc.x-xiloc.x);
rij.y = (xjloc.y-xiloc.y);
rij.z = (xjloc.z-xiloc.z);
__fkg_tmp4 = (rij.x*rij.x+eps2);
__fkg_tmp3 = (rij.y*rij.y+__fkg_tmp4);
r2_real = (rij.z*rij.z+__fkg_tmp3);
r2 = r2_real;
r2 = max(r2_real,rout2);
if((r2_real<rsearch2)){
__fkg_tmp1 = (nngb+1);
if((idi!=idj)){
__fkg_tmp0 = max(idngb,idj);
idngb = __fkg_tmp0;
}
nngb = __fkg_tmp1;
}
r_inv = rsqrt(r2);
__fkg_tmp5 = (r_inv*r_inv);
tmp = (3.0f - r2*__fkg_tmp5);
__fkg_tmp6 = (tmp*0.5f);
r_inv = (r_inv*__fkg_tmp6);
r2_inv = (r_inv*r_inv);
mr_inv = (mjloc*r_inv);
mr3_inv = (r2_inv*mr_inv);
af.x = (mr3_inv*rij.x+af.x);
af.y = (mr3_inv*rij.y+af.y);
af.z = (mr3_inv*rij.z+af.z);
pf = (pf-mr_inv);
} // loop of j

force[i+0].acc.x = (force[i+0].acc.x+af.x);
force[i+0].acc.y = (force[i+0].acc.y+af.y);
force[i+0].acc.z = (force[i+0].acc.z+af.z);
force[i+0].id_neighbor = max(idngb,force[i+0].id_neighbor);
force[i+0].neighbor = (force[i+0].neighbor+nngb);
force[i+0].phi = (force[i+0].phi+pf);
} // loop of i
} // end loop of reference 
} // Kernel_I1_J16 definition 
PIKG::F64 rsqrt(PIKG::F64 op){ return 1.0/std::sqrt(op); }
PIKG::F64 sqrt(PIKG::F64 op){ return std::sqrt(op); }
PIKG::F64 inv(PIKG::F64 op){ return 1.0/op; }
PIKG::F64 max(PIKG::F64 a,PIKG::F64 b){ return std::max(a,b);}
PIKG::F64 min(PIKG::F64 a,PIKG::F64 b){ return std::min(a,b);}
PIKG::F32 rsqrt(PIKG::F32 op){ return 1.f/std::sqrt(op); }
PIKG::F32 sqrt(PIKG::F32 op){ return std::sqrt(op); }
PIKG::F32 inv(PIKG::F32 op){ return 1.f/op; }
PIKG::S64 max(PIKG::S64 a,PIKG::S64 b){ return std::max(a,b);}
PIKG::S64 min(PIKG::S64 a,PIKG::S64 b){ return std::min(a,b);}
PIKG::S32 max(PIKG::S32 a,PIKG::S32 b){ return std::max(a,b);}
PIKG::S32 min(PIKG::S32 a,PIKG::S32 b){ return std::min(a,b);}
PIKG::F64 table(PIKG::F64 tab[],PIKG::S64 i){ return tab[i]; }
PIKG::F32 table(PIKG::F32 tab[],PIKG::S32 i){ return tab[i]; }
PIKG::F64 to_float(PIKG::U64 op){return (PIKG::F64)op;}
PIKG::F32 to_float(PIKG::U32 op){return (PIKG::F32)op;}
PIKG::F64 to_float(PIKG::S64 op){return (PIKG::F64)op;}
PIKG::F32 to_float(PIKG::S32 op){return (PIKG::F32)op;}
PIKG::S64   to_int(PIKG::F64 op){return (PIKG::S64)op;}
PIKG::S32   to_int(PIKG::F32 op){return (PIKG::S32)op;}
PIKG::U64  to_uint(PIKG::F64 op){return (PIKG::U64)op;}
PIKG::U32  to_uint(PIKG::F32 op){return (PIKG::U32)op;}
template<typename T> PIKG::F64 to_f64(const T& op){return (PIKG::F64)op;}
template<typename T> PIKG::F32 to_f32(const T& op){return (PIKG::F32)op;}
template<typename T> PIKG::S64 to_s64(const T& op){return (PIKG::S64)op;}
template<typename T> PIKG::S32 to_s32(const T& op){return (PIKG::S32)op;}
template<typename T> PIKG::U64 to_u64(const T& op){return (PIKG::U64)op;}
template<typename T> PIKG::U32 to_u32(const T& op){return (PIKG::U32)op;}
__m512 rsqrt(__m512 op){
return _mm512_rsqrt14_ps(op);
}
__m512 sqrt(__m512 op){ return _mm512_sqrt_ps(op); }
__m512 inv(__m512 op){
__m512 x1 = _mm512_rcp14_ps(op);
__m512 x2 = _mm512_fnmadd_ps(op,x1,_mm512_set1_ps(2.f));
x2 = _mm512_mul_ps(x2,x1);
__m512 ret = _mm512_fnmadd_ps(op,x2,_mm512_set1_ps(2.f));
ret = _mm512_mul_ps(ret,x2);
return ret;
}
__m512d rsqrt(__m512d op){
__m512d rinv = _mm512_rsqrt14_pd(op);
__m512d h = _mm512_mul_pd(op,rinv);
h = _mm512_fnmadd_pd(h,rinv,_mm512_set1_pd(1.0));
__m512d poly = _mm512_fmadd_pd(h,_mm512_set1_pd(0.375),_mm512_set1_pd(0.5));
poly = _mm512_mul_pd(poly,h);
return _mm512_fmadd_pd(rinv,poly,rinv);
}
__m512d max(__m512d a,__m512d b){ return _mm512_max_pd(a,b);}
__m512d min(__m512d a,__m512d b){ return _mm512_min_pd(a,b);}
__m512  max(__m512  a,__m512  b){ return _mm512_max_ps(a,b);}
__m512  min(__m512  a,__m512  b){ return _mm512_min_ps(a,b);}
__m512i max(__m512i a,__m512i b){ return _mm512_max_epi32(a,b);}
__m512i min(__m512i a,__m512i b){ return _mm512_min_epi32(a,b);}
__m512d table(__m512d tab,__m512i index){ return _mm512_permutexvar_pd(index,tab);}
__m512  table(__m512  tab,__m512i index){ return _mm512_permutexvar_ps(index,tab);}
__m512d to_double(__m512i op){ return _mm512_cvt_roundepi64_pd(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}
__m512  to_float(__m512i op){ return _mm512_cvt_roundepi32_ps(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}
__m512i  to_long(__m512d op){ return _mm512_cvt_roundpd_epi64(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}
__m512i  to_int(__m512  op){ return _mm512_cvt_roundps_epi32(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}
__m512i  to_ulong(__m512d op){ return _mm512_cvt_roundpd_epu64(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}
__m512i  to_uint(__m512  op){ return _mm512_cvt_roundps_epu32(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}
};// kernel functor definition 
