#include<pikg_vector.hpp>
#include<cmath>
#include<limits>
#include<chrono>

#include <pikg_avx512.hpp>
struct CalcForceLongEPSP{
PIKG::F32 eps2;
CalcForceLongEPSP(){}
CalcForceLongEPSP(PIKG::F32 eps2):eps2(eps2){}
void initialize(PIKG::F32 eps2_){
eps2 = eps2_;
}
int kernel_id = 0;
void operator()(const EPI_t* __restrict__ epi,const int ni,const SPJ_t* __restrict__ epj,const int nj,Force_t* __restrict__ force,const int kernel_select = 1){
static_assert(sizeof(EPI_t) == 48,"check consistency of EPI member variable definition between PIKG source and original source");
static_assert(sizeof(SPJ_t) == 80,"check consistency of EPJ member variable definition between PIKG source and original source");
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
void Kernel_I16_J1(const EPI_t* __restrict__ epi,const PIKG::S32 ni,const SPJ_t* __restrict__ epj,const PIKG::S32 nj,Force_t* __restrict__ force){
PIKG::S32 i;
PIKG::S32 j;
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_x[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_y[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_z[ni];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_x[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_y[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_z[nj];
PIKG::F32  __attribute__ ((aligned(64))) mjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_xxloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_yyloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_zzloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_xyloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_yzloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_zxloc_tmp[nj];
for(i = 0;i < ni;++i){
xiloc_tmp_x[i] = (epi[i].pos.x-epi[0].pos.x);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_y[i] = (epi[i].pos.y-epi[0].pos.y);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_z[i] = (epi[i].pos.z-epi[0].pos.z);
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
qj_xxloc_tmp[j] = epj[j].quad.xx;
} // loop of j
for(j = 0;j < nj;++j){
qj_yyloc_tmp[j] = epj[j].quad.yy;
} // loop of j
for(j = 0;j < nj;++j){
qj_zzloc_tmp[j] = epj[j].quad.zz;
} // loop of j
for(j = 0;j < nj;++j){
qj_xyloc_tmp[j] = epj[j].quad.xy;
} // loop of j
for(j = 0;j < nj;++j){
qj_yzloc_tmp[j] = epj[j].quad.yz;
} // loop of j
for(j = 0;j < nj;++j){
qj_zxloc_tmp[j] = epj[j].quad.xz;
} // loop of j
for(i = 0;i < (ni/16)*16;i += 16){
__m512x3 xiloc;

xiloc.v0 = _mm512_load_ps(((float*)&xiloc_tmp_x[i+0]));
xiloc.v1 = _mm512_load_ps(((float*)&xiloc_tmp_y[i+0]));
xiloc.v2 = _mm512_load_ps(((float*)&xiloc_tmp_z[i+0]));
__m512x3 af;

af.v0 = _mm512_set1_ps(0.0f);
af.v1 = _mm512_set1_ps(0.0f);
af.v2 = _mm512_set1_ps(0.0f);
__m512 pf;

pf = _mm512_set1_ps(0.0f);
for(j = 0;j < (nj/1)*1;++j){
__m512 mjloc;

mjloc = _mm512_set1_ps(mjloc_tmp[j+0]);
__m512 qj_xxloc;

qj_xxloc = _mm512_set1_ps(qj_xxloc_tmp[j+0]);
__m512 qj_xyloc;

qj_xyloc = _mm512_set1_ps(qj_xyloc_tmp[j+0]);
__m512 qj_yyloc;

qj_yyloc = _mm512_set1_ps(qj_yyloc_tmp[j+0]);
__m512 qj_yzloc;

qj_yzloc = _mm512_set1_ps(qj_yzloc_tmp[j+0]);
__m512 qj_zxloc;

qj_zxloc = _mm512_set1_ps(qj_zxloc_tmp[j+0]);
__m512 qj_zzloc;

qj_zzloc = _mm512_set1_ps(qj_zzloc_tmp[j+0]);
__m512x3 xjloc;

xjloc.v0 = _mm512_set1_ps(xjloc_tmp_x[j+0]);
xjloc.v1 = _mm512_set1_ps(xjloc_tmp_y[j+0]);
xjloc.v2 = _mm512_set1_ps(xjloc_tmp_z[j+0]);
__m512x3 rij;

__m512 __fkg_tmp1;

__m512 __fkg_tmp0;

__m512 r2;

__m512 r_inv;

__m512 __fkg_tmp2;

__m512 tmp;

__m512 __fkg_tmp3;

__m512 r2_inv;

__m512 r3_inv;

__m512 r4_inv;

__m512 r5_inv;

__m512 __fkg_tmp4;

__m512 tr;

__m512 qxx;

__m512 qyy;

__m512 qzz;

__m512 qxy;

__m512 qyz;

__m512 qzx;

__m512 __fkg_tmp5;

__m512 mtr;

__m512 __fkg_tmp7;

__m512 __fkg_tmp6;

__m512x3 qr;

__m512 __fkg_tmp9;

__m512 __fkg_tmp8;

__m512 __fkg_tmp11;

__m512 __fkg_tmp10;

__m512 __fkg_tmp13;

__m512 __fkg_tmp12;

__m512 rqr;

__m512 rqr_r4_inv;

__m512 meff;

__m512 __fkg_tmp14;

__m512 meff3;

__m512 __fkg_tmp15;

__m512 __fkg_tmp16;

__m512 __fkg_tmp17;

rij.v0 = _mm512_sub_ps(xjloc.v0,xiloc.v0);
rij.v1 = _mm512_sub_ps(xjloc.v1,xiloc.v1);
rij.v2 = _mm512_sub_ps(xjloc.v2,xiloc.v2);
__fkg_tmp1 = _mm512_fmadd_ps(rij.v0,rij.v0,_mm512_set1_ps(eps2));
__fkg_tmp0 = _mm512_fmadd_ps(rij.v1,rij.v1,__fkg_tmp1);
r2 = _mm512_fmadd_ps(rij.v2,rij.v2,__fkg_tmp0);
r_inv = rsqrt(r2);
__fkg_tmp2 = _mm512_mul_ps(r_inv,r_inv);
tmp = _mm512_fnmadd_ps(r2,__fkg_tmp2,_mm512_set1_ps(3.0f));
__fkg_tmp3 = _mm512_mul_ps(tmp,_mm512_set1_ps(0.5f));
r_inv = _mm512_mul_ps(r_inv,__fkg_tmp3);
r2_inv = _mm512_mul_ps(r_inv,r_inv);
r3_inv = _mm512_mul_ps(r2_inv,r_inv);
r4_inv = _mm512_mul_ps(r2_inv,r2_inv);
r5_inv = _mm512_mul_ps(r2_inv,r3_inv);
__fkg_tmp4 = _mm512_add_ps(qj_xxloc,qj_yyloc);
tr = _mm512_add_ps(__fkg_tmp4,qj_zzloc);
qxx = _mm512_fmsub_ps(_mm512_set1_ps(3.0f),qj_xxloc,tr);
qyy = _mm512_fmsub_ps(_mm512_set1_ps(3.0f),qj_yyloc,tr);
qzz = _mm512_fmsub_ps(_mm512_set1_ps(3.0f),qj_zzloc,tr);
qxy = _mm512_mul_ps(_mm512_set1_ps(3.0f),qj_xyloc);
qyz = _mm512_mul_ps(_mm512_set1_ps(3.0f),qj_yzloc);
qzx = _mm512_mul_ps(_mm512_set1_ps(3.0f),qj_zxloc);
__fkg_tmp5 = _mm512_mul_ps(_mm512_set1_ps(eps2),tr);
mtr = _mm512_sub_ps(_mm512_set1_ps((PIKG::F32)0.0),__fkg_tmp5);
__fkg_tmp7 = _mm512_mul_ps(qxx,rij.v0);
__fkg_tmp6 = _mm512_fmadd_ps(qxy,rij.v1,__fkg_tmp7);
qr.v0 = _mm512_fmadd_ps(qzx,rij.v2,__fkg_tmp6);
__fkg_tmp9 = _mm512_mul_ps(qyy,rij.v1);
__fkg_tmp8 = _mm512_fmadd_ps(qyz,rij.v2,__fkg_tmp9);
qr.v1 = _mm512_fmadd_ps(qxy,rij.v0,__fkg_tmp8);
__fkg_tmp11 = _mm512_mul_ps(qzz,rij.v2);
__fkg_tmp10 = _mm512_fmadd_ps(qzx,rij.v0,__fkg_tmp11);
qr.v2 = _mm512_fmadd_ps(qyz,rij.v1,__fkg_tmp10);
__fkg_tmp13 = _mm512_fmadd_ps(qr.v0,rij.v0,mtr);
__fkg_tmp12 = _mm512_fmadd_ps(qr.v1,rij.v1,__fkg_tmp13);
rqr = _mm512_fmadd_ps(qr.v2,rij.v2,__fkg_tmp12);
rqr_r4_inv = _mm512_mul_ps(rqr,r4_inv);
meff = _mm512_fmadd_ps(_mm512_set1_ps(0.5f),rqr_r4_inv,mjloc);
__fkg_tmp14 = _mm512_fmadd_ps(_mm512_set1_ps(2.5f),rqr_r4_inv,mjloc);
meff3 = _mm512_mul_ps(__fkg_tmp14,r3_inv);
pf = _mm512_fnmadd_ps(meff,r_inv,pf);
__fkg_tmp15 = _mm512_fnmadd_ps(r5_inv,qr.v0,af.v0);
af.v0 = _mm512_fmadd_ps(meff3,rij.v0,__fkg_tmp15);
__fkg_tmp16 = _mm512_fnmadd_ps(r5_inv,qr.v1,af.v1);
af.v1 = _mm512_fmadd_ps(meff3,rij.v1,__fkg_tmp16);
__fkg_tmp17 = _mm512_fnmadd_ps(r5_inv,qr.v2,af.v2);
af.v2 = _mm512_fmadd_ps(meff3,rij.v2,__fkg_tmp17);
} // loop of j

{
__m512 __fkg_tmp_accum;
alignas(32) int32_t index_gather_load0[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_gather_load0 = _mm512_load_epi32(index_gather_load0);
__fkg_tmp_accum = _mm512_i32gather_ps(vindex_gather_load0,((float*)&force[i+0].acc.x),4);
__fkg_tmp_accum = _mm512_add_ps(__fkg_tmp_accum,af.v0);
int32_t index_scatter_store0[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_scatter_store0 = _mm512_load_epi32(index_scatter_store0);
_mm512_i32scatter_ps(((float*)&force[i+0].acc.x),vindex_scatter_store0,__fkg_tmp_accum,4);
}

{
__m512 __fkg_tmp_accum;
alignas(32) int32_t index_gather_load1[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_gather_load1 = _mm512_load_epi32(index_gather_load1);
__fkg_tmp_accum = _mm512_i32gather_ps(vindex_gather_load1,((float*)&force[i+0].acc.y),4);
__fkg_tmp_accum = _mm512_add_ps(__fkg_tmp_accum,af.v1);
int32_t index_scatter_store1[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_scatter_store1 = _mm512_load_epi32(index_scatter_store1);
_mm512_i32scatter_ps(((float*)&force[i+0].acc.y),vindex_scatter_store1,__fkg_tmp_accum,4);
}

{
__m512 __fkg_tmp_accum;
alignas(32) int32_t index_gather_load2[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_gather_load2 = _mm512_load_epi32(index_gather_load2);
__fkg_tmp_accum = _mm512_i32gather_ps(vindex_gather_load2,((float*)&force[i+0].acc.z),4);
__fkg_tmp_accum = _mm512_add_ps(__fkg_tmp_accum,af.v2);
int32_t index_scatter_store2[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_scatter_store2 = _mm512_load_epi32(index_scatter_store2);
_mm512_i32scatter_ps(((float*)&force[i+0].acc.z),vindex_scatter_store2,__fkg_tmp_accum,4);
}

{
__m512 __fkg_tmp_accum;
alignas(32) int32_t index_gather_load3[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_gather_load3 = _mm512_load_epi32(index_gather_load3);
__fkg_tmp_accum = _mm512_i32gather_ps(vindex_gather_load3,((float*)&force[i+0].phi),4);
__fkg_tmp_accum = _mm512_add_ps(__fkg_tmp_accum,pf);
int32_t index_scatter_store3[16] = {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120};
__m512i vindex_scatter_store3 = _mm512_load_epi32(index_scatter_store3);
_mm512_i32scatter_ps(((float*)&force[i+0].phi),vindex_scatter_store3,__fkg_tmp_accum,4);
}

} // loop of i
{ // tail loop of reference 
for(;i < ni;++i){
PIKG::F32vec xiloc;

xiloc.x = xiloc_tmp_x[i+0];
xiloc.y = xiloc_tmp_y[i+0];
xiloc.z = xiloc_tmp_z[i+0];
PIKG::F32vec af;

af.x = 0.0f;
af.y = 0.0f;
af.z = 0.0f;
PIKG::F32 pf;

pf = 0.0f;
for(j = 0;j < nj;++j){
PIKG::F32 mjloc;

mjloc = mjloc_tmp[j+0];
PIKG::F32 qj_xxloc;

qj_xxloc = qj_xxloc_tmp[j+0];
PIKG::F32 qj_xyloc;

qj_xyloc = qj_xyloc_tmp[j+0];
PIKG::F32 qj_yyloc;

qj_yyloc = qj_yyloc_tmp[j+0];
PIKG::F32 qj_yzloc;

qj_yzloc = qj_yzloc_tmp[j+0];
PIKG::F32 qj_zxloc;

qj_zxloc = qj_zxloc_tmp[j+0];
PIKG::F32 qj_zzloc;

qj_zzloc = qj_zzloc_tmp[j+0];
PIKG::F32vec xjloc;

xjloc.x = xjloc_tmp_x[j+0];
xjloc.y = xjloc_tmp_y[j+0];
xjloc.z = xjloc_tmp_z[j+0];
PIKG::F32vec rij;

PIKG::F32 __fkg_tmp1;

PIKG::F32 __fkg_tmp0;

PIKG::F32 r2;

PIKG::F32 r_inv;

PIKG::F32 __fkg_tmp2;

PIKG::F32 tmp;

PIKG::F32 __fkg_tmp3;

PIKG::F32 r2_inv;

PIKG::F32 r3_inv;

PIKG::F32 r4_inv;

PIKG::F32 r5_inv;

PIKG::F32 __fkg_tmp4;

PIKG::F32 tr;

PIKG::F32 qxx;

PIKG::F32 qyy;

PIKG::F32 qzz;

PIKG::F32 qxy;

PIKG::F32 qyz;

PIKG::F32 qzx;

PIKG::F32 __fkg_tmp5;

PIKG::F32 mtr;

PIKG::F32 __fkg_tmp7;

PIKG::F32 __fkg_tmp6;

PIKG::F32vec qr;

PIKG::F32 __fkg_tmp9;

PIKG::F32 __fkg_tmp8;

PIKG::F32 __fkg_tmp11;

PIKG::F32 __fkg_tmp10;

PIKG::F32 __fkg_tmp13;

PIKG::F32 __fkg_tmp12;

PIKG::F32 rqr;

PIKG::F32 rqr_r4_inv;

PIKG::F32 meff;

PIKG::F32 __fkg_tmp14;

PIKG::F32 meff3;

PIKG::F32 __fkg_tmp15;

PIKG::F32 __fkg_tmp16;

PIKG::F32 __fkg_tmp17;

rij.x = (xjloc.x-xiloc.x);
rij.y = (xjloc.y-xiloc.y);
rij.z = (xjloc.z-xiloc.z);
__fkg_tmp1 = (rij.x*rij.x+eps2);
__fkg_tmp0 = (rij.y*rij.y+__fkg_tmp1);
r2 = (rij.z*rij.z+__fkg_tmp0);
r_inv = rsqrt(r2);
__fkg_tmp2 = (r_inv*r_inv);
tmp = (3.0f - r2*__fkg_tmp2);
__fkg_tmp3 = (tmp*0.5f);
r_inv = (r_inv*__fkg_tmp3);
r2_inv = (r_inv*r_inv);
r3_inv = (r2_inv*r_inv);
r4_inv = (r2_inv*r2_inv);
r5_inv = (r2_inv*r3_inv);
__fkg_tmp4 = (qj_xxloc+qj_yyloc);
tr = (__fkg_tmp4+qj_zzloc);
qxx = (3.0f*qj_xxloc-tr);
qyy = (3.0f*qj_yyloc-tr);
qzz = (3.0f*qj_zzloc-tr);
qxy = (3.0f*qj_xyloc);
qyz = (3.0f*qj_yzloc);
qzx = (3.0f*qj_zxloc);
__fkg_tmp5 = (eps2*tr);
mtr = -(__fkg_tmp5);
__fkg_tmp7 = (qxx*rij.x);
__fkg_tmp6 = (qxy*rij.y+__fkg_tmp7);
qr.x = (qzx*rij.z+__fkg_tmp6);
__fkg_tmp9 = (qyy*rij.y);
__fkg_tmp8 = (qyz*rij.z+__fkg_tmp9);
qr.y = (qxy*rij.x+__fkg_tmp8);
__fkg_tmp11 = (qzz*rij.z);
__fkg_tmp10 = (qzx*rij.x+__fkg_tmp11);
qr.z = (qyz*rij.y+__fkg_tmp10);
__fkg_tmp13 = (qr.x*rij.x+mtr);
__fkg_tmp12 = (qr.y*rij.y+__fkg_tmp13);
rqr = (qr.z*rij.z+__fkg_tmp12);
rqr_r4_inv = (rqr*r4_inv);
meff = (0.5f*rqr_r4_inv+mjloc);
__fkg_tmp14 = (2.5f*rqr_r4_inv+mjloc);
meff3 = (__fkg_tmp14*r3_inv);
pf = (pf - meff*r_inv);
__fkg_tmp15 = (af.x - r5_inv*qr.x);
af.x = (meff3*rij.x+__fkg_tmp15);
__fkg_tmp16 = (af.y - r5_inv*qr.y);
af.y = (meff3*rij.y+__fkg_tmp16);
__fkg_tmp17 = (af.z - r5_inv*qr.z);
af.z = (meff3*rij.z+__fkg_tmp17);
} // loop of j

force[i+0].acc.x = (force[i+0].acc.x+af.x);
force[i+0].acc.y = (force[i+0].acc.y+af.y);
force[i+0].acc.z = (force[i+0].acc.z+af.z);
force[i+0].phi = (force[i+0].phi+pf);
} // loop of i
} // end loop of reference 
} // Kernel_I16_J1 definition 
void Kernel_I1_J16(const EPI_t* __restrict__ epi,const PIKG::S32 ni,const SPJ_t* __restrict__ epj,const PIKG::S32 nj,Force_t* __restrict__ force){
PIKG::S32 i;
PIKG::S32 j;
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_x[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_y[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_z[ni];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_x[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_y[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_z[nj];
PIKG::F32  __attribute__ ((aligned(64))) mjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_xxloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_yyloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_zzloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_xyloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_yzloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_zxloc_tmp[nj];
for(i = 0;i < ni;++i){
xiloc_tmp_x[i] = (epi[i].pos.x-epi[0].pos.x);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_y[i] = (epi[i].pos.y-epi[0].pos.y);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_z[i] = (epi[i].pos.z-epi[0].pos.z);
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
qj_xxloc_tmp[j] = epj[j].quad.xx;
} // loop of j
for(j = 0;j < nj;++j){
qj_yyloc_tmp[j] = epj[j].quad.yy;
} // loop of j
for(j = 0;j < nj;++j){
qj_zzloc_tmp[j] = epj[j].quad.zz;
} // loop of j
for(j = 0;j < nj;++j){
qj_xyloc_tmp[j] = epj[j].quad.xy;
} // loop of j
for(j = 0;j < nj;++j){
qj_yzloc_tmp[j] = epj[j].quad.yz;
} // loop of j
for(j = 0;j < nj;++j){
qj_zxloc_tmp[j] = epj[j].quad.xz;
} // loop of j
for(i = 0;i < (ni/1)*1;++i){
__m512x3 xiloc;

xiloc.v0 = _mm512_set1_ps(xiloc_tmp_x[i+0]);
xiloc.v1 = _mm512_set1_ps(xiloc_tmp_y[i+0]);
xiloc.v2 = _mm512_set1_ps(xiloc_tmp_z[i+0]);
__m512x3 af;

af.v0 = _mm512_set1_ps(0.0f);
af.v1 = _mm512_set1_ps(0.0f);
af.v2 = _mm512_set1_ps(0.0f);
__m512 pf;

pf = _mm512_set1_ps(0.0f);
for(j = 0;j < (nj/16)*16;j += 16){
__m512 mjloc;

mjloc = _mm512_load_ps(((float*)&mjloc_tmp[j+0]));
__m512 qj_xxloc;

qj_xxloc = _mm512_load_ps(((float*)&qj_xxloc_tmp[j+0]));
__m512 qj_xyloc;

qj_xyloc = _mm512_load_ps(((float*)&qj_xyloc_tmp[j+0]));
__m512 qj_yyloc;

qj_yyloc = _mm512_load_ps(((float*)&qj_yyloc_tmp[j+0]));
__m512 qj_yzloc;

qj_yzloc = _mm512_load_ps(((float*)&qj_yzloc_tmp[j+0]));
__m512 qj_zxloc;

qj_zxloc = _mm512_load_ps(((float*)&qj_zxloc_tmp[j+0]));
__m512 qj_zzloc;

qj_zzloc = _mm512_load_ps(((float*)&qj_zzloc_tmp[j+0]));
__m512x3 xjloc;

xjloc.v0 = _mm512_load_ps(((float*)&xjloc_tmp_x[j+0]));
xjloc.v1 = _mm512_load_ps(((float*)&xjloc_tmp_y[j+0]));
xjloc.v2 = _mm512_load_ps(((float*)&xjloc_tmp_z[j+0]));
__m512x3 rij;

__m512 __fkg_tmp1;

__m512 __fkg_tmp0;

__m512 r2;

__m512 r_inv;

__m512 __fkg_tmp2;

__m512 tmp;

__m512 __fkg_tmp3;

__m512 r2_inv;

__m512 r3_inv;

__m512 r4_inv;

__m512 r5_inv;

__m512 __fkg_tmp4;

__m512 tr;

__m512 qxx;

__m512 qyy;

__m512 qzz;

__m512 qxy;

__m512 qyz;

__m512 qzx;

__m512 __fkg_tmp5;

__m512 mtr;

__m512 __fkg_tmp7;

__m512 __fkg_tmp6;

__m512x3 qr;

__m512 __fkg_tmp9;

__m512 __fkg_tmp8;

__m512 __fkg_tmp11;

__m512 __fkg_tmp10;

__m512 __fkg_tmp13;

__m512 __fkg_tmp12;

__m512 rqr;

__m512 rqr_r4_inv;

__m512 meff;

__m512 __fkg_tmp14;

__m512 meff3;

__m512 __fkg_tmp15;

__m512 __fkg_tmp16;

__m512 __fkg_tmp17;

rij.v0 = _mm512_sub_ps(xjloc.v0,xiloc.v0);
rij.v1 = _mm512_sub_ps(xjloc.v1,xiloc.v1);
rij.v2 = _mm512_sub_ps(xjloc.v2,xiloc.v2);
__fkg_tmp1 = _mm512_fmadd_ps(rij.v0,rij.v0,_mm512_set1_ps(eps2));
__fkg_tmp0 = _mm512_fmadd_ps(rij.v1,rij.v1,__fkg_tmp1);
r2 = _mm512_fmadd_ps(rij.v2,rij.v2,__fkg_tmp0);
r_inv = rsqrt(r2);
__fkg_tmp2 = _mm512_mul_ps(r_inv,r_inv);
tmp = _mm512_fnmadd_ps(r2,__fkg_tmp2,_mm512_set1_ps(3.0f));
__fkg_tmp3 = _mm512_mul_ps(tmp,_mm512_set1_ps(0.5f));
r_inv = _mm512_mul_ps(r_inv,__fkg_tmp3);
r2_inv = _mm512_mul_ps(r_inv,r_inv);
r3_inv = _mm512_mul_ps(r2_inv,r_inv);
r4_inv = _mm512_mul_ps(r2_inv,r2_inv);
r5_inv = _mm512_mul_ps(r2_inv,r3_inv);
__fkg_tmp4 = _mm512_add_ps(qj_xxloc,qj_yyloc);
tr = _mm512_add_ps(__fkg_tmp4,qj_zzloc);
qxx = _mm512_fmsub_ps(_mm512_set1_ps(3.0f),qj_xxloc,tr);
qyy = _mm512_fmsub_ps(_mm512_set1_ps(3.0f),qj_yyloc,tr);
qzz = _mm512_fmsub_ps(_mm512_set1_ps(3.0f),qj_zzloc,tr);
qxy = _mm512_mul_ps(_mm512_set1_ps(3.0f),qj_xyloc);
qyz = _mm512_mul_ps(_mm512_set1_ps(3.0f),qj_yzloc);
qzx = _mm512_mul_ps(_mm512_set1_ps(3.0f),qj_zxloc);
__fkg_tmp5 = _mm512_mul_ps(_mm512_set1_ps(eps2),tr);
mtr = _mm512_sub_ps(_mm512_set1_ps((PIKG::F32)0.0),__fkg_tmp5);
__fkg_tmp7 = _mm512_mul_ps(qxx,rij.v0);
__fkg_tmp6 = _mm512_fmadd_ps(qxy,rij.v1,__fkg_tmp7);
qr.v0 = _mm512_fmadd_ps(qzx,rij.v2,__fkg_tmp6);
__fkg_tmp9 = _mm512_mul_ps(qyy,rij.v1);
__fkg_tmp8 = _mm512_fmadd_ps(qyz,rij.v2,__fkg_tmp9);
qr.v1 = _mm512_fmadd_ps(qxy,rij.v0,__fkg_tmp8);
__fkg_tmp11 = _mm512_mul_ps(qzz,rij.v2);
__fkg_tmp10 = _mm512_fmadd_ps(qzx,rij.v0,__fkg_tmp11);
qr.v2 = _mm512_fmadd_ps(qyz,rij.v1,__fkg_tmp10);
__fkg_tmp13 = _mm512_fmadd_ps(qr.v0,rij.v0,mtr);
__fkg_tmp12 = _mm512_fmadd_ps(qr.v1,rij.v1,__fkg_tmp13);
rqr = _mm512_fmadd_ps(qr.v2,rij.v2,__fkg_tmp12);
rqr_r4_inv = _mm512_mul_ps(rqr,r4_inv);
meff = _mm512_fmadd_ps(_mm512_set1_ps(0.5f),rqr_r4_inv,mjloc);
__fkg_tmp14 = _mm512_fmadd_ps(_mm512_set1_ps(2.5f),rqr_r4_inv,mjloc);
meff3 = _mm512_mul_ps(__fkg_tmp14,r3_inv);
pf = _mm512_fnmadd_ps(meff,r_inv,pf);
__fkg_tmp15 = _mm512_fnmadd_ps(r5_inv,qr.v0,af.v0);
af.v0 = _mm512_fmadd_ps(meff3,rij.v0,__fkg_tmp15);
__fkg_tmp16 = _mm512_fnmadd_ps(r5_inv,qr.v1,af.v1);
af.v1 = _mm512_fmadd_ps(meff3,rij.v1,__fkg_tmp16);
__fkg_tmp17 = _mm512_fnmadd_ps(r5_inv,qr.v2,af.v2);
af.v2 = _mm512_fmadd_ps(meff3,rij.v2,__fkg_tmp17);
} // loop of j

if(j<nj){ // tail j loop
__m512x3 __fkg_tmp18;

__fkg_tmp18.v0 = af.v0;
__fkg_tmp18.v1 = af.v1;
__fkg_tmp18.v2 = af.v2;
__m512 __fkg_tmp19;

__fkg_tmp19 = pf;
for(;j < nj;++j){
__m512 mjloc;

mjloc = _mm512_set1_ps(mjloc_tmp[j+0]);
__m512 qj_xxloc;

qj_xxloc = _mm512_set1_ps(qj_xxloc_tmp[j+0]);
__m512 qj_xyloc;

qj_xyloc = _mm512_set1_ps(qj_xyloc_tmp[j+0]);
__m512 qj_yyloc;

qj_yyloc = _mm512_set1_ps(qj_yyloc_tmp[j+0]);
__m512 qj_yzloc;

qj_yzloc = _mm512_set1_ps(qj_yzloc_tmp[j+0]);
__m512 qj_zxloc;

qj_zxloc = _mm512_set1_ps(qj_zxloc_tmp[j+0]);
__m512 qj_zzloc;

qj_zzloc = _mm512_set1_ps(qj_zzloc_tmp[j+0]);
__m512x3 xjloc;

xjloc.v0 = _mm512_set1_ps(xjloc_tmp_x[j+0]);
xjloc.v1 = _mm512_set1_ps(xjloc_tmp_y[j+0]);
xjloc.v2 = _mm512_set1_ps(xjloc_tmp_z[j+0]);
__m512x3 rij;

__m512 __fkg_tmp1;

__m512 __fkg_tmp0;

__m512 r2;

__m512 r_inv;

__m512 __fkg_tmp2;

__m512 tmp;

__m512 __fkg_tmp3;

__m512 r2_inv;

__m512 r3_inv;

__m512 r4_inv;

__m512 r5_inv;

__m512 __fkg_tmp4;

__m512 tr;

__m512 qxx;

__m512 qyy;

__m512 qzz;

__m512 qxy;

__m512 qyz;

__m512 qzx;

__m512 __fkg_tmp5;

__m512 mtr;

__m512 __fkg_tmp7;

__m512 __fkg_tmp6;

__m512x3 qr;

__m512 __fkg_tmp9;

__m512 __fkg_tmp8;

__m512 __fkg_tmp11;

__m512 __fkg_tmp10;

__m512 __fkg_tmp13;

__m512 __fkg_tmp12;

__m512 rqr;

__m512 rqr_r4_inv;

__m512 meff;

__m512 __fkg_tmp14;

__m512 meff3;

__m512 __fkg_tmp15;

__m512 __fkg_tmp16;

__m512 __fkg_tmp17;

rij.v0 = _mm512_sub_ps(xjloc.v0,xiloc.v0);
rij.v1 = _mm512_sub_ps(xjloc.v1,xiloc.v1);
rij.v2 = _mm512_sub_ps(xjloc.v2,xiloc.v2);
__fkg_tmp1 = _mm512_fmadd_ps(rij.v0,rij.v0,_mm512_set1_ps(eps2));
__fkg_tmp0 = _mm512_fmadd_ps(rij.v1,rij.v1,__fkg_tmp1);
r2 = _mm512_fmadd_ps(rij.v2,rij.v2,__fkg_tmp0);
r_inv = rsqrt(r2);
__fkg_tmp2 = _mm512_mul_ps(r_inv,r_inv);
tmp = _mm512_fnmadd_ps(r2,__fkg_tmp2,_mm512_set1_ps(3.0f));
__fkg_tmp3 = _mm512_mul_ps(tmp,_mm512_set1_ps(0.5f));
r_inv = _mm512_mul_ps(r_inv,__fkg_tmp3);
r2_inv = _mm512_mul_ps(r_inv,r_inv);
r3_inv = _mm512_mul_ps(r2_inv,r_inv);
r4_inv = _mm512_mul_ps(r2_inv,r2_inv);
r5_inv = _mm512_mul_ps(r2_inv,r3_inv);
__fkg_tmp4 = _mm512_add_ps(qj_xxloc,qj_yyloc);
tr = _mm512_add_ps(__fkg_tmp4,qj_zzloc);
qxx = _mm512_fmsub_ps(_mm512_set1_ps(3.0f),qj_xxloc,tr);
qyy = _mm512_fmsub_ps(_mm512_set1_ps(3.0f),qj_yyloc,tr);
qzz = _mm512_fmsub_ps(_mm512_set1_ps(3.0f),qj_zzloc,tr);
qxy = _mm512_mul_ps(_mm512_set1_ps(3.0f),qj_xyloc);
qyz = _mm512_mul_ps(_mm512_set1_ps(3.0f),qj_yzloc);
qzx = _mm512_mul_ps(_mm512_set1_ps(3.0f),qj_zxloc);
__fkg_tmp5 = _mm512_mul_ps(_mm512_set1_ps(eps2),tr);
mtr = _mm512_sub_ps(_mm512_set1_ps((PIKG::F32)0.0),__fkg_tmp5);
__fkg_tmp7 = _mm512_mul_ps(qxx,rij.v0);
__fkg_tmp6 = _mm512_fmadd_ps(qxy,rij.v1,__fkg_tmp7);
qr.v0 = _mm512_fmadd_ps(qzx,rij.v2,__fkg_tmp6);
__fkg_tmp9 = _mm512_mul_ps(qyy,rij.v1);
__fkg_tmp8 = _mm512_fmadd_ps(qyz,rij.v2,__fkg_tmp9);
qr.v1 = _mm512_fmadd_ps(qxy,rij.v0,__fkg_tmp8);
__fkg_tmp11 = _mm512_mul_ps(qzz,rij.v2);
__fkg_tmp10 = _mm512_fmadd_ps(qzx,rij.v0,__fkg_tmp11);
qr.v2 = _mm512_fmadd_ps(qyz,rij.v1,__fkg_tmp10);
__fkg_tmp13 = _mm512_fmadd_ps(qr.v0,rij.v0,mtr);
__fkg_tmp12 = _mm512_fmadd_ps(qr.v1,rij.v1,__fkg_tmp13);
rqr = _mm512_fmadd_ps(qr.v2,rij.v2,__fkg_tmp12);
rqr_r4_inv = _mm512_mul_ps(rqr,r4_inv);
meff = _mm512_fmadd_ps(_mm512_set1_ps(0.5f),rqr_r4_inv,mjloc);
__fkg_tmp14 = _mm512_fmadd_ps(_mm512_set1_ps(2.5f),rqr_r4_inv,mjloc);
meff3 = _mm512_mul_ps(__fkg_tmp14,r3_inv);
pf = _mm512_fnmadd_ps(meff,r_inv,pf);
__fkg_tmp15 = _mm512_fnmadd_ps(r5_inv,qr.v0,af.v0);
af.v0 = _mm512_fmadd_ps(meff3,rij.v0,__fkg_tmp15);
__fkg_tmp16 = _mm512_fnmadd_ps(r5_inv,qr.v1,af.v1);
af.v1 = _mm512_fmadd_ps(meff3,rij.v1,__fkg_tmp16);
__fkg_tmp17 = _mm512_fnmadd_ps(r5_inv,qr.v2,af.v2);
af.v2 = _mm512_fmadd_ps(meff3,rij.v2,__fkg_tmp17);
} // loop of j
af.v0 = _mm512_mask_blend_ps(_cvtu32_mask16(0b00000001),__fkg_tmp18.v0,af.v0);
af.v1 = _mm512_mask_blend_ps(_cvtu32_mask16(0b00000001),__fkg_tmp18.v1,af.v1);
af.v2 = _mm512_mask_blend_ps(_cvtu32_mask16(0b00000001),__fkg_tmp18.v2,af.v2);
pf = _mm512_mask_blend_ps(_cvtu32_mask16(0b00000001),__fkg_tmp19,pf);
} // if of j tail loop

((float*)&force[i+0].acc.x)[0] += _mm512_reduce_add_ps(af.v0);

((float*)&force[i+0].acc.y)[0] += _mm512_reduce_add_ps(af.v1);

((float*)&force[i+0].acc.z)[0] += _mm512_reduce_add_ps(af.v2);

((float*)&force[i+0].phi)[0] += _mm512_reduce_add_ps(pf);

} // loop of i
{ // tail loop of reference 
for(;i < ni;++i){
PIKG::F32vec xiloc;

xiloc.x = xiloc_tmp_x[i+0];
xiloc.y = xiloc_tmp_y[i+0];
xiloc.z = xiloc_tmp_z[i+0];
PIKG::F32vec af;

af.x = 0.0f;
af.y = 0.0f;
af.z = 0.0f;
PIKG::F32 pf;

pf = 0.0f;
for(j = 0;j < nj;++j){
PIKG::F32 mjloc;

mjloc = mjloc_tmp[j+0];
PIKG::F32 qj_xxloc;

qj_xxloc = qj_xxloc_tmp[j+0];
PIKG::F32 qj_xyloc;

qj_xyloc = qj_xyloc_tmp[j+0];
PIKG::F32 qj_yyloc;

qj_yyloc = qj_yyloc_tmp[j+0];
PIKG::F32 qj_yzloc;

qj_yzloc = qj_yzloc_tmp[j+0];
PIKG::F32 qj_zxloc;

qj_zxloc = qj_zxloc_tmp[j+0];
PIKG::F32 qj_zzloc;

qj_zzloc = qj_zzloc_tmp[j+0];
PIKG::F32vec xjloc;

xjloc.x = xjloc_tmp_x[j+0];
xjloc.y = xjloc_tmp_y[j+0];
xjloc.z = xjloc_tmp_z[j+0];
PIKG::F32vec rij;

PIKG::F32 __fkg_tmp1;

PIKG::F32 __fkg_tmp0;

PIKG::F32 r2;

PIKG::F32 r_inv;

PIKG::F32 __fkg_tmp2;

PIKG::F32 tmp;

PIKG::F32 __fkg_tmp3;

PIKG::F32 r2_inv;

PIKG::F32 r3_inv;

PIKG::F32 r4_inv;

PIKG::F32 r5_inv;

PIKG::F32 __fkg_tmp4;

PIKG::F32 tr;

PIKG::F32 qxx;

PIKG::F32 qyy;

PIKG::F32 qzz;

PIKG::F32 qxy;

PIKG::F32 qyz;

PIKG::F32 qzx;

PIKG::F32 __fkg_tmp5;

PIKG::F32 mtr;

PIKG::F32 __fkg_tmp7;

PIKG::F32 __fkg_tmp6;

PIKG::F32vec qr;

PIKG::F32 __fkg_tmp9;

PIKG::F32 __fkg_tmp8;

PIKG::F32 __fkg_tmp11;

PIKG::F32 __fkg_tmp10;

PIKG::F32 __fkg_tmp13;

PIKG::F32 __fkg_tmp12;

PIKG::F32 rqr;

PIKG::F32 rqr_r4_inv;

PIKG::F32 meff;

PIKG::F32 __fkg_tmp14;

PIKG::F32 meff3;

PIKG::F32 __fkg_tmp15;

PIKG::F32 __fkg_tmp16;

PIKG::F32 __fkg_tmp17;

rij.x = (xjloc.x-xiloc.x);
rij.y = (xjloc.y-xiloc.y);
rij.z = (xjloc.z-xiloc.z);
__fkg_tmp1 = (rij.x*rij.x+eps2);
__fkg_tmp0 = (rij.y*rij.y+__fkg_tmp1);
r2 = (rij.z*rij.z+__fkg_tmp0);
r_inv = rsqrt(r2);
__fkg_tmp2 = (r_inv*r_inv);
tmp = (3.0f - r2*__fkg_tmp2);
__fkg_tmp3 = (tmp*0.5f);
r_inv = (r_inv*__fkg_tmp3);
r2_inv = (r_inv*r_inv);
r3_inv = (r2_inv*r_inv);
r4_inv = (r2_inv*r2_inv);
r5_inv = (r2_inv*r3_inv);
__fkg_tmp4 = (qj_xxloc+qj_yyloc);
tr = (__fkg_tmp4+qj_zzloc);
qxx = (3.0f*qj_xxloc-tr);
qyy = (3.0f*qj_yyloc-tr);
qzz = (3.0f*qj_zzloc-tr);
qxy = (3.0f*qj_xyloc);
qyz = (3.0f*qj_yzloc);
qzx = (3.0f*qj_zxloc);
__fkg_tmp5 = (eps2*tr);
mtr = -(__fkg_tmp5);
__fkg_tmp7 = (qxx*rij.x);
__fkg_tmp6 = (qxy*rij.y+__fkg_tmp7);
qr.x = (qzx*rij.z+__fkg_tmp6);
__fkg_tmp9 = (qyy*rij.y);
__fkg_tmp8 = (qyz*rij.z+__fkg_tmp9);
qr.y = (qxy*rij.x+__fkg_tmp8);
__fkg_tmp11 = (qzz*rij.z);
__fkg_tmp10 = (qzx*rij.x+__fkg_tmp11);
qr.z = (qyz*rij.y+__fkg_tmp10);
__fkg_tmp13 = (qr.x*rij.x+mtr);
__fkg_tmp12 = (qr.y*rij.y+__fkg_tmp13);
rqr = (qr.z*rij.z+__fkg_tmp12);
rqr_r4_inv = (rqr*r4_inv);
meff = (0.5f*rqr_r4_inv+mjloc);
__fkg_tmp14 = (2.5f*rqr_r4_inv+mjloc);
meff3 = (__fkg_tmp14*r3_inv);
pf = (pf - meff*r_inv);
__fkg_tmp15 = (af.x - r5_inv*qr.x);
af.x = (meff3*rij.x+__fkg_tmp15);
__fkg_tmp16 = (af.y - r5_inv*qr.y);
af.y = (meff3*rij.y+__fkg_tmp16);
__fkg_tmp17 = (af.z - r5_inv*qr.z);
af.z = (meff3*rij.z+__fkg_tmp17);
} // loop of j

force[i+0].acc.x = (force[i+0].acc.x+af.x);
force[i+0].acc.y = (force[i+0].acc.y+af.y);
force[i+0].acc.z = (force[i+0].acc.z+af.z);
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
